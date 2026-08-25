# **************************************************************************
# *
# * Authors:     Eduardo Rivas Tortuero
# *              Blanca Pueche (blanca.pueche@cnb.csis.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
protocol_metal_placer.py  ?  ProtMetalPlacer
=============================================
Scipion protocol (Path 1): places metal ions by homology identification
followed by MetalKB positioning.

Steps

  1. extractSequenceStep  ? extract amino-acid sequence from the query chain
  2. searchHomologsStep   ? RCSB PDB: find a homologue (?30% id, ?1e-5 evalue)
                            that contains metal(s); read metal types + counts
                            from metadata (NO PDB download or superimposition)
                            = if no homologue found: WARNING, no output
  3. runMetalKBStep       ? run MetalKB once per metal type using MESPEUS
                            statistical potentials; take N best sites per metal
                            (N = count in homologue); warn if fewer found
  4. createOutputStep     ? write pseudo_holo.pdb with ALL placed ions and
                            register output

Outputs

  outputStructure : AtomStruct ? pseudo-holo PDB with all placed ions
  extra/placer_result.json    ? homologue, metals, placements (positions +
                                coordinating residues per metal type)

Multi-metal support
-------------------
  If metalFilter = ANY and the homologue contains multiple metal types
  (e.g. ZN + CU in SOD1), MetalKB is run once per type and all ions are
  written to a single pseudo-holo.  Use metalFilter = <METAL> to restrict
  placement to a single type.
"""

import json
import os
import shutil
import subprocess

from pyworkflow.protocol import Protocol, params
from pyworkflow import BETA
import pyworkflow.object as pwobj
import pwem.objects as emobj
from Bio.PDB import PDBParser, MMCIFParser
from pwem.convert import cifToPdb
from Bio.SeqUtils import seq1

from ..scripts.homology_search import search_homologs
from biometall import Plugin
from biometall.constants import METALS_SUPPORTED
from pwem.protocols import EMProtocol

_STD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE', 'HSD', 'HSE', 'HSP',
}


class ProtMetalPlacer(EMProtocol):
    """
    Predicts metal-binding sites in a protein structure using BioMetAll and
    MetalKB guided by homologous metalloproteins.

    The protocol identifies homologous proteins with experimentally determined
    metal ions in the Protein Data Bank (PDB), infers the most likely metal
    type(s) present in the query protein, and predicts their positions using
    MetalKB statistical potentials.

    Workflow
    --------
    1. Receive a protein structure without metal ions.
    2. Convert the structure to PDB format if necessary.
    3. Extract the amino acid sequence from the selected chain.
    4. Search the Protein Data Bank for homologous metalloproteins satisfying
       the user-defined sequence identity and e-value thresholds.
    5. Determine the metal type(s) and number of ions present in the best
       homolog.
    6. Run MetalKB independently for each predicted metal type.
    7. Select the highest-scoring metal-binding sites according to the expected
       number of ions and discard duplicate placements.
    8. Generate a pseudo-holo structure by adding the predicted metal ions to
       the original protein structure.

    Input
    -----
    - inputStructure:
        Protein structure for which metal ions will be predicted.

        The structure should correspond to the apo form (without metal ions).

    - chainId:
        Protein chain used for sequence extraction and homology search.

    Parameters
    ----------
    Homologue search

    - Min sequence identity:
        Minimum sequence identity required for homologues retrieved from the
        Protein Data Bank.

        The default value is 30%, a commonly accepted threshold for structural
        conservation of protein folds and many metal-binding sites.

    - E-value cutoff:
        Maximum acceptable sequence alignment E-value.

    - Max homologues to retrieve:
        Maximum number of homologous structures inspected before selecting the
        first metalloprotein satisfying the search criteria.

    - Metal:
        Metal type to place.

        *ANY* predicts all metal species identified in the selected homolog,
        while selecting a specific metal restricts the prediction to that metal.

    MetalKB

    - Energy threshold:
        MetalKB energy threshold used to accept predicted metal-binding sites.

        More negative values produce stricter predictions, whereas values closer
        to zero allow weaker candidate sites.

    Output
    ------
    - outputStructure:
        Pseudo-holo protein structure containing the predicted metal ions
        inserted as HETATM records.

    - placer_result.json:
        Summary of the prediction containing:

        - Selected homolog
        - Sequence identity
        - Predicted metal type(s)
        - Expected and placed ions
        - Coordinates of the predicted metal ions
        - Coordinating residues associated with each predicted site

    Summary
    -------
    The protocol reports:

    - Selected homolog and sequence identity.
    - Number of predicted metal-binding sites for each metal type.
    - Expected versus successfully placed metal ions.

    Use Cases
    ---------
    - Restoring missing metal ions in predicted protein structures.
    - Building pseudo-holo models from apo structures.
    - Predicting biologically relevant metal-binding sites through homology.
    - Preparing metalloprotein models for structural analysis or simulation.

    Notes
    -----
    Metal identities are inferred from experimentally determined homologous
    structures deposited in the Protein Data Bank.

    Metal positions are predicted using MetalKB statistical potentials rather
    than copied directly from the homologous structure.

    If no suitable metalloprotein homolog is found, or if MetalKB cannot
    identify energetically favorable binding sites, no output structure is
    generated.
    """
    _label     = 'homology-based metal identification and MetalKB placement'


    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', params.PointerParam,
                      pointerClass='AtomStruct',
                      label='Protein structure: ',
                      important=True,
                      help='PDB structure of the protein without metal ions.')
        form.addParam('chainId', params.StringParam,
                      label='Chain: ',
                      help='Protein chain for the homologue search.')

        group = form.addGroup('Homologue search')
        group.addParam('identityCutoff', params.FloatParam,
                      default=0.30,
                      label='Min sequence identity: ',
                      help='Minimum fractional sequence identity (0.30 = 30 %).\n\n'
                           'At 30 % identity the 3D fold and metal-binding site '
                           'are almost always conserved (Chothia & Lesk, 1986).\n\n'
                           'If no homologue is found above this threshold, the '
                           'protocol stops with a warning and no output is produced.')
        group.addParam('evalueCutoff', params.FloatParam,
                      default=1e-5, expertLevel=params.LEVEL_ADVANCED,
                      label='E-value cutoff: ',
                      help='Maximum acceptable alignment e-value.\n'
                           '1e-5 = 1 in 100,000 chance of a random match.')
        group.addParam('maxHomologs', params.IntParam,
                      default=5,
                      label='Max homologues to retrieve: ',
                      help='Number of RCSB hits to fetch; the first with a '
                           'confirmed metal is used.')
        group.addParam('metalFilter', params.EnumParam,
                      choices=['ANY'] + METALS_SUPPORTED,
                      default=0,
                      label='Metal: ',
                      help='ANY: place all metal types found in the homologue '
                           '(multi-metal support).\n'
                           'Specific metal: restrict placement to that type only.')

        metal = form.addGroup('MetalKB params')
        metal.addParam('energyThreshold', params.FloatParam,
                      default=-1.7,
                      label='Energy threshold: ',
                      help='MetalKB energy threshold (kcal/mol) (negative value; more negative '
                           '= stricter).  -1.7 kcal/mol is the standard threshold '
                           'for typical metal sites.  Values closer to 0 '
                           '(e.g. -1.5) also detect weaker or unusual sites.')

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('extractSequenceStep')
        self._insertFunctionStep('searchHomologsStep')
        self._insertFunctionStep('runMetalKBStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        filePath = self.inputStructure.get().getFileName()
        proteinName = os.path.splitext(
            os.path.basename(filePath)
        )[0]
        ext = os.path.splitext(filePath)[1].lower()
        outputFile = self._getExtraPath(
            f'{proteinName}.pdb'
        )
        if ext == '.pdb':
            shutil.copy(filePath, outputFile)
        elif ext == '.cif':
            cifToPdb(filePath, outputFile)

    def extractSequenceStep(self):
        originalFile = self.inputStructure.get().getFileName()
        proteinName = os.path.splitext(os.path.basename(originalFile))[0]
        structureFile = self._getExtraPath(f'{proteinName}.pdb'
                                            )
        try:
            chainInfo = json.loads(self.chainId.get())
            chainId = chainInfo['chain']
        except (json.JSONDecodeError, TypeError, KeyError):
            chainId = self.chainId.get().strip()

        ext = os.path.splitext(structureFile)[1].lower()

        if ext in ('.cif', '.mmcif'):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure('query', structureFile)
        model = structure[0]

        if chainId not in model:
            raise RuntimeError(
                f'Chain {chainId} not found in input structure.'
            )

        chain = model[chainId]

        residues = [
            r for r in chain
            if r.get_resname().strip().upper() in _STD_AA
        ]

        if not residues:
            raise RuntimeError(
                f'No standard amino acids in chain {chainId}.'
            )

        sequence = ''.join(
            seq1(r.get_resname().strip().upper())
            for r in residues
        )

        state = {
            'sequence': sequence,
            'chainId': chainId,
            'pdb_path': structureFile,
            'metalsToPlace': {},
            'homolog': None,
            'identity': None,
            'placements': [],
        }

        with open(self._getTmpPath('placer_state.json'), 'w') as f:
            json.dump(state, f)

        self._log.info(
            f'Sequence: {len(sequence)} aa, chain {chainId}.'
        )


    def searchHomologsStep(self):
        with open(self._getTmpPath('placer_state.json')) as f:
            state = json.load(f)

        choices      = ['ANY'] + METALS_SUPPORTED
        metalFilter = choices[self.metalFilter.get()]
        metals       = (set(METALS_SUPPORTED)
                        if metalFilter == 'ANY' else {metalFilter})

        self._log.info(
            f'Searching RCSB: identity ? {self.identityCutoff.get():.0%}, '
            f'e-value ? {self.evalueCutoff.get():.0e}, metal = {metalFilter}'
        )
        try:
            hits = search_homologs(
                sequence=state['sequence'],
                identityCutoff=self.identityCutoff.get(),
                evalueCutoff=self.evalueCutoff.get(),
                maxResults=self.maxHomologs.get(),
                metals=metals,
            )
        except RuntimeError as e:
            self._log.warning(f'RCSB error: {e}')
            hits = []

        if not hits:
            self._log.warning(
                f'No homologue with {self.identityCutoff.get():.0%} identity '
                f'and a metal ion found in the PDB.\n'
                f'The protein may be a novel metalloprotein with no close '
                f'structural relatives. Manual inspection is recommended.'
            )
            with open(self._getTmpPath('placer_state.json'), 'w') as f:
                json.dump(state, f)
            return

        best         = hits[0]
        metalCounts = best.get('metalCounts', {})

        # Restrict to requested metal(s)
        if metalFilter == 'ANY':
            metalsToPlace = metalCounts
        else:
            metalsToPlace = {metalFilter: metalCounts.get(metalFilter, 1)}

        self._log.info(
            f'Best homologue: {best["pdb_id"]} '
            f'(identity = {best["identity"]:.0%}, '
            f'metals = {metalsToPlace})'
        )
        state.update({
            'metalsToPlace': metalsToPlace,
            'homolog':         best['pdb_id'],
            'identity':        best['identity'],
        })
        with open(self._getTmpPath('placer_state.json'), 'w') as f:
            json.dump(state, f)


    def runMetalKBStep(self):
        with open(self._getTmpPath('placer_state.json')) as f:
            state = json.load(f)

        metalsToPlace = state.get('metalsToPlace', {})
        if not metalsToPlace:
            return

        threshold = float(self.energyThreshold.get())
        metalkb   = self._resolveMetalKBBin()

        if not os.path.isfile(metalkb):
            raise RuntimeError(
                f'MetalKB binary not found: {metalkb}\n'
            )

        inputPdb = self._getTmpPath('input.pdb')
        shutil.copy(state['pdb_path'], inputPdb)

        placements        = []
        allPlacedSoFar = []   # positions of all ions placed in previous iterations

        for metal, expectedCount in metalsToPlace.items():
            self._log.info(
                f'Running MetalKB: {metal}  '
                f'threshold = {threshold} kcal/mol  '
                f'expected = {expectedCount}'
            )
            try:
                subprocess.run(
                    [metalkb, 'input.pdb', metal, str(threshold)],
                    check=True, capture_output=True, timeout=120,
                    cwd=self._getTmpPath(),
                )
            except subprocess.CalledProcessError as e:
                self._log.warning(
                    f'MetalKB error for {metal}:\n'
                    f'{e.stderr.decode(errors="replace")[:300]}'
                )
                placements.append({
                    'metal':           metal,
                    'count_expected':  expectedCount,
                    'count_placed':    0,
                    'positions':       [],
                    'coordResidues':  {},
                })
                continue
            except subprocess.TimeoutExpired:
                self._log.warning(f'MetalKB timed out (120 s) for {metal}.')
                placements.append({
                    'metal':          metal,
                    'count_expected': expectedCount,
                    'count_placed':   0,
                    'positions':      [],
                    'coordResidues': {},
                })
                continue

            # Rename outputs before next iteration to avoid overwrite
            rawPdb = self._getTmpPath('out.pdb')
            rawDat = self._getTmpPath('out.dat')
            metalPdb = self._getTmpPath(f'out_{metal}.pdb')
            metalDat = self._getTmpPath(f'out_{metal}.dat')
            if os.path.exists(rawPdb):
                shutil.move(rawPdb, metalPdb)
            if os.path.exists(rawDat):
                shutil.move(rawDat, metalDat)

            allPositions  = self._parse_metalkb_positions(metalPdb, metal)
            coordResidues = self._parse_metalkb_coordresidues(metalDat)


            selectedOrigIdx = []
            positions         = []
            for origI, pos in enumerate(allPositions, start=1):
                if len(positions) >= expectedCount:
                    break
                if any(self._dist3(pos, placed) < 1.5 for placed in allPlacedSoFar):
                    self._log.info(
                        f'{metal} site {origI} skipped ? within 1.5 Å of an '
                        f'already-placed ion.'
                    )
                    continue
                selectedOrigIdx.append(origI)
                positions.append(pos)

            allPlacedSoFar.extend(positions)

            # Keep only coordResidues for the placed sites, renumbered 1..N
            coordResiduesPlaced = {}
            for rank, origI in enumerate(selectedOrigIdx, start=1):
                key = str(origI)
                if key in coordResidues:
                    coordResiduesPlaced[str(rank)] = coordResidues[key]

            # Assign each placed position to the nearest protein chain so the
            # viewer can build chain-specific PyMOL selections (chain A and resi ...)
            chains = [self._assign_chain(state['pdb_path'], pos) for pos in positions]

            if not positions:
                self._log.warning(
                    f'MetalKB found no site for {metal} '
                    f'above {threshold} kcal/mol.'
                )
            elif len(positions) < expectedCount:
                self._log.warning(
                    f'MetalKB: {len(positions)}/{expectedCount} '
                    f'{metal} site(s) found ? placing {len(positions)}.'
                )
            else:
                self._log.info(
                    f'{metal}: {len(positions)} site(s) placed ? '
                    + ', '.join(str([round(v, 3) for v in p])
                                for p in positions)
                )

            placements.append({
                'metal':          metal,
                'count_expected': expectedCount,
                'count_placed':   len(positions),
                'positions':      positions,
                'chains':         chains,
                'coordResidues': coordResiduesPlaced,
            })

        state['placements'] = placements
        with open(self._getTmpPath('placer_state.json'), 'w') as f:
            json.dump(state, f)


    def createOutputStep(self):
        with open(self._getTmpPath('placer_state.json')) as f:
            state = json.load(f)

        placements = state.get('placements', [])
        placed     = [p for p in placements if p['count_placed'] > 0]

        if not state.get('metalsToPlace'):
            self._log.warning(
                'No output: no homologue with metal found.\n'
                'Consider lowering the identity cutoff or running '
                'ProtMetalScreener first.'
            )
            return

        if not placed:
            self._log.warning(
                'No output: MetalKB found no sites for any metal.\n'
                f'Try a lower energy threshold (current: '
                f'{self.energyThreshold.get()} kcal/mol).'
            )
            return

        outPdb = self._getExtraPath('pseudo_holo.pdb')
        self._write_pseudo_holo(state['pdb_path'], placed, outPdb)

        summary = {
            'method':      'homology_metalkb',
            'homolog_pdb': state.get('homolog'),
            'identity':    state.get('identity'),
            'placements':  placements,
        }
        with open(self._getExtraPath('placer_result.json'), 'w') as f:
            json.dump(summary, f, indent=2)

        outStruct = emobj.AtomStruct(filename=outPdb)
        self._defineOutputs(outputStructure=outStruct)
        self._defineRelation(pwobj.RELATION_SOURCE,
                             self.inputStructure, outStruct)

        for p in placed:
            posStr = '; '.join(
                str([round(v, 3) for v in pos]) for pos in p['positions'])
            print(
                f'  {p["metal"]} × {p["count_placed"]}: {posStr}'
            )
        print(f'Pseudo-holo written: {outPdb}')


    def _resolveMetalKBBin(self):
        return Plugin.getMetalKBBin()

    def _validate(self):
        errors = []
        if self.inputStructure.get() is None:
            errors.append('An input protein structure is required.')
        if not (0.20 <= float(self.identityCutoff.get()) <= 1.0):
            errors.append('Min sequence identity must be between 0.20 and 1.0.')
        if float(self.evalueCutoff.get()) <= 0:
            errors.append('E-value cutoff must be positive.')
        return errors

    def _summary(self):
        result_json = self._getExtraPath('placer_result.json')
        if not os.path.exists(result_json):
            return ['Protocol not yet executed.']
        with open(result_json) as f:
            r = json.load(f)
        lines = [
            'Homologue: *%s*  (%.0f%% identity)'
            % (r.get('homolog_pdb', '?'), (r.get('identity') or 0) * 100)
        ]
        for p in r.get('placements', []):
            lines.append(
                '*%s*: %d/%d site(s) placed'
                % (p['metal'], p['count_placed'], p['count_expected'])
            )
        return lines

    def _citations(self):
        return ['Berman2000', 'Cock2009', 'Zhao2026', 'Lin2024', 'Chothia1986']

    def _methods(self):
        return [
            'Metal ion placement (Path 1) used homology-based metal '
            'identification followed by MetalKB positioning.  A homologous '
            'PDB entry with ? %.0f %% sequence identity (e-value ? %.0e) '
            'containing experimentally determined metal(s) was found via the '
            'RCSB REST API.  MetalKB placed each ion type using MESPEUS '
            'statistical potentials specific to the identified metal.'
            % (self.identityCutoff.get() * 100, self.evalueCutoff.get())
        ]


# ?? Module helpers ????????????????????????????????????????????????????????????

    def _parse_metalkb_positions(self,outPdb, metal):
        """
        Read ALL ions of *metal* from MetalKB's out.pdb.
        MetalKB outputs them sorted by energy (most favourable first).
        Returns list of [x, y, z].
        """
        positions = []
        if not os.path.isfile(outPdb):
            return positions
        with open(outPdb) as fh:
            for line in fh:
                if line.startswith('HETATM') and metal in line[17:20]:
                    try:
                        positions.append([float(line[30:38]),
                                          float(line[38:46]),
                                          float(line[46:54])])
                    except ValueError:
                        pass
        return positions


    def _parse_metalkb_coordresidues(self,out_dat):
        """
        Parse MetalKB's out.dat into {site_idx: [(resname, resnum), ...]}.
        Format: <site_idx>   <RESNAME>   <resnum>
        Keys are strings (consistent with JSON round-trip).
        Used for the text summary; viewer uses geometric selection instead.
        """
        coord = {}
        if not os.path.isfile(out_dat):
            return coord
        with open(out_dat) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        idx     = str(int(parts[0]))   # string key, consistent with JSON
                        resname = parts[1].upper()
                        resnum  = parts[2]
                        coord.setdefault(idx, []).append((resname, resnum))
                    except ValueError:
                        pass
        return coord


    def _dist3(self,a, b):
        """Euclidean distance between two [x, y, z] points."""
        return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) ** 0.5


    def _assign_chain(self,pdb_path, position):
        """
        Return the chain ID of the CA atom nearest to *position* in *pdb_path*.
        Used to associate a placed metal position with a protein chain so that
        PyMOL coord-residue selections can be restricted to the correct chain.
        """
        min_d2   = float('inf')
        best_ch  = 'A'
        x0, y0, z0 = position
        with open(pdb_path) as fh:
            for line in fh:
                if not line.startswith('ATOM'):
                    continue
                if line[12:16].strip() != 'CA':
                    continue
                try:
                    x  = float(line[30:38])
                    y  = float(line[38:46])
                    z  = float(line[46:54])
                    d2 = (x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2
                    if d2 < min_d2:
                        min_d2  = d2
                        best_ch = line[21].strip() or 'A'
                except ValueError:
                    pass
        return best_ch


    def _write_pseudo_holo(self,query_pdb, placements, out_path):
        """
        Write the query PDB + one HETATM line per placed ion.
        *placements* is a list of dicts with keys 'metal' and 'positions'.
        """
        with open(query_pdb) as fin, open(out_path, 'w') as fout:
            for line in fin:
                if line.startswith(('END', 'CONECT')):
                    continue
                fout.write(line)
            atom_idx = 1
            for placement in placements:
                metal = placement['metal']
                for x, y, z in placement['positions']:
                    fout.write(
                        f"HETATM{atom_idx:5d} {metal:<4s} {metal:<3s} A{atom_idx:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          "
                        f"{metal:>2s}\n"
                    )
                    atom_idx += 1
            fout.write("END\n")