# **************************************************************************
# *
# * Authors:     Eduardo Rivas Tortuero
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

from ..scripts.homology_search import search_homologs
from biometall import Plugin
from biometall.constants import METALS_SUPPORTED

_STD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE', 'HSD', 'HSE', 'HSP',
}


class ProtMetalPlacer(Protocol):
    """
    Path 1: homology-based metal identification + MetalKB placement.

    Identifies metal type(s) from an experimentally determined homologue
    in the RCSB PDB, then uses MetalKB (MESPEUS statistical potentials)
    to place each ion type in the query protein geometry.
    """
    _label     = 'place metal ion(s) ? homology + MetalKB (Path 1)'
    _devStatus = BETA

    def __init__(self, **kwargs):
        Protocol.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', params.PointerParam,
                      pointerClass='AtomStruct',
                      label='Apo protein structure',
                      important=True,
                      help='PDB structure of the protein without metal ions.')
        form.addParam('chainId', params.StringParam,
                      default='A',
                      label='Chain ID',
                      help='Protein chain for the homologue search.  Default: A.')

        form.addSection(label='Homologue search')
        form.addParam('identityCutoff', params.FloatParam,
                      default=0.30,
                      label='Min sequence identity',
                      help='Minimum fractional sequence identity (0.30 = 30 %).\n\n'
                           'At ? 30 % identity the 3D fold and metal-binding site '
                           'are almost always conserved (Chothia & Lesk, 1986).\n\n'
                           'If no homologue is found above this threshold, the '
                           'protocol stops with a warning ? no output is produced.')
        form.addParam('evalueCutoff', params.FloatParam,
                      default=1e-5,
                      label='E-value cutoff',
                      help='Maximum acceptable alignment e-value.\n'
                           '1e-5 = 1 in 100,000 chance of a random match.')
        form.addParam('maxHomologs', params.IntParam,
                      default=5,
                      label='Max homologues to retrieve',
                      help='Number of RCSB hits to fetch; the first with a '
                           'confirmed metal is used.')
        form.addParam('metalFilter', params.EnumParam,
                      choices=['ANY'] + METALS_SUPPORTED,
                      default=0,
                      label='Metal to search for',
                      help='ANY: place all metal types found in the homologue '
                           '(multi-metal support).\n'
                           'Specific metal: restrict placement to that type only.')

        form.addSection(label='MetalKB')
        form.addParam('energyThreshold', params.FloatParam,
                      default=-1.7,
                      label='Energy threshold (kcal/mol)',
                      help='MetalKB energy threshold (negative value; more negative '
                           '= stricter).  -1.7 kcal/mol is the standard threshold '
                           'for typical metal sites.  Values closer to 0 '
                           '(e.g. -1.5) also detect weaker or unusual sites.')

    def _insertAllSteps(self):
        self._insertFunctionStep('extractSequenceStep',
                                 self.inputStructure.get().getObjId())
        self._insertFunctionStep('searchHomologsStep')
        self._insertFunctionStep('runMetalKBStep')
        self._insertFunctionStep('createOutputStep')

    # ?? STEP 1 ????????????????????????????????????????????????????????????????

    def extractSequenceStep(self, inputId):
        pdb_path = self.inputStructure.get().getFileName()
        chain_id = (self.chainId.get() or 'A').strip().upper()
        try:
            from Bio.PDB import PDBParser
            from Bio.SeqUtils import seq1
        except ImportError:
            raise RuntimeError('BioPython is required.')

        parser    = PDBParser(QUIET=True)
        structure = parser.get_structure('query', pdb_path)
        model     = structure[0]
        available = [c.id for c in model]
        if chain_id not in available:
            self._log.warning(f'Chain {chain_id} not found; using {available[0]}')
            chain_id = available[0]

        chain    = model[chain_id]
        residues = [r for r in chain
                    if r.get_resname().strip().upper() in _STD_AA]
        if not residues:
            raise RuntimeError(f'No standard amino acids in chain {chain_id}.')

        sequence = ''.join(seq1(r.get_resname().strip()) for r in residues)
        state    = {
            'sequence':        sequence,
            'chain_id':        chain_id,
            'pdb_path':        pdb_path,
            'metals_to_place': {},
            'homolog':         None,
            'identity':        None,
            'placements':      [],
        }
        with open(self._getTmpPath('placer_state.json'), 'w') as f:
            json.dump(state, f)
        self._log.info(f'Sequence: {len(sequence)} aa, chain {chain_id}.')

    # ?? STEP 2 ????????????????????????????????????????????????????????????????

    def searchHomologsStep(self):
        with open(self._getTmpPath('placer_state.json')) as f:
            state = json.load(f)

        choices      = ['ANY'] + METALS_SUPPORTED
        metal_filter = choices[self.metalFilter.get()]
        metals       = (set(METALS_SUPPORTED)
                        if metal_filter == 'ANY' else {metal_filter})

        self._log.info(
            f'Searching RCSB: identity ? {self.identityCutoff.get():.0%}, '
            f'e-value ? {self.evalueCutoff.get():.0e}, metal = {metal_filter}'
        )
        try:
            hits = search_homologs(
                sequence=state['sequence'],
                identity_cutoff=self.identityCutoff.get(),
                evalue_cutoff=self.evalueCutoff.get(),
                max_results=self.maxHomologs.get(),
                metals=metals,
            )
        except RuntimeError as e:
            self._log.warning(f'RCSB error: {e}')
            hits = []

        if not hits:
            self._log.warning(
                f'No homologue with ? {self.identityCutoff.get():.0%} identity '
                f'and a metal ion found in the PDB.\n'
                f'The protein may be a novel metalloprotein with no close '
                f'structural relatives. Manual inspection is recommended.'
            )
            with open(self._getTmpPath('placer_state.json'), 'w') as f:
                json.dump(state, f)
            return

        best         = hits[0]
        metal_counts = best.get('metal_counts', {})

        # Restrict to requested metal(s)
        if metal_filter == 'ANY':
            metals_to_place = metal_counts
        else:
            metals_to_place = {metal_filter: metal_counts.get(metal_filter, 1)}

        self._log.info(
            f'Best homologue: {best["pdb_id"]} '
            f'(identity = {best["identity"]:.0%}, '
            f'metals = {metals_to_place})'
        )
        state.update({
            'metals_to_place': metals_to_place,
            'homolog':         best['pdb_id'],
            'identity':        best['identity'],
        })
        with open(self._getTmpPath('placer_state.json'), 'w') as f:
            json.dump(state, f)

    # ?? STEP 3 ????????????????????????????????????????????????????????????????

    def runMetalKBStep(self):
        with open(self._getTmpPath('placer_state.json')) as f:
            state = json.load(f)

        metals_to_place = state.get('metals_to_place', {})
        if not metals_to_place:
            return

        threshold = float(self.energyThreshold.get())
        metalkb   = self._resolveMetalKBBin()

        if not os.path.isfile(metalkb):
            raise RuntimeError(
                f'MetalKB binary not found: {metalkb}\n'
                'Set the path in the MetalKB section, or install via:\n'
                '  scipion3 installb metalkb'
            )

        input_pdb = self._getTmpPath('input.pdb')
        shutil.copy(state['pdb_path'], input_pdb)

        placements        = []
        all_placed_so_far = []   # positions of all ions placed in previous iterations

        for metal, expected_count in metals_to_place.items():
            self._log.info(
                f'Running MetalKB: {metal}  '
                f'threshold = {threshold} kcal/mol  '
                f'expected = {expected_count}'
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
                    'count_expected':  expected_count,
                    'count_placed':    0,
                    'positions':       [],
                    'coord_residues':  {},
                })
                continue
            except subprocess.TimeoutExpired:
                self._log.warning(f'MetalKB timed out (120 s) for {metal}.')
                placements.append({
                    'metal':          metal,
                    'count_expected': expected_count,
                    'count_placed':   0,
                    'positions':      [],
                    'coord_residues': {},
                })
                continue

            # Rename outputs before next iteration to avoid overwrite
            raw_pdb = self._getTmpPath('out.pdb')
            raw_dat = self._getTmpPath('out.dat')
            metal_pdb = self._getTmpPath(f'out_{metal}.pdb')
            metal_dat = self._getTmpPath(f'out_{metal}.dat')
            if os.path.exists(raw_pdb):
                shutil.move(raw_pdb, metal_pdb)
            if os.path.exists(raw_dat):
                shutil.move(raw_dat, metal_dat)

            all_positions  = _parse_metalkb_positions(metal_pdb, metal)
            coord_residues = _parse_metalkb_coordresidues(metal_dat)


            selected_orig_idx = []
            positions         = []
            for orig_i, pos in enumerate(all_positions, start=1):
                if len(positions) >= expected_count:
                    break
                if any(_dist3(pos, placed) < 1.5 for placed in all_placed_so_far):
                    self._log.info(
                        f'{metal} site {orig_i} skipped ? within 1.5 Å of an '
                        f'already-placed ion.'
                    )
                    continue
                selected_orig_idx.append(orig_i)
                positions.append(pos)

            all_placed_so_far.extend(positions)

            # Keep only coord_residues for the placed sites, renumbered 1..N
            coord_residues_placed = {}
            for rank, orig_i in enumerate(selected_orig_idx, start=1):
                key = str(orig_i)
                if key in coord_residues:
                    coord_residues_placed[str(rank)] = coord_residues[key]

            # Assign each placed position to the nearest protein chain so the
            # viewer can build chain-specific PyMOL selections (chain A and resi ...)
            chains = [_assign_chain(state['pdb_path'], pos) for pos in positions]

            if not positions:
                self._log.warning(
                    f'MetalKB found no site for {metal} '
                    f'above {threshold} kcal/mol.'
                )
            elif len(positions) < expected_count:
                self._log.warning(
                    f'MetalKB: {len(positions)}/{expected_count} '
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
                'count_expected': expected_count,
                'count_placed':   len(positions),
                'positions':      positions,
                'chains':         chains,
                'coord_residues': coord_residues_placed,
            })

        state['placements'] = placements
        with open(self._getTmpPath('placer_state.json'), 'w') as f:
            json.dump(state, f)

    # ?? STEP 4 ????????????????????????????????????????????????????????????????

    def createOutputStep(self):
        with open(self._getTmpPath('placer_state.json')) as f:
            state = json.load(f)

        placements = state.get('placements', [])
        placed     = [p for p in placements if p['count_placed'] > 0]

        if not state.get('metals_to_place'):
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

        out_pdb = self._getExtraPath('pseudo_holo.pdb')
        _write_pseudo_holo(state['pdb_path'], placed, out_pdb)

        summary = {
            'method':      'homology_metalkb',
            'homolog_pdb': state.get('homolog'),
            'identity':    state.get('identity'),
            'placements':  placements,
        }
        with open(self._getExtraPath('placer_result.json'), 'w') as f:
            json.dump(summary, f, indent=2)

        out_struct = emobj.AtomStruct(filename=out_pdb)
        self._defineOutputs(outputStructure=out_struct)
        self._defineRelation(pwobj.RELATION_SOURCE,
                             self.inputStructure, out_struct)

        for p in placed:
            pos_str = '; '.join(
                str([round(v, 3) for v in pos]) for pos in p['positions'])
            self._log.info(
                f'  {p["metal"]} × {p["count_placed"]}: {pos_str}'
            )
        self._log.info(f'Pseudo-holo written: {out_pdb}')

    # ?? Helpers ???????????????????????????????????????????????????????????????

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

def _parse_metalkb_positions(out_pdb, metal):
    """
    Read ALL ions of *metal* from MetalKB's out.pdb.
    MetalKB outputs them sorted by energy (most favourable first).
    Returns list of [x, y, z].
    """
    positions = []
    if not os.path.isfile(out_pdb):
        return positions
    with open(out_pdb) as fh:
        for line in fh:
            if line.startswith('HETATM') and metal in line[17:20]:
                try:
                    positions.append([float(line[30:38]),
                                      float(line[38:46]),
                                      float(line[46:54])])
                except ValueError:
                    pass
    return positions


def _parse_metalkb_coordresidues(out_dat):
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


def _dist3(a, b):
    """Euclidean distance between two [x, y, z] points."""
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) ** 0.5


def _assign_chain(pdb_path, position):
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


def _write_pseudo_holo(query_pdb, placements, out_path):
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