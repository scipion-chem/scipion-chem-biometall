# **************************************************************************
# *
# * Authors:   Blanca Pueche (blanca.pueche@cnb.csis.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import os
import glob
import shutil
from pwem.protocols import EMProtocol
from pyworkflow.protocol import Protocol, params
from pyworkflow import BETA
import pwem.objects as emobj
from biometall import Plugin
from biometall.constants import BIOMETALL_DIC

from pwem.convert import cifToPdb

from pwchem.objects import SetOfStructROIs, StructROI


class ProtBioMetAll(EMProtocol):
    _label = 'predict metal binding sites'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', params.PointerParam,
                      pointerClass='AtomStruct',
                      label='Input structure: ',allowsNull=False,
                      help='PDB structure on which BioMetAll will search'
                            'for metal coordination sites.')

        basic = form.addGroup('Basic options')
        basic.addParam('residues', params.StringParam,
                      default = 'HIS,CYS,ASP,GLU',
                      label='Coordinating residues: ',
                      help='Residue types BioMetAll will consider as potential'
                           'metal coordinators. Comma-separated three letter codes.'
                           'Default: HIS,CYS,ASP,GLU')
        basic.addParam('minCoordinators', params.IntParam,
                      default= 3,
                      label='Minimum coordinating residues: ',
                      help='Minimum number of coordinating residues'
                            'required to consider a site valid.'
                            'Higher values = higher confidence sites.')
        basic.addParam('cutoff', params.FloatParam,
                      default= 0.0,
                      label='Cutoff fraction: ',
                      help='Filters results, keeping only the top'
                            'fraction by probe count.\n'
                            '0.0 = shows all results. \n'
                           '0.5 = shows only the top 50% of sites.')

        bb = form.addGroup('Backbone')
        bb.addParam('useBackbone', params.BooleanParam,
                      default=False,
                      label='Include backbone oxygens: ',
                      help= 'By default BioMetAll considers sidechain atoms.'
                            'Enable this option to include backbone carbonyl oxygens'
                            'as potential coordination donors.')
        bb.addParam('backboneResidues', params.StringParam,
                      default='ALL', expertLevel=params.LEVEL_ADVANCED,
                      condition='useBackbone',
                      label='Backbone residues: ',
                      help='Residues whose backbone oxygens will be included.\n'
                           'ALL = every residue.\n'
                           'Or a comma-separated list of residues: HIS,CYS,ASP')

        sZone = form.addGroup('Search zone', expertLevel=params.LEVEL_ADVANCED)
        sZone.addParam('useSearchZone', params.BooleanParam,
                      default=False, expertLevel=params.LEVEL_ADVANCED,
                      label='Restrict search to a sphere: ',
                      help= 'By default BioMetAll builds a probe grid covering'
                            'the entire protein. \n'
                            'Enable this to restrict search to a sphere'
                            'defined by a centre and radius.')
        sZone.addParam('center', params.StringParam,
                      default='', expertLevel=params.LEVEL_ADVANCED,
                      condition='useSearchZone',
                      label='Sphere center: ',
                      help= 'XYZ coordinates of the sphere centre.\n'
                            'Format: x,y,z (no spaces).\n'
                            'Example: 84.98,42.82,16.04')
        sZone.addParam('radius', params.FloatParam,
                      default=10.0, expertLevel=params.LEVEL_ADVANCED,
                      condition='useSearchZone',
                      label='Sphere radius (Angstrom): ',
                      help= 'Radius in Angstroms of the search sphere.\n'
                            'Default: 10.0 A. Aprox 3-4 residues from the centre.')
        sZone.addParam('grid', params.FloatParam,
                      default=1.0, expertLevel=params.LEVEL_ADVANCED,
                      label= 'Grid spacing (Angstrom): ',
                      help='Distance in Angstroms between probe grid points.\n'
                           '1.0 A = standard resolution (default).\n'
                           '0.5 A = higher detail, longer computation.\n')

        motif = form.addGroup('Motif')
        motif.addParam('useMotif', params.BooleanParam,
                      default=False,
                      label='Search for a specific motif: ',
                      help='Instead of finding any coordination site, search'
                           'for specific pattern of coordinating residues.')
        motif.addParam('motif', params.StringParam,
                      default='',
                      condition='useMotif',
                      label= 'Motif [AA,AA...]: ',
                      help= 'Coordination motif to search for.\n'
                            'Format: three letter codes comma-separated.\n'
                            'Examples:\n'
                            'HIS,HIS,ASP\n'
                            'HIS/CYS,HIS,ASP position 1 can be HIS or CYS.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.runBioMetAllStep)
        self._insertFunctionStep(self.createOutputStep)

    def convertInputStep(self):
        filePath = self.inputStructure.get().getFileName()
        ext = os.path.splitext(filePath)[1]
        if ext == '.pdb':
            shutil.copy(filePath, self._getExtraPath('biometall.pdb'))
        elif ext == '.cif':
            inpFile = os.path.abspath(self._getExtraPath('biometall.pdb'))
            cifToPdb(filePath, inpFile)

    def runBioMetAllStep(self):
        args = os.path.abspath(self._getExtraPath('biometall.pdb')) + ' --pdb'
        residueList = self.residues.get().strip().split(',')
        if residueList:
            args += ' --residues [%s]' % ','.join(residueList)
        args += ' --min_coordinators %d' % self.minCoordinators.get()

        if self.cutoff.get() > 0.0:
            args += ' --cutoff %.2f' % self.cutoff.get()

        if self.useBackbone.get():
            backbone = self.backboneResidues.get().strip()
            if backbone:
                args += ' --backbone %s' % ','.join(backbone.split())

        if self.useSearchZone.get():
            center = self.center.get().strip()
            if center:
                args += ' --center [%s]' % center
            args += ' --radius %.2f' % self.radius.get()
        args += ' --grid %.2f' % self.grid.get()

        if self.useMotif.get():
            motif = self.motif.get().strip()
            if motif:
                args += ' --motif [%s]' % motif

        Plugin.runCondaCommand(
            self,
            args=args,
            condaDic=BIOMETALL_DIC,
            program="biometall",
            cwd=os.path.abspath(Plugin.getVar(BIOMETALL_DIC['home']))
        )

        probesFiles = glob.glob(os.path.join(self._getExtraPath(), 'probes_*.pdb'))
        if probesFiles:
            shutil.move(probesFiles[0], self._getPath('probes_biometall.pdb'))
        self.splitBiometallSols()

    def createOutputStep(self):
        probesDir = os.path.join(self._getPath(), "probes")
        structRois = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
        proteinFile = self.inputStructure.get().getFileName()

        for pdbFile in sorted(glob.glob(os.path.join(probesDir, "*.pdb"))):
            structRoi = StructROI()
            structRoi.setFileName(pdbFile)
            structRoi.setProteinFile(proteinFile)
            vol = structRoi.getPocketVolume()
            structRoi.setVolume(vol)
            structRois.append(structRoi)

        self._defineOutputs(outputStructROIs=structRois)

    def _validate(self):
        errors = []
        if self.minCoordinators.get() < 2:
            errors.append('Minimum coordinators must be at least 2.')

        if not (0.0 <= self.cutoff.get() < 1.0):
            errors.append('Cutoff must be between 0.0 (inclusive) and 1.0 (exclusive).')

        if self.grid.get() <= 0.0:
            errors.append('Grid spacing must be greater than 0.')

        if self.useSearchZone.get():
            center = self.center.get().strip()
            if not center:
                errors.append('If search zone is enabled, centre coordinates are required.')
            elif len(center.split(',')) != 3:
                errors.append('Centre must be in x,y,z format (e.g. 84.98,42.82,16.04).')

        if self.useMotif.get() and not self.motif.get().strip():
            errors.append('If motif search is enabled, a motif must be specified.')

        return errors

    def _summary(self):
        summary = []
        if self.inputStructure.get() is not None:
            summary.append('Structure: *%s*' % self.inputStructure.get().getFileName())
        summary.append('Residues: %s' % self.residues.get())
        summary.append('Min coordinators: %d' % self.minCoordinators.get())
        if self.cutoff.get() > 0.0:
            summary.append('Cutoff: top %.0f%%' % (self.cutoff.get() * 100))
        if self.useBackbone.get():
            summary.append('Backbone: %s' % self.backboneResidues.get())
        if self.useSearchZone.get():
            summary.append('Zone: centre %s, radius %.1f A'
                           % (self.center.get(), self.radius.get()))
        if self.useMotif.get():
            summary.append('Motif: %s' % self.motif.get())
        return summary

    def _citations(self):
        return ['SanchezAparicio2021']


    def _methods(self):
        msg = (
                'Metal coordination sites were predicted with BioMetAll '
                '[SanchezAparicio2021] using a minimum of %d coordinating residues '
                '(%s), a grid spacing of %.1f Å.'
                % (self.minCoordinators.get(),
                   self.residues.get(),
                   self.grid.get())
        )

        if self.cutoff.get() > 0.0:
            msg += (
                    ' Results were filtered keeping only the top %.0f%% of sites '
                    'by probe density (cutoff=%.2f).'
                    % (self.cutoff.get() * 100, self.cutoff.get())
            )

        if self.useBackbone.get():
            msg += (
                    ' Backbone carbonyl oxygens of %s residues were also considered '
                    'as potential coordination donors.'
                    % self.backboneResidues.get()
            )

        if self.useSearchZone.get():
            msg += (
                    ' The search was restricted to a sphere of radius %.1f Å '
                    'centred at coordinates %s.'
                    % (self.radius.get(), self.center.get())
            )

        if self.useMotif.get():
            msg += (
                    ' The search was restricted to the coordination motif %s.'
                    % self.motif.get()
            )

        return [msg]

    def splitBiometallSols(self):
        probesDir = os.path.join(self._getPath(), "probes")
        os.makedirs(probesDir, exist_ok=True)
        solutions = {}
        inputPdb = self._getPath('probes_biometall.pdb')

        with open(inputPdb, 'r') as f:
            for line in f:
                if not line.startswith(('ATOM', 'HETATM')):
                    continue
                solution_id = int(line[22:26].strip())
                solutions.setdefault(solution_id, []).append(line)

        for solution_id, lines in solutions.items():
            output_file = os.path.join(
                self._getPath('probes'),
                f'biometall_solution_{solution_id}.pdb'
            )
            with open(output_file, 'w') as f:
                f.writelines(lines)
