#!/usr/bin/python3

import os
from pyworkflow.protocol import Protocol, params
from pyworkflow import BETA
import pwem.objects as emobj
from scipionbiometall import Plugin


class ProtBioMetAll(Protocol):
    _label = 'predict metal binding sites'
    _devStatus = BETA

    def  __init__(self, **kwargs):
        Protocol.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', params.PointerParam,
                      pointerClass='AtomStruct',
                      label='Input protein structure',
                      important=True,
                      help='PDB structure on which BioMetAll will search'
                            'for metal coordination sites.')

        form.addSection(label='Basic options')
        form.addParam('residues', params.StringParam,
                      default = 'HIS,CYS,ASP,GLU',
                      label='Coordinating residues',
                      help='Residue types BioMetAll will consider as potential'
                           'metal coordinators. Comma-separated three letter codes.'
                           'Default: HIS,CYS,ASP,GLU')
        form.addParam('minCoordinators', params.IntParam,
                      default= 3,
                      label='Minimum coordinating residues',
                      help='Minimum number of coordinating residues'
                            'required to consider a site valid.'
                            'Higher values = higher confidence sites.')
        form.addParam('cutoff', params.FloatParam,
                      default= 0.0,
                      label='Cutoff fraction',
                      help='Filters results, keeping only the top'
                            'fraction by probe count.\n'
                            '0.0 = shows all results. \n'
                           '0.5 = shows only the top 50% of sites.')

        form.addSection(label='Backbone')
        form.addParam('useBackbone', params.BooleanParam,
                      default=False,
                      label='Include backbone oxygens',
                      help= 'By default BioMetAll considers sidechain atoms.'
                            'Enable this option to include backbone carbonyl oxygens'
                            'as potential coordination donors.')
        form.addParam('backboneResidues', params.StringParam,
                      default='ALL',
                      condition='useBackbone',
                      label='Backbone residues',
                      help='Residues whose backbone oxygens will be included.\n'
                           'ALL = every residue.\n'
                           'Or a comma-separated list of residues: HIS,CYS,ASP')

        form.addSection(label='Search zone')
        form.addParam('useSearchZone', params.BooleanParam,
                      default=False,
                      label='Restrict search to a sphere',
                      help= 'By default BioMetAll builds a probe grid covering'
                            'the entire protein. \n'
                            'Enable this to restrict search to a sphere'
                            'defined by a centre and radius.')
        form.addParam('center', params.StringParam,
                      default='',
                      condition='useSearchZone',
                      label='Sphere center',
                      help= 'XYZ coordinates of the sphere centre.\n'
                            'Format: x,y,z (no spaces).\n'
                            'Example: 84.98,42.82,16.04')
        form.addParam('radius', params.FloatParam,
                      default=10.0,
                      condition='useSearchZone',
                      label='Sphere radius (Angstrom)',
                      help= 'Radius in Angstroms of the search sphere.\n'
                            'Default: 10.0 A. Aprox 3-4 residues from the centre.')
        form.addParam('grid', params.FloatParam,
                      default=1.0,
                      label= 'Grid spacing (Angstrom)',
                      help='Distance in Angstroms between probe grid points.\n'
                           '1.0 A = standard resolution (default).\n'
                           '0.5 A = higher detail, longer computation.\n')

        form.addSection(label= 'Motif')
        form.addParam('useMotif', params.BooleanParam,
                      default=False,
                      label='Search for a specific motif',
                      help='Instead of finding any coordination site, search'
                           'for specific pattern of coordinating residues.')
        form.addParam('motif', params.StringParam,
                      default='',
                      condition='useMotif',
                      label= 'Motif [AA,AA...]',
                      help= 'Coordination motif to search for.\n'
                            'Format: three letter codes comma-separated.\n'
                            'Examples:\n'
                            'HIS,HIS,ASP\n'
                            'HIS/CYS,HIS,ASP position 1 can be HIS or CYS.')

    def _insertAllSteps(self):
        self._insertFunctionStep(
            'convertInputStep',
            self.inputStructure.get().getObjId()
        )
        self._insertFunctionStep('runBioMetAllStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self, inputId):
        import shutil
        pdbPath = self.inputStructure.get().getFileName()
        shutil.copy(pdbPath, self._getTmpPath('input.pdb'))

    def runBioMetAllStep(self):
        args = os.path.abspath(self._getTmpPath('input.pdb')) + ' --pdb'
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

        Plugin.runBioMetAll(self, args, cwd=self._getTmpPath())

        Plugin.runBioMetAll(self, args, cwd=self._getTmpPath())

        import glob
        import shutil
        probesFiles = glob.glob(os.path.join(self._getTmpPath(), 'probes_*.pdb'))
        if probesFiles:
            shutil.move(probesFiles[0], self._getTmpPath('probes_input.pdb'))

    def createOutputStep(self):
        probesFile = self._getTmpPath('probes_input.pdb')

        if not os.path.exists(probesFile):
            raise FileNotFoundError(
                'BioMetAll did not generate the expected probes file: %s\n'
                'Check the log for error details.' % probesFile
            )

        outputProbes = self._getExtraPath('probes.pdb')
        import shutil
        shutil.copy(probesFile, outputProbes)

        probesAtomStruct = emobj.AtomStruct(filename=outputProbes)
        self._defineOutputs(outputProbes=probesAtomStruct)

    def _validate(self):
        errors = []

        if self.inputStructure.get() is None:
            errors.append('Please select an input protein structure.')

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


