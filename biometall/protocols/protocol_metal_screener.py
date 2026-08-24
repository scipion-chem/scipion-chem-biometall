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
protocol_metal_screener.py  ?  ProtMetalScreener
=================================================
Scipion protocol (Path 0): structural metalloprotein binary screener.

Determines whether an input protein structure is likely to bind a metal ion.
Uses BioMetAll to detect candidate coordination sites and a calibrated Random
Forest classifier trained on 62 3D structural features (3 concentric shells +
geometry + physicochemistry + BioMetAll probe density).

Outputs
  isMetal        : Boolean  ? True if predicted metalloprotein
  pMetal         : Float    ? protein-level probability (0?1)
  nSitesScored   : Integer  ? BioMetAll coordination-site candidates scored
  outputStructure: AtomStruct ? input structure passed through (only if isMetal)
"""

import json
import os, shutil
import subprocess

from pyworkflow.protocol import Protocol, params
from pyworkflow import BETA
import pyworkflow.object as pwobj
from pwem.convert import cifToPdb

from biometall import Plugin
from biometall.constants import BIOMETALL_DIC
from pwem.protocols import EMProtocol
from pwchem.objects import SetOfStructROIs, StructROI


class ProtMetalScreener(EMProtocol):

    _label = 'metal-binding sites screening'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam(
            'inputSites',
            params.PointerParam,
            pointerClass='SetOfStructROIs',
            label='BioMetAll sites: ',
            important=True,
            help=(
                'SetOfStructROIs produced by ProtBioMetAll. '
                'The corresponding protein structure will be recovered '
                'from the BioMetAll output.'
            )
        )

        form.addParam(
            'tauProtein',
            params.FloatParam,
            default=0.60,
            label='Decision threshold: ',
            help=(
                'Probability threshold used to classify the protein '
                'as a metalloprotein.'
            )
        )

    # ------------------------------------------------------------------
    # Steps
    # ------------------------------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runScreenerStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        sites = self.inputSites.get()
        firstSite = sites.getFirstItem()
        filePath = firstSite.getProteinFile()
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

    def runScreenerStep(self):
        sites = self.inputSites.get()
        firstSite = sites.getFirstItem()
        proteinFile = firstSite.getProteinFile()

        proteinFile = os.path.abspath(proteinFile)

        probesDir = self._getExtraPath('biometall_sites')
        os.makedirs(
            probesDir,
            exist_ok=True
        )
        for i, site in enumerate(sites):
            source = site.getFileName()
            if not source:
                continue
            destination = os.path.join(
                probesDir,
                'site_%04d.pdb' % i
            )
            shutil.copy2(
                source,
                destination
            )

        probeFiles = [f for f in os.listdir(probesDir) if f.endswith('.pdb')]

        binaryModel = Plugin.getBinaryModelPath()
        if binaryModel is None:
            raise RuntimeError(
                'Binary model not found in: %s\n'
                'Expected: binary_model_*_calibrated.pkl'
                % Plugin.getModelsDir()
            )

        scriptPath = os.path.join(
            Plugin.getBaseDir(),
            'biometall',
            'ml_path0',
            'predict_structural.py'
        )

        jsonOut = self._getPath(
            'screener_result.json'
        )

        args = [
            scriptPath,
            '--pdb', os.path.abspath(self._getExtraPath('*.pdb')),
            '--probes-dir', os.path.abspath(probesDir),
            '--binary-model', binaryModel,
            '--tau', str(self.tauProtein.get()),
            '--output', os.path.abspath(jsonOut)
        ]

        Plugin.runCondaCommand(
            self,
            args=' '.join(args),
            condaDic=BIOMETALL_DIC,
            program='python',
            cwd=os.path.abspath(self._getPath()
            )
        )

        if not os.path.exists(jsonOut):
            raise RuntimeError(
                'Metal Screener did not produce the expected output: %s'
                % jsonOut
            )

    # ------------------------------------------------------------------
    # Outputs
    # ------------------------------------------------------------------
    def createOutputStep(self):
        jsonPath = self._getPath('screener_result.json')
        with open(jsonPath) as fh:
            res = json.load(fh)

        topSites = res.get('top_sites', [])
        inputRois = self.inputSites.get()
        outRois = inputRois.createCopy(self._getPath(), copyInfo=True)

        for roi, site in zip(inputRois, topSites):
            structRoi = roi.clone()
            structRoi.pMetal = pwobj.Float(
                float(site['p_metal'])
            )
            outRois.append(structRoi)

        self._defineOutputs(outputStructROIs=outRois)


    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    def _summary(self):
        summary = []
        jsonPath = self._getPath(
            'screener_result.json'
        )
        with open(jsonPath) as fh:
            res = json.load(fh)
        isMetal = res.get(
            'is_metalloprotein',
            False
        )
        pMetal = res.get(
            'p_metal',
            0.0
        )
        verdict = (
            '*Metalloprotein*'
            if isMetal
            else '*Non-metalloprotein*'
        )
        summary.append(
            'Result: %s' % verdict
        )
        summary.append(
            'pMetal: %.4f' % pMetal
        )
        return summary

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    def _validate(self):
        errors = []
        return errors

    # ------------------------------------------------------------------
    # Citations
    # ------------------------------------------------------------------
    def _citations(self):
        return [
            'SanchezAparicio2021',
            'Pedregosa2011',
            'Putignano2018',
            'Cock2009',
        ]

    # ------------------------------------------------------------------
    # Methods
    # ------------------------------------------------------------------
    def _methods(self):
        return [
            (
                'Candidate metal-binding sites were generated using '
                'BioMetAll. The resulting coordination sites were '
                'subsequently evaluated using a calibrated Random Forest '
                'classifier based on 62 structural features including '
                'local amino-acid composition, coordination geometry, '
                'physicochemical environment, and BioMetAll probe density. '
                'Site probabilities were aggregated to a protein-level '
                'score using the 95th percentile (p95) strategy. '
                'The protein was classified as a metalloprotein when the '
                'resulting probability exceeded %.2f.'
            )
            % self.tauProtein.get()
        ]