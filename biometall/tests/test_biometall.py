# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import setupTestProject, DataSet, BaseTest

# Scipion chem imports
from pwchem.protocols import  ProtChemPrepareReceptor
from pwchem.utils import assertHandle
from pwem.protocols import ProtImportPdb
from ..protocols import ProtBioMetAll, ProtMetalPlacer


class TestBioMetAll(BaseTest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)

        cls._runImportPDB()
        cls._runPrepareTarget()
        cls._runBioMetAll()

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0, pdbId='1ca2')
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runPrepareTarget(cls):
        protPrepRec = cls.newProtocol(
            ProtChemPrepareReceptor,
            inputAtomStruct=cls.protImportPDB.outputPdb)

        cls.proj.launchProtocol(protPrepRec, wait=True)
        cls.protPrepRec = protPrepRec

    @classmethod
    def _runBioMetAll(cls):
        protBiometall = cls.newProtocol(
            ProtBioMetAll,
            inputStructure=cls.protPrepRec.outputStructure,
            cutoff=0.5
        )

        cls.proj.launchProtocol(protBiometall, wait=True)
        cls.protBiometall = protBiometall

    def test(self):
        self._waitOutput(self.protBiometall, 'outputStructROIs', sleepTime=10)
        assertHandle(self.assertIsNotNone, getattr(self.protBiometall, 'outputStructROIs', None))

class TestMetalPlacer(BaseTest):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)

        cls._runImportPDB()
        cls._runPrepareTarget()
        cls._runPlacer()

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0,
            pdbId='1ca2'
        )
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runPrepareTarget(cls):
        protPrepRec = cls.newProtocol(
            ProtChemPrepareReceptor,
            inputAtomStruct=cls.protImportPDB.outputPdb
        )
        cls.proj.launchProtocol(protPrepRec, wait=True)
        cls.protPrepRec = protPrepRec

    @classmethod
    def _runPlacer(cls):
        protPlacer = cls.newProtocol(
            ProtMetalPlacer,
            inputStructure=cls.protPrepRec.outputStructure,
            chainId='{"model": 0, "chain": "A", "residues": 256}',
            metalFilter=1
        )

        cls.proj.launchProtocol(protPlacer, wait=True)
        cls.protPlacer = protPlacer

    def test(self):
        self._waitOutput(self.protPlacer,'outputStructure', sleepTime=10)
        assertHandle(self.assertIsNotNone,getattr(self.protPlacer, 'outputStructure', None))