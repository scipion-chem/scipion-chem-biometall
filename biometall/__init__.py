# **************************************************************************
# *
# * Authors:  Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

from scipion.install.funcs import InstallHelper
import os
import glob

from pwchem import Plugin as pwchemPlugin
from .constants import *
from pwchem.constants import RDKIT_DIC

_references = ['']

_BASE_DIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)



class Plugin(pwchemPlugin):
    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(BIOMETALL_DIC['home'], cls.getEnvName(BIOMETALL_DIC))
        cls._defineEmVar(METALKB_DIC['home'], cls.getEnvName(METALKB_DIC))

    @classmethod
    def defineBinaries(cls, env):
        cls.addBioMetAllPackage(env)
        cls.addMetalKBPackage(env)

    @classmethod
    def addBioMetAllPackage(cls, env, default=True):
        installer = InstallHelper(
            BIOMETALL_DIC['name'],
            packageHome=cls.getVar(BIOMETALL_DIC['home']),
            packageVersion=BIOMETALL_DIC['version']
        )
        installer.getCondaEnvCommand(
            BIOMETALL_DIC['name'],
            binaryVersion=BIOMETALL_DIC['version'],
            pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(BIOMETALL_DIC)} && "
            "pip install biometall==1.0 freesasa pandas joblib scikit-learn",
            f"{BIOMETALL_DIC['name']}_installed"
        )
        installer.addPackage(
            env,
            dependencies=['conda', 'pip', 'git'],
            default=default
        )

    @classmethod
    def addMetalKBPackage(cls, env, default=True):
        installer = InstallHelper(
            METALKB_DIC['name'],
            packageHome=cls.getVar(METALKB_DIC['home']),
            packageVersion=METALKB_DIC['version'],
        )
        installer.addCommand(
            "git clone --depth=1 "
            "https://github.com/huang-laboratory/MetalKB.git MetalKB_repo && "
            "cp MetalKB_repo/MetalKB ./MetalKB_binary && "
            "chmod +x ./MetalKB_binary",
            f"{METALKB_DIC['name']}_installed",
        )
        installer.addPackage(
            env,
            dependencies=['conda', 'git'],
            default=default,
        )

    @classmethod
    def getBioMetAllEnvName(cls):
        """Return the conda environment name used for BioMetAll."""
        return cls.getEnvName(BIOMETALL_DIC)

    @classmethod
    def getBaseDir(cls):
        """Root of metalplacer_entrega ? ml_path0 is importable from here."""
        return os.environ.get('METALPLACER_BASE_DIR', _BASE_DIR)

    @classmethod
    def getML0Dir(cls):
        """Directory containing ml_path0 source code (config.yaml lives here)."""
        return os.environ.get(
            'METALPLACER_ML0_DIR',
            os.path.join(cls.getBaseDir(), 'biometall/ml_path0'),
        )

    @classmethod
    def getModelsDir(cls):
        """Directory containing trained model .pkl files and thresholds JSON."""
        return os.environ.get(
            'METALPLACER_MODELS_DIR',
            os.path.join(cls.getML0Dir(), 'models'),
        )

    @classmethod
    def getMetalKBBin(cls):
        """Return the full path to the MetalKB_binary executable."""
        return os.path.join(cls.getVar(METALKB_DIC['home']), 'MetalKB_binary')

    @classmethod
    def getConfigPath(cls):
        """Absolute path to ml_path0/config.yaml."""
        return os.path.join(cls.getML0Dir(), 'config.yaml')

    @classmethod
    def getBinaryModelPath(cls):
        """Calibrated binary screener model.

        Returns the canonical ``model_path0_calibrated.pkl``.
        Falls back to any ``model_path0*_calibrated.pkl`` (timestamped retrains),
        then to any ``*_calibrated.pkl`` as last resort.
        """
        models_dir = cls.getModelsDir()
        # Canonical name (current)
        canonical = os.path.join(models_dir, 'model_path0_calibrated.pkl')
        if os.path.isfile(canonical):
            return canonical
        matches = sorted(glob.glob(os.path.join(models_dir, 'model_path0*_calibrated.pkl')))
        if matches:
            return matches[-1]
        matches = sorted(glob.glob(os.path.join(models_dir, '*_calibrated.pkl')))
        return matches[-1] if matches else None

    @classmethod
    def getDefaultTau(cls):
        """Read tau_protein from thresholds.json; fallback to 0.5."""
        models_dir = cls.getModelsDir()
        # Canonical name (current)
        canonical = os.path.join(models_dir, 'thresholds.json')
        if os.path.isfile(canonical):
            try:
                with open(canonical) as fh:
                    return float(json.load(fh).get('tau_protein', 0.5))
            except Exception:
                pass
        matches = sorted(glob.glob(os.path.join(models_dir, 'thresholds_*.json')))
        for path in reversed(matches):
            try:
                with open(path) as fh:
                    return float(json.load(fh).get('tau_protein', 0.5))
            except Exception:
                continue
        return 0.5


from pyworkflow.plugin import Domain

Domain.registerPlugin(__name__)



