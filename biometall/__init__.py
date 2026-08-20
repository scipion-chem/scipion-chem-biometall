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

from pwchem import Plugin as pwchemPlugin
from .constants import *
from pwchem.constants import RDKIT_DIC

_references = ['']


class Plugin(pwchemPlugin):
    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(BIOMETALL_DIC['home'], cls.getEnvName(BIOMETALL_DIC))

    @classmethod
    def defineBinaries(cls, env):
        cls.addBioMetAllPackage(env)

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
            "pip install biometall==1.0 freesasa pandas",
            f"{BIOMETALL_DIC['name']}_installed"
        )
        installer.addPackage(
            env,
            dependencies=['conda', 'pip', 'git'],
            default=default
        )






