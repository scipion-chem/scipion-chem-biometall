#!/usr/bin/python3

from scipion.install.funcs import InstallHelper
from .constants import *
from pwchem import Plugin as pwchemPlugin

_logo = "logo.png"
_references = ['SanchezAparicio2021']

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


