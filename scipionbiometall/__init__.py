#!/usr/bin/python3

import os
import pwem
import pyworkflow.utils as pwutils
from .constants import *

_logo = "logo.png"
_references = ['SanchezAparicio2021']

class Plugin(pwem.Plugin):
    _homeVar = BIOMETALL_HOME
    _supportedVersions = [V1_0]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(BIOMETALL_HOME, 'biometall-1.0')

    @classmethod
    def getEnviron(cls):
        environ = pwutils.Environ(os.environ.copy())
        environ.update({
            'PATH': '/home/vboxuser/miniconda3/envs/scipion3/bin'
        }, position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def defineBinaries(cls, env):
        pass

    @classmethod
    def runBioMetAll (cls, protocol, args, cwd=None):
        protocol.runJob('biometall', args, env=cls.getEnviron(), cwd=cwd)




