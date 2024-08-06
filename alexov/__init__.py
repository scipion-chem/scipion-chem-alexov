# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Natl. Center of Biotechnology CSIC
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
import pyworkflow.utils as pwutils
import pwem

from alexov.constants import *

__version__ = "0.1"  # plugin version
_logo = "icon.png"
_references = []


class Plugin(pwem.Plugin):
    _url = "https://github.com/scipion-chem/scipion-chem-alexov"
    _supportedVersions = [V1]  # binary version

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(SAAMBE_BINARY, "saambe")
        cls._defineEmVar(SAAMBE_HOME, f"saambe-{V1}")

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch my program. """
        environ = pwutils.Environ(os.environ)
        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. """
        neededProgs = []

        return neededProgs

    @classmethod
    def runSAAMBE(cls, protocol, args):
        """ Run saambe command from a given protocol. """
        protocol.runJob("python",os.path.join(cls.getVar(SAAMBE_HOME),"saambe-3d.py")+" "+args)

    @classmethod
    def defineBinaries(cls, env):
    	pass
        #installCmds = [("make -j 4", "")]  # replace the target "" with e.g. "bin/myprogram"
        #env.addPackage('alexov', version=V1,
        #               tar='void.tgz',
        #               commands=installCmds,
        #               neededProgs=cls.getDependencies(),
        #               default=True)
