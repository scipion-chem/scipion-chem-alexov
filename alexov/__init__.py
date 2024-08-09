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

# General imports
import os

# Scipion em imports
import pyworkflow.utils as pwutils
import pwem
from scipion.install.funcs import InstallHelper

# Plugin imports
from .constants import *

__version__ = ALEXOV_VERSION  # plugin version
_logo = "alexov_icon.png"
_references = ['Pahari2020']


class Plugin(pwem.Plugin):
    """
    Definition of class variables. For each package, a variable will be created.
    _<packageNameInLowercase>Home will contain the full path of the package, ending with a folder whose name will be <packageNameFirstLetterLowercase>-<defaultPackageVersion> variable.
    
    Inside that package, for each binary, there will also be another variable.
    _<binaryNameInLowercase>Binary will be a folder inside _<packageNameInLowercase>Home and its name will be <binaryName>.
    """
    _url = "https://github.com/scipion-chem/scipion-chem-alexov"
    _supportedVersions = V1_0

    saambeDefaultVersion = SAAMBE_DEFAULT_VERSION  
    _saambeHome = os.path.join(pwem.Config.EM_ROOT, f'saambe-{saambeDefaultVersion}')
    _saambeBinary = os.path.join(_saambeHome, 'saambe')


    @classmethod
    def _defineVariables(cls):
        """
        Return and write a home and conda enviroment variable in the config file.
        Each package will have a variable called <packageNameInUppercase>_HOME, and another called <packageNameInUppercase>_ENV
        <packageNameInUppercase>_HOME will contain the path to the package installation."
        <packageNameInUppercase>_ENV will contain the name of the conda enviroment for that package."
        """
        cls._defineEmVar(SAAMBE_HOME, cls._saambeHome)
        cls._defineEmVar(SAAMBE_BINARY, cls._saambeBinary)

        cls._defineVar('SAAMBE_ENV', f'saambe-{cls.saambeDefaultVersion}')

    @classmethod
    def defineBinaries(cls, env):
        """
        This function defines the binaries for each protocol.
        """
        cls.addSAAMBE(env)

    @classmethod    
    def addSAAMBE(cls, env):
        """
        This function provides the neccessary commands for installing SAAMBE-3D.
        """
        # Defining protocol variables
        packageName = 'saambe'

        # Instanciating installer
        installer = InstallHelper(packageName, packageHome=cls.getVar(SAAMBE_HOME),packageVersion=cls.saambeDefaultVersion)

        # Installing protocol    
        installer.getCloneCommand('https://github.com/delphi001/SAAMBE-3D', binaryFolderName=packageName)\
            .getCondaEnvCommand(pythonVersion='3.11.5', binaryPath=cls._saambeBinary, requirementsFile=False, 
                                requirementList=['numpy', 'prody', 'xgboost'])\
            .addPackage(env, dependencies=['git', 'conda'])
    
    @classmethod
    def runSAAMBE(cls, protocol, args):
        """ Run saambe command from a given protocol. """
        protocol.runJob("python",os.path.join(cls.getVar(SAAMBE_BINARY),"saambe-3d.py")+" "+args)




    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def getProtocolEnvName(cls, protocolName, repoName=None):
        """
        This function returns the env name for a given protocol and repo.
        """
        return f"{repoName if repoName else protocolName}-{getattr(cls, protocolName + 'DefaultVersion')}"
    
    @classmethod
    def getProtocolActivationCommand(cls, protocolName, repoName=None):
        """
        Returns the conda activation command for the given protocol.
        """
        return f"conda activate {cls.getProtocolEnvName(protocolName, repoName)}"

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch my program. """
        environ = pwutils.Environ(os.environ)
        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')
        
        return neededProgs

