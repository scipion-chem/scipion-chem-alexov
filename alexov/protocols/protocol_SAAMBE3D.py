# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Natalia del Rey
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

"""
Wrapper around the SAAMBE3D method from http://compbio.clemson.edu/saambe_webserver/
"""

import os

from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct
from pwchem.utils import cleanPDB

from alexov import Plugin

class ProtocolSAAMBE3D(EMProtocol):
    """
    This protocol computes the change in free energy at the interface between two proteins
    when there is a mutation in one of the proteins.
    """
    _label = 'SAAMBE3D'
    _devStatus = BETA

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct",
                      label='Atomic structure', allowsNull=False,
                      help='The atomic structure should have the two proteins')

        form.addParam('toMutateList', params.TextParam, width=70,
                      default='', label='List of mutations:',
                      help='The syntax of a mutation is "[Chain] [Position] [aaFrom] [aaTo]", for instance, '
                            'the mutation A 182 C Y, mutates position 182 of chain A that is a (C)ystein to a '
                            't(Y)rosine. For the one-letter aminoacid code see https://en.wikipedia.org/wiki/Amino_acid')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.computeDDG)

    def computeDDG(self):
        fnPDB = self._getExtraPath("atomicStructure.pdb")
        cleanPDB(self.inputAtomStruct.get().getFileName(),fnPDB)
        fnMut = self._getExtraPath("mutations.txt")
        with open(fnMut,"w") as fh:
            fh.write(self.toMutateList.get()+"\n")
        args="-i %s -d 1 -o %s -f %s"%(fnPDB, self._getExtraPath("ddg.txt"), fnMut)
        Plugin.runSAAMBE(self, args=args)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = ['We calculated the ddG for the mutations using the SAAMBE method described in [Pahari2020]']
        return methods
