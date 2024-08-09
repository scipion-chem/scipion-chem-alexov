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
import pwem.convert as emconv
from pwchem.utils.utils import cleanPDB

from alexov import Plugin
from alexov.constants import AA_THREE_TO_ONE

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
                            't(Y)rosine. For the one-letter aminoacid code see https://en.wikipedia.org/wiki/Amino_acid.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.computeDDG)

    def computeDDG(self):
        fnPDB = self._getExtraPath("atomicStructure.pdb")
        cleanPDB(self.inputAtomStruct.get().getFileName(),fnPDB)
        fnMut = self._getExtraPath("mutations.txt")
        with open(fnMut,"w") as fh:
            fh.write(self.toMutateList.get().upper()+"\n")
        args="-i %s -d 1 -o %s -f %s"%(fnPDB, self._getExtraPath("ddg.txt"), fnMut)
        Plugin.runSAAMBE(self, args=args)

        os.remove(fnPDB)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
    
        structureHandler = emconv.AtomicStructHandler()
        structureHandler.read(self.inputAtomStruct.get().getFileName())
        structureHandler.getStructure()
        modelsLength, modelsFirstResidue = structureHandler.getModelsChains()
        
        validChains = set()
        chainResidues = {}

        for modelID, chains in modelsFirstResidue.items():
            for chainID, residues in chains.items():
                filtered_residues = [res for res in residues if res[1] != 'HOH']
                validChains.add(chainID)
                if chainID not in chainResidues:
                    chainResidues[chainID] = filtered_residues

        if not self.toMutateList.get().strip():
            errors.append('You have not added any mutation to the list.')
        else:
            for i, line in enumerate(self.toMutateList.get().strip().split('\n')):
                parts = line.upper().split()

                if len(parts) != 4:
                    errors.append(f'The mutation "{line}" does not have the 4 necessary parameters. ' 
                                  'Mutation format must be "[Chain] [Position] [aaFrom] [aaTo]".')
                else:
                    chain, position, aaFrom, aaTo = parts

                    if chain not in validChains:
                        errors.append(f'The chain "{chain}" of the mutation "{line}" is not present in the PDB file. '
                                      f'The PDB file contains the following chains: {", ".join(validChains)}.')
                    
                    elif not position.isdigit():
                        errors.append(f'The position of the mutation "{line}" must be an integer.')
                    
                    elif aaFrom not in AA_THREE_TO_ONE.values():
                        errors.append(f'The wild-type aminoacid of the mutation "{line}" does not ' 
                                      'exist or is not written with its one-letter code.')
                    
                    elif aaTo not in AA_THREE_TO_ONE.values():
                        errors.append(f'The mutant aminoacid of the mutation "{line}" does not ' 
                                       'exist or is not written with its one-letter code.')

                    else:   
                        position = int(position)             
                        residues_dict = {res[0]: res[1] for res in chainResidues[chain]}
                        if position not in residues_dict.keys():
                            first_residue = list(residues_dict.keys())[0]
                            last_residue = list(residues_dict.keys())[-1]
                            errors.append(f'Position "{position}" in chain "{chain}" for mutation "{line}" is out of range. '
                                          f'The chain "{chain}" has positions from {first_residue} to {last_residue}.')
                        
                        elif AA_THREE_TO_ONE[residues_dict[position]] != aaFrom:
                            errors.append(f'The wild-type aminoacid "{aaFrom}" at position "{position}" in chain "{chain}" '
                                          f'for mutation "{line}" does not match the PDB file. The aminoacid at that position '
                                          f'is {residues_dict[position]} ({AA_THREE_TO_ONE[residues_dict[position]]}).')
       
        return errors

    def _summary(self):
        summary = []
        ddgFile = self._getExtraPath('ddg.txt')
        if os.path.exists(ddgFile):
            with open(ddgFile) as f:
              summary.append(f.read())    
        return summary

    def _methods(self):
        methods = []
        methods.append("Prediction of the binding free energy change (ΔΔG) for protein-protein interactions "
                       "due to a point mutation in an aminoacid using the SAAMBE-3D method.")
        return methods
    
    def _citations(self):
        return ['Pahari2020']