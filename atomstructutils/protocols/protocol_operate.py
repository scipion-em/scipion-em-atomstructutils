# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import sys
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em import AtomStruct
from pyworkflow.protocol.params import (EnumParam,
                                        IntParam,
                                        MultiPointerParam,
                                        PointerParam,
                                        StringParam)

class ProtAtomStrucOperate(EMProtocol):
    """Utilities for handling PDB/mmcif atomic structure files.
Current plugin utilities: (A) extract a chain from an atom structure (pdb/cif file),
(B) perform union of several atomic structures"""
    operationsDict = {0: 'addChain', 1: 'extractChain'}
    operationsDictInv = {value:key for key, value in operationsDict.items()}
    _label = 'operator'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('pdbFileToBeRefined', PointerParam, pointerClass="AtomStruct",
                       label='Atomic Structure 1:', allowsNull=True,
                       important=True,
                       help="Input the reference atomic structure.")
        form.addParam('Operation', EnumParam,
                         choices=[val for key, val in sorted(self.operationsDict.items())],
                         label="Operation:", default=0,
                         help="Select operation to be performed")
        form.addParam('InputAtomStruct2', MultiPointerParam, pointerClass="AtomStruct",
                       label='Atomic Structure 2:', allowsNull=True,
                       condition="Operation == %d" % self.operationsDictInv['addChain'],
                       important=True,
                       help="Input the atomic structures to be added.")
        form.addParam('inputStructureChain', StringParam,
                      condition="Operation == %d" % self.operationsDictInv['extractChain'],
                      label="Chain ", allowsNull=True,
                      help="Select a particular chain of the atomic "
                           "structure.")
        form.addParam('start', IntParam,
                      condition="Operation == %d" % self.operationsDictInv['extractChain'],
                      label="Start at residue #", allowsNull=True,
                      default=-1,
                      help="Extract Chain starting at this number of residue. "
                           "-1 = first residue")
        form.addParam('end', IntParam,
                      condition="Operation == %d" % self.operationsDictInv['extractChain'],
                      label="End at residue #", allowsNull=True,
                      default=-1,
                      help="Extract Chain ending at this number of residue."
                           " -1 -> last residue")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        if self.Operation == self.operationsDictInv['addChain']:
            listStructFileName = []
            for aStruct in self.InputAtomStruct2:
                listStructFileName.append(aStruct.get().getFileName())
            self._insertFunctionStep('addChainStep',
                                     self.pdbFileToBeRefined.get().getFileName(),
                                     listStructFileName
                                     )
        elif self.Operation == self.operationsDictInv['extractChain']:
            self._insertFunctionStep('extractChainStep',
                                     self.pdbFileToBeRefined.get().getFileName())
        else:
            raise Exception("ERROR: Invalid operation *%s* I quit" % self.Operation)

    def addChainStep(self, structFileName, listStructFileName):
        # Horrible hack to release this plugin before scipion next version.
        # TODO: remove when possible
        from pyworkflow import LAST_VERSION, VERSION_2_0
        if LAST_VERSION == VERSION_2_0:
            from pyworkflow.utils import importFromPlugin
            AtomicStructHandler = importFromPlugin('chimera.atom_struct', 'AtomicStructHandler')
        else:
            from pyworkflow.em.convert.atom_struct import AtomicStructHandler

        outFileName = self._getExtraPath("atomStruct_addChain.cif")
        aStruct1 = AtomicStructHandler(structFileName)
        print "Adding to Atomic Struct {}".format(structFileName)
        for fileName in listStructFileName:
            print "AddingStruct {}".format(fileName)
            sys.stdout.flush()
            aStruct1.addStruct(fileName, outFileName)
        #aStruct1.write(outFileName)
        self.createOutputStep(outFileName, twoRelations=True)

    def extractChainStep(self, structFileName):
        # Horrible hack to release this plugin before scipion next version.
        # TODO: remove when possible
        from pyworkflow import LAST_VERSION, VERSION_2_0
        if LAST_VERSION == VERSION_2_0:
            from pyworkflow.utils import importFromPlugin
            AtomicStructHandler = importFromPlugin('chimera.atom_struct', 'AtomicStructHandler')
        else:
            from pyworkflow.em.convert.atom_struct import AtomicStructHandler

        import json
        outFileName = self._getExtraPath("atomStruct_extractChain.cif")
        aStruct1 = AtomicStructHandler(structFileName)
        chainIdDict = json.loads(self.inputStructureChain.get())
        end = self.end.get()
        if end == -1:
            end = sys.maxint

        aStruct1.extractChain(chainID=chainIdDict['chain'],
                              start=self.start.get(),
                              end=end,
                              modelID=chainIdDict['model'],
                              filename=outFileName)
        #aStruct1.write(outFileName)
        self.createOutputStep(outFileName)

    def createOutputStep(self, outFileName, twoRelations=False):
        pdb = AtomStruct()
        pdb.setFileName(outFileName)
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.pdbFileToBeRefined, pdb)
        if twoRelations:
            self._defineSourceRelation(self.InputAtomStruct2, pdb)

    # --------------------------- UTILS functions ------------------

    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Cock2009']
