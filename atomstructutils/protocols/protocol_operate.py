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
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct
from pwem.viewers import Chimera
from pyworkflow.protocol.params import (EnumParam,
                                        IntParam,
                                        MultiPointerParam,
                                        PointerParam,
                                        StringParam)
from pwem.convert.atom_struct import AtomicStructHandler, fromCIFTommCIF
import os


class ProtAtomStrucOperate(EMProtocol):
    """Utilities for handling PDB/mmcif atomic structure files.
Current plugin utilities: (A) extract a chain from an atom structure (pdb/cif file),
(B) perform union of several atomic structures"""
    operationsDict = {0: 'addChain',
                      1: 'extractChain',
                      2: 'reNumberChain',
                      3: 'reNameChain',
                      4: 'extractAllChains'}
    operationsDictInv = {value:key for key, value in operationsDict.items()}
    # operationsDictInv = {value:key for key, value in list(operationsDict.items())}
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
                      condition="Operation == %d or Operation == %d" %
                                (self.operationsDictInv['extractChain'],
                                self.operationsDictInv['reNameChain']),
                      label="Chain ", allowsNull=True,
                      help="Select a particular chain of the atomic "
                           "structure.")
        form.addParam('start', IntParam,
                      condition="Operation == %d" %
                                (self.operationsDictInv['extractChain']),
                      label="Start at residue #", allowsNull=True,
                      default=-1,
                      help="Extract Chain starting at this number of residue. "
                           "-1 = first residue")
        form.addParam('offset', IntParam,
                      condition="Operation == %d" %
                                (self.operationsDictInv['reNumberChain']),
                      label="Offset residue by", allowsNull=True,
                      default=-1,
                      help="renumber residues by adding offset value")
        form.addParam('end', IntParam,
                      condition="Operation == %d" % self.operationsDictInv['extractChain'],
                      label="End at residue #", allowsNull=True,
                      default=-1,
                      help="Extract Chain ending at this number of residue."
                           " -1 -> last residue")
        form.addParam('chainName', StringParam,
                      condition="Operation == %d " %
                                (self.operationsDictInv['reNameChain']),
                      label="New Chain name", allowsNull=True,
                      default=-1,
                      help="Give Chain this new name")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        dim = 150.
        sampling = 1.

        bildFileName = os.path.abspath(self._getExtraPath(
            "axis_output.bild"))
        Chimera.createCoordinateAxisFile(dim,
                                 bildFileName=bildFileName,
                                 sampling=sampling)
        fnCmd = self._getExtraPath("chimera_output.cxc")
        f = open(fnCmd, 'w')
        f.write("open %s\n" % bildFileName)
        f.write("cofr 0,0,0\n")  # set center of coordinates
        f.close()
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
        elif self.Operation == self.operationsDictInv['extractAllChains']:
            self._insertFunctionStep('extractAllChainsStep',
                                     self.pdbFileToBeRefined.get().getFileName())
        elif self.Operation == self.operationsDictInv['reNumberChain']:
            self._insertFunctionStep('reNumberChainStep',
                                     self.pdbFileToBeRefined.get().getFileName())
        elif self.Operation == self.operationsDictInv['reNameChain']:
            self._insertFunctionStep('reNameChainStep',
                                     self.pdbFileToBeRefined.get().getFileName()
                                     )
        else:
            raise Exception("ERROR: Invalid operation *%s* I quit" % self.Operation)

    def reNumberChainStep(self, structFileName):
        import json
        outFileName = self._getExtraPath("atomStruct_reNumberedChain.cif")
        aStruct1 = AtomicStructHandler(structFileName)
        chainIdDict = json.loads(self.inputStructureChain.get())
        aStruct1.renumberChain(chainID=chainIdDict['chain'],
                              offset=self.offset.get(),
                              modelID=chainIdDict['model'],
                              filename=outFileName)
        #aStruct1.write(outFileName)
        self.createOutputStep(outFileName)

    def reNameChainStep(self, structFileName):
        import json
        outFileName = self._getExtraPath("atomStruct_reNamedChain.cif")
        aStruct1 = AtomicStructHandler(structFileName)
        chainIdDict = json.loads(self.inputStructureChain.get())
        aStruct1.renameChain(chainID=chainIdDict['chain'],
                              newChainName=self.chainName.get(),
                              modelID=chainIdDict['model'],
                              filename=outFileName)
        #aStruct1.write(outFileName)
        self.createOutputStep(outFileName)

    def addChainStep(self, structFileName, listStructFileName):

        outFileName = self._getExtraPath("atomStruct_addChain.cif")
        aStruct1 = AtomicStructHandler(structFileName)
        print("Adding to Atomic Struct {}".format(structFileName))
        for fileName in listStructFileName:
            print("AddingStruct {}".format(fileName))
            sys.stdout.flush()
            aStruct1.addStruct(fileName, outFileName)
        #aStruct1.write(outFileName)
        self.createOutputStep(outFileName, twoRelations=True)

    def extractChainStep(self, structFileName):

        import json
        outFileName = self._getExtraPath("atomStruct_extractChain.cif")
        aStruct1 = AtomicStructHandler(structFileName)
        chainIdDict = json.loads(self.inputStructureChain.get())
        end = self.end.get()
        if end == -1:
            end = sys.maxsize

        aStruct1.extractChain(chainID=chainIdDict['chain'],
                              start=self.start.get(),
                              end=end,
                              modelID=chainIdDict['model'],
                              filename=outFileName)
        self.createOutputStep(outFileName)

    def extractAllChainsStep(self, structFileName):
        import json
        outFileName = self._getExtraPath("atomStruct_extractChain_%s.cif")
        aStruct1 = AtomicStructHandler(structFileName)
        listOfChains, _ = aStruct1.getModelsChains()
        for model, chainDic in listOfChains.items():
            for chainID, lenResidues in chainDic.items():
                chainIdDict = json.loads('{"model": %d, "chain": "%s", "residues": %d}' % (model, str(chainID), lenResidues))
                chainIDStr=chainIdDict['chain']
                aStruct1.extractChain(modelID=chainIdDict['model'], chainID=chainIDStr,
                                      start=-1, end=sys.maxsize,
                                      filename=outFileName%chainIDStr)
                self.createOutputStep(outFileName%chainIDStr,suffix=chainIDStr)

    def createOutputStep(self, outFileName, twoRelations=False, suffix=''):
        outFileName = os.path.abspath(outFileName)
        fnCmd = self._getExtraPath("chimera_output.cxc")
        f = open(fnCmd, 'a+')
        f.write("open %s\n" % outFileName)
        f.close()

        pdb = AtomStruct()
        pdb.setFileName(outFileName)
        # MM: to get appropriate cif files to be visualize with Chimera
        # Transform the output cif file in mmcif
        log = self._log
        fromCIFTommCIF(outFileName, outFileName, log)

        if suffix=="":
            self._defineOutputs(outputPdb=pdb)
        else:
            outputDict = {'outputPdb_chain%s'%suffix: pdb}
            self._defineOutputs(**outputDict)
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
