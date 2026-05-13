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
    """
    Provides a collection of utilities for manipulating and reorganizing
    atomic structure files in PDB or mmCIF format. The protocol is intended
    to simplify common structural editing tasks that frequently arise during
    cryo-EM and structural biology workflows, including combining structures,
    isolating chains, renumbering residues, renaming chains, and separating
    complete models into individual chain files.

    AI Generated:

    Atomic Structure Operator (ProtAtomStrucOperate) - User Manual
        Overview

        The Atomic Structure Operator protocol is designed to facilitate
        practical editing and preparation of macromolecular atomic models.
        In structural biology projects, researchers often need to manipulate
        specific chains, merge structural fragments, or standardize residue
        numbering before downstream refinement, fitting, validation, or
        visualization. This protocol provides a convenient environment for
        carrying out these structural organization tasks while preserving
        compatibility with cryo-EM and molecular modeling workflows.

        The protocol operates on atomic structures represented as PDB or
        mmCIF files and produces new atomic models suitable for continued
        analysis. Its role is not to refine coordinates or alter molecular
        geometry, but rather to reorganize structural information in ways
        that support interpretation, comparison, or integration of models.

        Biological Context and Typical Applications

        In many cryo-EM studies, atomic structures originate from multiple
        sources such as homology modeling, AlphaFold predictions, crystal
        structures, or previously refined cryo-EM models. Before these
        structures can be interpreted together, users frequently need to
        combine chains into a unified assembly, isolate regions of interest,
        or adapt chain identifiers and residue numbering conventions.

        This protocol becomes especially useful in integrative modeling
        workflows, multi-subunit complexes, comparative analyses, and
        iterative refinement pipelines. For example, users may extract a
        flexible domain for independent analysis, merge several subunits
        into a composite model, or renumber residues to match canonical
        sequence annotations used in publications or databases.

        Chain Extraction

        One of the most common operations is extracting a specific chain
        from a larger structure. This is particularly important when working
        with oligomeric complexes, ribosomes, membrane assemblies, or
        heterogeneous particles where only a subset of chains is relevant
        for a given biological question.

        The protocol allows extraction of complete chains or selected
        residue ranges within a chain. This capability is valuable when
        focusing on catalytic domains, flexible regions, ligand-binding
        interfaces, or experimentally resolved fragments. Limiting the
        extracted region to biologically meaningful residues can simplify
        visualization and improve downstream processing efficiency.

        When extracting residue intervals, users should ensure that the
        selected range corresponds to meaningful structural regions.
        Arbitrary truncation may generate incomplete domains that are
        difficult to interpret biologically.

        Extraction of All Chains

        In some workflows, it is useful to automatically separate an entire
        assembly into independent chain files. This protocol supports the
        extraction of all chains individually, generating separate outputs
        for each chain present in the structure.

        This operation is particularly valuable when analyzing symmetry-
        related subunits, comparing homologous chains independently, or
        preparing chain-specific fitting and refinement tasks. It can also
        facilitate automated pipelines where each chain must be processed
        separately.

        Merging Atomic Structures

        The protocol also supports combining several atomic structures into
        a single integrated model. This functionality is especially useful
        for assembling multicomponent complexes from independently derived
        structures or predicted domains.

        In practical cryo-EM modeling, users often fit multiple subunits
        individually into a density map and later combine them into a final
        composite assembly. The protocol simplifies this process by creating
        a unified structure suitable for visualization, validation, and
        deposition workflows.

        Biological care is important when merging structures originating
        from different coordinate systems or conformational states.
        Structures should already be approximately aligned before merging,
        otherwise the resulting assembly may not represent a meaningful
        biological configuration.

        Chain Renaming and Residue Renumbering

        Structural datasets obtained from different software packages or
        databases often use inconsistent chain identifiers or residue
        numbering schemes. The protocol provides utilities to standardize
        these annotations.

        Renaming chains is particularly useful when preparing structures
        for molecular dynamics, comparative analyses, or deposition
        pipelines where unique and biologically interpretable chain labels
        are required. Consistent naming also simplifies communication among
        collaborators and improves readability during visualization.

        Residue renumbering allows users to shift residue indices to match
        canonical protein sequences, experimental constructs, or published
        numbering systems. This operation is especially important when
        comparing structures from different organisms, truncated constructs,
        or engineered variants.

        Visualization and Workflow Integration

        The protocol is designed to integrate naturally with molecular
        visualization environments. Generated structures can be immediately
        inspected visually, facilitating rapid confirmation that the
        desired chains or residue ranges were processed correctly.

        Because the outputs remain in standard structural biology formats,
        they can be directly used in downstream refinement, docking,
        validation, flexible fitting, or comparative structural analysis
        protocols. This makes the protocol a convenient utility component
        within larger cryo-EM processing pipelines.

        Practical Recommendations

        Before performing chain extraction or merging operations, users
        should verify chain identifiers carefully, especially in large
        assemblies containing repeated or symmetry-related chains.
        Incorrect chain selection is one of the most common sources of
        confusion in structural preparation workflows.

        When merging structures, it is advisable to confirm that all models
        are already expressed in compatible coordinate systems. Visual
        inspection after merging is strongly recommended to detect possible
        overlaps or misplaced domains.

        For residue renumbering, maintaining correspondence with known
        biological sequence annotations improves reproducibility and
        simplifies interpretation during publication or deposition.

        Final Perspective

        Structural manipulation and organization are essential components
        of modern cryo-EM and integrative structural biology workflows.
        Although these operations do not directly refine structural quality,
        they strongly influence the clarity, usability, and biological
        interpretability of atomic models. Careful organization of chains,
        residue numbering, and structural assemblies enables more reliable
        downstream refinement, visualization, and biological analysis.
    """
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
        f.write("view\n")
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
