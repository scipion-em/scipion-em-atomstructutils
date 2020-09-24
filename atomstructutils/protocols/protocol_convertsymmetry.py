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
from pwem.protocols import EMProtocol

from pwem.objects import AtomStruct
from pwem.convert.symmetry import  Icosahedron
from pyworkflow.protocol.params import (EnumParam,
                                        IntParam,
                                        MultiPointerParam,
                                        PointerParam,
                                        StringParam,
                                        LEVEL_ADVANCED)


from pwem.constants import (SYM_I222, SYM_I222r, SYM_In25, SYM_In25r,
                                     SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r,
                                     SCIPION_SYM_NAME)
import numpy as np
from pwem.convert.atom_struct import AtomicStructHandler, fromCIFTommCIF

LOCAL_SYM_NAME = {}
LOCAL_SYM_NAME[SYM_I222] = 'I1'
LOCAL_SYM_NAME[SYM_I222r] = 'I2'
LOCAL_SYM_NAME[SYM_In25] = 'I3'
LOCAL_SYM_NAME[SYM_In25r] = 'I4'
LOCAL_SYM_NAME[SYM_I2n3] = 'I5'
LOCAL_SYM_NAME[SYM_I2n3r] = 'I6'
LOCAL_SYM_NAME[SYM_I2n5] = 'I7'
LOCAL_SYM_NAME[SYM_I2n5r] = 'I8'

class ProtAtomStrucConvertSymmetry(EMProtocol):
    _label = 'convert_sym'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('pdbFileToBeRefined', PointerParam, pointerClass="AtomStruct",
                       label='Atomic Structure:', allowsNull=True,
                       important=True,
                       help="Input the atomic structure to be modified.")

        form.addParam('originSymmetryGroup', EnumParam,
                       choices=[LOCAL_SYM_NAME[SYM_I222] +
                                " (" + SCIPION_SYM_NAME[SYM_I222] + ")",
                                LOCAL_SYM_NAME[SYM_I222r] +
                                " (" + SCIPION_SYM_NAME[SYM_I222r] + ")",
                                LOCAL_SYM_NAME[SYM_In25] +
                                " (" + SCIPION_SYM_NAME[SYM_In25] + ")",
                                LOCAL_SYM_NAME[SYM_In25r] +
                                " (" + SCIPION_SYM_NAME[SYM_In25r] + ")",
                                LOCAL_SYM_NAME[SYM_I2n3] +
                                " (" + SCIPION_SYM_NAME[SYM_I2n3] + ")",
                                LOCAL_SYM_NAME[SYM_I2n3r] +
                                " (" + SCIPION_SYM_NAME[SYM_I2n3r] + ")",
                                LOCAL_SYM_NAME[SYM_I2n5] +
                                " (" + SCIPION_SYM_NAME[SYM_I2n5] + ")",
                                LOCAL_SYM_NAME[SYM_I2n5r] +
                                " (" + SCIPION_SYM_NAME[SYM_I2n5r] + ")"
                                ],
                       default=SYM_I222r - SYM_I222,
                       label="Symmetry",
                       help="Select the current symmetry of your atomic structure./n"
                            "See https://github.com/I2PC/xmipp-portal/wiki/Symmetry"
                            "Symmetry for a description of the symmetry groups "
                       )

        form.addParam('targetSymmetryGroup', EnumParam,
                  choices=[LOCAL_SYM_NAME[SYM_I222] +
                           " (" + SCIPION_SYM_NAME[SYM_I222] + ")",
                           LOCAL_SYM_NAME[SYM_I222r] +
                           " (" + SCIPION_SYM_NAME[SYM_I222r] + ")",
                           LOCAL_SYM_NAME[SYM_In25] +
                           " (" + SCIPION_SYM_NAME[SYM_In25] + ")",
                           LOCAL_SYM_NAME[SYM_In25r] +
                           " (" + SCIPION_SYM_NAME[SYM_In25r] + ")",
                           LOCAL_SYM_NAME[SYM_I2n3] +
                           " (" + SCIPION_SYM_NAME[SYM_I2n3] + ")",
                           LOCAL_SYM_NAME[SYM_I2n3r] +
                           " (" + SCIPION_SYM_NAME[SYM_I2n3r] + ")",
                           LOCAL_SYM_NAME[SYM_I2n5] +
                           " (" + SCIPION_SYM_NAME[SYM_I2n5] + ")",
                           LOCAL_SYM_NAME[SYM_I2n5r] +
                           " (" + SCIPION_SYM_NAME[SYM_I2n5r] + ")"
                           ],
                  default=SYM_I222 - SYM_I222,
                  label="Symmetry",
                  help="Select the desired symmetry of your atomic structure.\n"
                       "See https://github.com/I2PC/xmipp-portal/wiki/Symmetry"
                       "Symmetry for a description of the symmetry groups "
                  )

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        origSym = self.originSymmetryGroup.get() + SYM_I222
        targetSym = self.targetSymmetryGroup.get() + SYM_I222
        self._insertFunctionStep('applyTransformationStep',origSym, targetSym)

    def applyTransformationStep(self, origSym, targetSym):
        matrix = self.getTransformationMatrix(origSym, targetSym)

        originAtomStructFn = self.pdbFileToBeRefined.get().getFileName()
        targetAtomStructFn = self._getExtraPath("outPutAtomStruct.cif")
        self.rotateAtomStruct(inAtomStructFn=originAtomStructFn,
                              outAtomStructFn=targetAtomStructFn,
                              matrix=matrix)
        self.createOutputStep(targetSym, targetAtomStructFn)

    def getTransformationMatrix(self, origin, target):
        "compute matrix that converts between two icosahedal simmetries"
        ico = Icosahedron(orientation = SCIPION_SYM_NAME[origin][1:])
        matrix = ico.coordinateSystemTransform(SCIPION_SYM_NAME[origin][1:],
                                               SCIPION_SYM_NAME[target][1:])
        return np.array(matrix)

    def rotateAtomStruct(self, inAtomStructFn, outAtomStructFn, matrix):
        "apply rotation matrix to input atomic structure"

        atSH = AtomicStructHandler(inAtomStructFn)
        atSH.transform(matrix)
        atSH.write(outAtomStructFn)

    def createOutputStep(self, targetSym, targetAtomStructFn):
        """ save new atomic structure"""
        pdb = AtomStruct()
        pdb.setFileName(targetAtomStructFn)
        # MM: to get appropriate cif files to be visualize with Chimera
        # Transform the output cif file in mmcif
        log = self._log
        fromCIFTommCIF(pdb.getFileName(), pdb.getFileName(), log)
        self._defineOutputs(rotatedAtomStruct=pdb)
        self._defineSourceRelation(self.pdbFileToBeRefined, pdb)

    # --------------------------- UTILS functions ------------------

    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        return summary

#    def _citations(self):
#        return ['Cock2009']