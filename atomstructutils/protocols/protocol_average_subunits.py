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
import os
from pwem.convert import Ccp4Header

from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pwem.viewers import Chimera

from pwem.objects import Float, String
from pyworkflow.protocol.params import (IntParam,
                                        MultiPointerParam,
                                        PointerParam, FloatParam)
import numpy as np
from pwem.objects import Volume

from pwem.convert.atom_struct import AtomicStructHandler, fromCIFTommCIF

class ProtAverageSubunits(EMProtocol):
    """Average densities related with the given atomic structures.
    For example hexons in a unit cell. Alpha carbon are used
    to compute the transformation matrix"""
    _label = 'average_sub_unit'
    _program = ""
    DEBUG = True
    DEBUGX = not DEBUG

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.listTransformationMatrices = []
        self.listRMS = []
        self.listObjID = []

    def outputFileName(self):
        return self._getExtraPath("outPut.mrc")

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      important=True,
                      help="Volume to process")
        form.addParam('atomStructReference', PointerParam, pointerClass="AtomStruct",
                      label='Atomic Structure 1:', allowsNull=True,
                      important=True,
                      help="Input the reference atomic structure.")
        form.addParam('otherAtomStructs', MultiPointerParam,
                      pointerClass="AtomStruct", allowsNull=True,
                      label='Atomic structures',
                      help="Other atom Structs")
        form.addParam('rangeStart', IntParam,
                      label="First Residue", allowsNull=True,
                      default=-1,
                      help="First residue to be taken into account"
                           " -1 -> first residue")
        form.addParam('rangeEnd', IntParam,
                      label="Last Residue", allowsNull=True,
                      default=-1,
                      help="Last residue to be taken into account"
                           " -1 -> last residue")
        form.addParam('expand', FloatParam,
                      label="expand BB", allowsNull=True,
                      default=3,
                      help="Expand bounding box by this"
                           " number of A.")

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutput')

    def computeTransformationMatrices(self, origin, sampling):
        """ compute transformation matrix between reference
        atom structure (self.atomStructReference) and the other
        atom structures (self.otherAtomStructs)
         read atomic structures
         identify common AA with first structure
         compute transformation matrix in PDB system
         function fills lists
            * listTransformationMatrices
            * listRMS (errors)
            * listObjID (atom structure names)
         transform matrix to unitcell system of equations
         shifts in px
        """
        # reference atom struct
        atomStructReferenceFn = \
            self.atomStructReference.get().getFileName()
        self.aStructReferenceHa = AtomicStructHandler(
            atomStructReferenceFn)

        # each one of the atom structs different from the reference
        for aStruct in self.otherAtomStructs:
            atomStructFn = aStruct.get().getFileName()
            matrix, rms = self.aStructReferenceHa.getTransformMatrix(
                atomStructFn, self.rangeStart.get(), self.rangeEnd.get())
            # getTransformMatrix
            # [[ 9.95590478e-01 -3.44402543e-03 -9.37429380e-02 -5.47417e+01]
            #  [-2.80233170e-03  9.97787836e-01 -6.64197386e-02  2.35492e+01]
            #  [ 9.37643145e-02  6.63895581e-02  9.93378417e-01  7.36214e+00]
            #  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.00000e+00]]

            # transform matrix in atom struct (p) reference, shifts A
            # to unit cell (p') reference, shifts px
            convertAtomicModel2Unit = \
                self.convertAtomicModel2UnitCellMatrix(sampling, origin)
            # origin [-318.23999023  -53.04000092  212.16000366    1.]

            convertUnitCell2Atomic = \
                self.convertUnitCell2AtomicModelMatrix(sampling, origin)

            # @ -> matrix_ multiplication
            matrixUnitCell = \
                convertAtomicModel2Unit @ matrix @ convertUnitCell2Atomic
            #  matrixUnitCell
            #  [[ 9.95590478e-01 -3.44402543e-03 -9.37429380e-02  5.3709e+01]
            #   [-2.80233170e-03  9.97787836e-01 -6.64197386e-02 -7.6961e+00]
            #   [ 9.37643145e-02  6.63895581e-02  9.93378417e-01  2.0149e+01]
            #   [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.0000e+00]]

            # scipion uses the inverse matrix
            matrixUnitCellInv = np.linalg.inv(matrixUnitCell)
            # matrixUnitCellInv
            # [[9.95590478e-01 -2.80233170e-03  9.37643145e-02 -5.53830e+01]
            # [-3.44402543e-03  9.97787836e-01  6.63895581e-02  6.52637e+00]
            # [-9.37429380e-02 -6.64197386e-02  9.93378417e-01 -1.54925e+01]
            # [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.00000e+00]]

            self.listTransformationMatrices.append(matrixUnitCellInv)
            self.listRMS.append(Float(rms))  # root mean square error
            self.listObjID.append(String(aStruct.getNameId())) # filenames

    def convertUnitCell2AtomicModelMatrix(self, sampling, origin):
        """ Compute matrix that transform points in unit cell (px)
            to atom struct (A) system"""
        return np.linalg.inv(
            self.convertAtomicModel2UnitCellMatrix(sampling, origin))

    def convertAtomicModel2UnitCellMatrix(self, sampling, origin):
        """ Compute matrix that transform points in atomic structs (A)
            to unit cell (px) system. p' = (p - origin ) /sampling"""
        shiftScale = np.eye(4, dtype=float)
        shiftScale[:, 3] = - origin
        shiftScale = shiftScale / sampling
        shiftScale[3, 3] = 1.
        return shiftScale

    def computeAveragedSubUnit(self):
        """Average the pixel values of the pixels related by symmetry"""
        # 1) get sampling and origin
        volume = self.inputVolume.get()
        origin = np.ones(4)
        origin[:3] = volume.getShiftsFromOrigin()
        sampling = volume.getSamplingRate()

        # 2) get transformation matrices that relate otherAtomStructs
        # with reference atomStructReference, matrix unit cell (p')
        # system of coordinates
        self.computeTransformationMatrices(origin, sampling)

        # 3) get matrix with data (voxel values)
        # careful here, index (0,0,) is top left corner
        # relate to PDB coordinates (p) as:
        # are p' = (p - origin)/sampling
        ih = ImageHandler()
        img = ih.createImage()
        img.read(volume.getLocation())
        matrix = img.getData()  # get matrix with voxels
        # matrix.shape (220, 274, 455) (z,y,x)
        # scipion reports  455 274 220 (x,y,z)

        # 4) compute atomic struct bounding box in reference
        # atom structure in A (p)
        boundingBoxA = np.ones((2,4))
        boundingBoxA[:,:3] = self.aStructReferenceHa.\
            getBoundingBox(expand=self.expand.get())
        self.boundingBoxA = boundingBoxA  # self.bounding needed in viewer

        # compute bounding box in unit cell reference and pixels (p')
        boundingBoxPxOrigin = np.ones((2,4))
        convert = self.convertAtomicModel2UnitCellMatrix(sampling, origin)
        for i in range(2):
            boundingBoxPxOrigin[i,:] = convert @ boundingBoxA[i,:]
            # convert bounding box to closest integers
            boundingBoxPxOrigin = boundingBoxPxOrigin.astype(int)
        # bounding box A [[-176.99899292  -48.58000183  322.22198486    1.]
        #  [ -84.15100098   47.80899811  436.43399048    1.        ]]
        # boundingBoxPxOrigin [[103   3  80   1]
        #  [172  74 164   1]]

        # 5) loop through points and interpolated
        # grid points to average
        l = boundingBoxPxOrigin  # bounding box in pixels (int) around origin
        # values for tests (px - float)
        # 00032:   l [[103   3  80   1]
        # 00033:    [172  74 164   1]]

        # Coords in p' (unit cell system of reference)
        # this x_, y_, z_ are auxiliary ranges that will be
        # used to calculate a meshgrid. The purpose of meshgrid is to
        # create a rectangular grid out of an array of x_, y_ and z_
        # values.
        #
        # So, for example, if we want to create a grid where we have
        # a point at each integer value between 0 and 4 in x, y and z
        # directions. To create a rectangular grid, we need every
        # combination of the x, y and z points.
        x_ = np.arange(l[0,0], l[1,0]+1, 1); xdim = np.shape(x_)[0]
        y_ = np.arange(l[0,1], l[1,1]+1, 1); ydim = np.shape(y_)[0]
        z_ = np.arange(l[0,2], l[1,2]+1, 1); zdim = np.shape(z_)[0]
        # x_ x_.shape [103 104 105 106] [169 170 171 172] (70,)
        # y_ y_.shape [3 4 5 6] [71 72 73 74] (72,)
        # z_ z_.shape [80 81 82 83] [161 162 163 164] (85,)
        # mesh grid with all the points that need to be evaluated

        coords = np.meshgrid(x_, y_, z_, indexing='ij')
        # coords coords.shape [array([[[103, 103, 103, ..., 103, 103, 103],
        #         [103, 103, 103, ..., 103, 103, 103],
        #         [103, 103, 103, ..., 103, 103, 103],
        #         ...,
        #         ...,
        #         [ 80,  81,  82, ..., 162, 163, 164],
        #         [ 80,  81,  82, ..., 162, 163, 164],
        #         [ 80,  81,  82, ..., 162, 163, 164]]])]
        #  (3, 70, 72, 85)

        # stack the meshgrid to position vectors
        # meshgrid are several lists with x, y and z coordinates. We want
        # this values as vectors so we can multiply them by the transformation
        # matrix
        # Since I use generalized coordinates, vectors should have
        # 4 dimensions and the
        # last dimension should be always 1 (integer 1!!).
        xyz=np.vstack([coords[0].reshape(-1),
                       coords[1].reshape(-1),
                       coords[2].reshape(-1),
                       np.ones((xdim, ydim, zdim), dtype=int).reshape(-1)])
        #xyz xyz.shape [[103. 103. 103. ... 172. 172. 172.]
        # [  3.   3.   3. ...  74.  74.  74.]
        # [ 80.  81.  82. ... 162. 163. 164.]
        # [  1.   1.   1. ...   1.   1.   1.]] (4, 428400=  70 *  72 * 85)
        # let us create a volume with the same size than the input volume
        # to store the output
        outmatrix = np.zeros_like(matrix)

        # I can not use xyz as indexes for numpy since xys is in
        # generalized coordiantes
        x = xyz[0, :]
        y = xyz[1, :]
        z = xyz[2, :]
        # copy bounding box region to output volume
        # note that
        outmatrix[z,y,x] = matrix[z,y,x]

        # now copy the region related by the transformation matrices
        # interpolating the points if needed
        for counter, transformationMatrix in enumerate(self.listTransformationMatrices):
            transformed_xyz = transformationMatrix @ xyz
            # transformationMatrix
            # [[ 9.95590478e-01 -2.80233170e-03  9.37643145e-02  3.96147e+01]
            #  [-3.44402543e-03  9.97787836e-01  6.63895581e-02 -1.77753e+01]
            #  [-9.37429380e-02 -6.64197386e-02  9.93378417e-01 -8.00067e+00]
            #  [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  7.35294e-01]]
            # transformed_xyz transformed_xyz.shape
            # [[165.42162908 165.51539339 165.60915771 ... 241.60708031 241.70084462
            # 00183:     241.79460894]
            # 00184:    [  1.42342186   1.48981142   1.55620098 ...  77.47266421  77.53905377
            # 00185:      77.60544333]
            # 00186:    [ 85.10806994  86.10144836  87.09482677 ... 155.38103599 156.37441441
            # 00187:     157.36779283]
            # 00188:    [  1.           1.           1.         ...   1.           1.
            # 00189:       1.        ]] (4, 428400)

            # matrix(xp, yp and zp) should be averaged with matrix(x,y,z)
            xp = transformed_xyz[0, :]
            yp = transformed_xyz[1, :]
            zp = transformed_xyz[2, :]
            outmatrix[z,y,x] += self.trilinearInterpolation(xp, yp, zp, matrix)

        # divide sum by the number averaged regions
        factor = len(self.listTransformationMatrices) + 1.
        outmatrix /= factor
        img.setData(outmatrix)
        kkFn = self.outputFileName()
        img.write(kkFn)
        origin = volume.getOrigin(force=True).getShifts()
        Ccp4Header.fixFile(kkFn, kkFn, origin, sampling,
                           Ccp4Header.START)

    def createOutput(self):
        self.computeAveragedSubUnit()
        # save results
        inVol = self.inputVolume.get()
        outVol = Volume()
        origin = inVol.getOrigin(force=True)
        sampling = inVol.getSamplingRate()
        outVol.setSamplingRate(sampling)
        outVol.setOrigin(origin)
        outVol.setFileName(self.outputFileName())
        self._defineOutputs(outVol=outVol)
        with open("chimera_output.cxc", 'w') as f:
            dim = 150.
            bildFileName = os.path.abspath(self._getExtraPath(
                "axis_output.bild"))
            Chimera.createCoordinateAxisFile(dim,
                                             bildFileName=bildFileName,
                                             sampling=sampling)
            fnCmd = self._getExtraPath("chimera_output.cxc")
            f = open(fnCmd, 'w')
            # axis
            f.write("open %s\n" % bildFileName)
            f.write("cofr 0,0,0\n")  # set center of coordinates
            # 3D map
            f.write("open %s\n" % os.path.abspath(inVol.getFileName()))
            x_input, y_input, z_input = inVol.getShiftsFromOrigin()
            f.write("volume #2  voxelSize %f origin %0.2f,%0.2f,%0.2f\n" %
                    (sampling, x_input, y_input, z_input))
            f.write("open %s\n" % os.path.abspath(
                self.outputFileName()))
            x_input, y_input, z_input = outVol.getShiftsFromOrigin()
            f.write("volume #3  voxelSize %f origin %0.2f,%0.2f,%0.2f\n" %
                    (sampling, x_input, y_input, z_input))

            # reference pdb
            f.write("open %s\n" % os.path.abspath(
                self.atomStructReference.get().getFileName()))
            # other pdb
            for aStruct in self.otherAtomStructs:
                f.write("open %s\n" % os.path.abspath(
                    aStruct.get().getFileName()))

            bildFileFn = os.path.abspath(
                self._getExtraPath("chimera_output.bild"))
            with open(bildFileFn, 'w') as ff:
                ff.write(".transparency 0.8\n")
                # bounding box reference
                l = self.boundingBoxA
                ff.write(".color red\n")
                ff.write(".box %f %f %f %f %f %f\n" %
                         (l[0][0], l[0][1], l[0][2],
                          l[1][0], l[1][1], l[1][2]))
                ff.write(".transparency 0.\n")
            f.write("open %s\n" % bildFileFn)
            f.write("view\n")
            f.close()

    def _validate(self):
        errors = []
        # check number of chains in each entry
        atomStructReferenceFn = self.atomStructReference.get().getFileName()
        aStructReferenceHa = AtomicStructHandler(atomStructReferenceFn)
        listOfChains, listOfResidues = \
            aStructReferenceHa.getModelsChains()

        for aStruct in self.otherAtomStructs:
            atomStructFn = aStruct.get().getFileName()
            aStructReferenceHa = AtomicStructHandler(atomStructFn)
            listOfChains2, listOfResidues2 = \
                aStructReferenceHa.getModelsChains()
            if len(listOfChains[0]) != len(listOfChains2[0]):
                errors.append("Number of chains in reference atomic struct"
                              " and atomic struct called %s is differente"
                              " %d != %d" % (atomStructFn, len(listOfChains[0]),
                                             len(listOfChains2[0])))
        return errors

    def trilinearInterpolation(self, xf, yf, zf, input_array):
        """ Interpolates the value of point "indices" in matrix
            input_array.
            return the interpolated value
        """
        #output = np.empty(indices[0].shape)
        #x_indices = indices[0]
        #y_indices = indices[1]
        #z_indices = indices[2]

        x0 = xf.astype(np.integer)
        y0 = yf.astype(np.integer)
        z0 = zf.astype(np.integer)
        x1 = x0 + 1
        y1 = y0 + 1
        z1 = z0 + 1

        # Check if xyz1 is beyond array boundary:
        x1[np.where(x1 == input_array.shape[0])] = x0.max()
        y1[np.where(y1 == input_array.shape[1])] = y0.max()
        z1[np.where(z1 == input_array.shape[2])] = z0.max()

        x = xf - x0
        y = yf - y0
        z = zf - z0
        output = (input_array[z0, y0, x0] * (1 - x) * (1 - y) * (1 - z) +
                  input_array[z1, y0, x0] * x * (1 - y) * (1 - z) +
                  input_array[z0, y1, x0] * (1 - x) * y * (1 - z) +
                  input_array[z0, y0, x1] * (1 - x) * (1 - y) * z +
                  input_array[z1, y0, x1] * x * (1 - y) * z +
                  input_array[z0, y1, x1] * (1 - x) * y * z +
                  input_array[z1, y1, x0] * x * y * (1 - z) +
                  input_array[z1, y1, x1] * x * y * z)

        return output

    def _summary(self):
        summary = []
        for rms, objId in zip(self.listRMS, self.listObjID):
            summary.append("Obj: %d has rms %d", objId, rms)
        return summary

    def _citations(self):
        return []

    def _methods(self):
        return []
