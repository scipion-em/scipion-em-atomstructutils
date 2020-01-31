# **************************************************************************
# *
# * Authors:  Roberto Marabini (roberto@cnb.csic.es), May 2013
# *           Marta Martinez (mmmtnez@cnb.csic.es)
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

from pwem.emlib.image import ImageHandler
from atomstructutils.protocols.protocol_convertsymmetry import \
    ProtAtomStrucConvertSymmetry
from pwem.viewers import Chimera
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer
from pwem import Domain


class ProtAtomStrucConvertSymmetryViewer(Viewer):
    """ Visualize the output of protocol protocol_convertsymmetry """
    _environments = [DESKTOP_TKINTER]
    _label = 'convertsymmetry viewer'
    _targets = [ProtAtomStrucConvertSymmetry]

    def _visualize(self, obj, **args):
        # The input pdb is a parameter from the protocol
        # and from the parent protocol.
        inputAtomStruct = self.protocol.pdbFileToBeRefined.get()

        # To show pdbs only
        dim = 150.
        sampling = 1.

        bildFileName = os.path.abspath(self.protocol._getTmpPath(
            "axis_output.bild"))
        Chimera.createCoordinateAxisFile(dim,
                                 bildFileName=bildFileName,
                                 sampling=sampling)
        fnCmd = self.protocol._getTmpPath("chimera_output.cmd")
        f = open(fnCmd, 'w')
        f.write("open %s\n" % bildFileName)
        f.write("cofr 0,0,0\n")  # set center of coordinates
        f.write("open %s\n"
                % os.path.abspath(inputAtomStruct.getFileName()))
        if self.protocol.hasAttribute('rotatedAtomStruct'):
            outputAtomStruct = self.protocol.rotatedAtomStruct.getFileName()
            f.write("open %s\n" % os.path.abspath(outputAtomStruct))

        f.close()

        # run in the background
        chimeraPlugin = Domain.importFromPlugin('chimera', 'Plugin', doRaise=True)
        chimeraPlugin.runChimeraProgram(chimeraPlugin.getProgram(), fnCmd + "&")
        return []