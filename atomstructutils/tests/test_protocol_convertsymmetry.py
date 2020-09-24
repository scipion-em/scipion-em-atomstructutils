# ***************************************************************************
# * Authors:    Roberto Marabini (roberto@cnb.csic.es)
# *
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
# ***************************************************************************/

# protocol to test the operation on PDB files
from atomstructutils.protocols import ProtAtomStrucConvertSymmetry
from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols.protocol_import import ProtImportPdb
from pwem.convert.symmetry import Icosahedron
from pwem.constants import (SYM_I222, SYM_I222r, SYM_In25, SYM_In25r,
                                     SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r,
                                     SCIPION_SYM_NAME)
from pwem.convert.atom_struct import AtomicStructHandler

class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

class TestConvertSymmetry(TestImportBase):

    def createPDBFile(self):
        "Create test CIF file with 12 Atoms in icosahedron vertexes"
        from Bio.PDB.Structure import Structure
        from Bio.PDB.Model import Model
        from Bio.PDB.Chain import Chain
        from Bio.PDB.Residue import Residue
        from Bio.PDB.Atom import Atom
        from Bio.PDB.mmcifio import MMCIFIO
        import os
        CIFFILENAME = "/tmp/out.cif"

        # create atom struct with ico simmety (i222r)
        icosahedron = Icosahedron(circumscribed_radius = 100,
                                  orientation='222r')
        pentomVectorI222r = icosahedron.getVertices()

        # create biopython object
        structure = Structure('result')  # structure_id
        model=Model(1,1)  # model_id,serial_num
        structure.add(model)
        chain=Chain('A')  # chain Id
        model.add(chain)
        for i, v in enumerate(pentomVectorI222r, 1):
            res_id = (' ', i, ' ')  # first arg ' ' -> aTOm else heteroatom
            res_name = "ALA" #+ str(i)  # define name of residue
            res_segid = '    '
            residue = Residue(res_id, res_name, res_segid)
            chain.add(residue)
            # ATOM name, coord, bfactor, occupancy, altloc, fullname, serial_number,
            #             element=None)
            atom=Atom('CA',v,    0.,            1.,     " ",    " CA ",  i,
                      "C")
            residue.add(atom)

        io = MMCIFIO()
        io.set_structure(structure)
        # delete file if exists
        if os.path.exists(CIFFILENAME):
            os.remove(CIFFILENAME)
        io.save(CIFFILENAME)
        return CIFFILENAME

    def _importStructurePDB(self, pdbFileName):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbId': pdbFileName,
                'pdbFile': pdbFileName
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n %s' % pdbFileName)
        self.launchProtocol(protImportPDB)
        return protImportPDB.outputPdb

    def compare(self, pdb, sym):

        args = {'pdbFileToBeRefined': pdb,
                'originSymmetryGroup': SYM_I222r - SYM_I222,
                'targetSymmetryGroup': sym - SYM_I222
                }

        protAtomStrucOperate = self.newProtocol(ProtAtomStrucConvertSymmetry, **args)
        protAtomStrucOperate.setObjLabel('rotate atom structs, to %s'%SCIPION_SYM_NAME[sym])
        self.launchProtocol(protAtomStrucOperate)

        aSH = AtomicStructHandler(protAtomStrucOperate.rotatedAtomStruct.getFileName())
        atoms_coord = [atom.coord for atom in aSH.getStructure().get_atoms()]
        icosahedron = Icosahedron(circumscribed_radius = 100,
                                  orientation=SCIPION_SYM_NAME[sym][1:])
        pentomVector = icosahedron.getVertices()
        for atom, vertex in zip(atoms_coord, pentomVector):
            for a,v in zip(atom, vertex):
                self.assertAlmostEqual(a, v, places=2)

    def testConvertSym(self):

        CIFFILENAME = self.createPDBFile()
        pdb = self._importStructurePDB(CIFFILENAME)

        self.compare(pdb, SYM_I222)
        self.compare(pdb, SYM_I222r)
        self.compare(pdb, SYM_In25)
        self.compare(pdb, SYM_In25r)
        self.compare(pdb, SYM_I2n3)
        self.compare(pdb, SYM_I2n3r)
        self.compare(pdb, SYM_I2n5)
        self.compare(pdb, SYM_I222r)