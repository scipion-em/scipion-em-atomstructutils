# ***************************************************************************
# * Authors:    Marta Martinez (mmmtnez@cnb.csic.es)
# *             Roberto Marabini (roberto@cnb.csic.es)
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
from os.path import exists
from collections import Counter

from atomstructutils.protocols import ProtAtomStrucOperate
from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols.protocol_import import ProtImportPdb
from pwem.convert.atom_struct import AtomicStructHandler


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

class TestOperate(TestImportBase):
    #import 2 proteins

    def _importStructurePDB(self, pdbID):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': pdbID
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import structure\n %s' % pdbID)
        self.launchProtocol(protImportPDB)
        self.assertTrue(protImportPDB.outputPdb.getFileName())
        return protImportPDB.outputPdb

    def testAddChain(self):
        pdb1 = self._importStructurePDB('1P30') # A
        pdb2 = self._importStructurePDB('5NI1') # A, B, C, D
        pdb3 = self._importStructurePDB('1J77') # A
        _dictOperations = ProtAtomStrucOperate.operationsDictInv
        args = {'pdbFileToBeRefined': pdb1,
                'InputAtomStruct2': [pdb2, pdb3],
                'Operation': _dictOperations['addChain']
                }

        protAtomStrucOperate = self.newProtocol(ProtAtomStrucOperate, **args)
        protAtomStrucOperate.setObjLabel('add atom structs')
        self.launchProtocol(protAtomStrucOperate)
        outPutPDB = protAtomStrucOperate.outputPdb.getFileName()

        # check file exists
        self.assertTrue(exists(protAtomStrucOperate.outputPdb.getFileName()),
                        "Filename {} does not exists".format(outPutPDB))

        aSH = AtomicStructHandler(outPutPDB)
        chains = [chain.id for chain in aSH.getStructure().get_chains()]
        goal = ['A', 'A002', 'B', 'C', 'D', 'A003']
        self.assertTrue(Counter(chains) == Counter(goal),
                        "{} != {}".format(chains, goal))

        # atoms are OK
        aSH1 = AtomicStructHandler(pdb1.getFileName())
        aSH2 = AtomicStructHandler(pdb2.getFileName())
        aSH3 = AtomicStructHandler(pdb3.getFileName())
        #
        atomsNum1 = len([atom.id for atom in aSH1.getStructure().get_atoms()])
        atomsNum2 = len([atom.id for atom in aSH2.getStructure().get_atoms()])
        atomsNum3 = len([atom.id for atom in aSH3.getStructure().get_atoms()])
        atomsNumT = len([atom.id for atom in aSH.getStructure().get_atoms()])
        self.assertEqual(atomsNum1 + atomsNum2 + atomsNum3, atomsNumT)


    def testExtractChain(self):
        pdb1 = self._importStructurePDB('1P30') # A, B, C
        _dictOperations = ProtAtomStrucOperate.operationsDictInv
        args = {'pdbFileToBeRefined': pdb1,
                'Operation': _dictOperations['extractChain'],
                'inputStructureChain': '{"model": 0, "chain": "A", "residues": 891}',
                'start': 0,
                'end': 20
                }

        protAtomStrucExtractChain = self.newProtocol(ProtAtomStrucOperate, **args)
        self.launchProtocol(protAtomStrucExtractChain)
        aSH2 = AtomicStructHandler()
        outFileName = protAtomStrucExtractChain.outputPdb.getFileName()
        aSH2.read(outFileName)
        atomsNumT = len([atom.id for atom in aSH2.getStructure().get_atoms()])
        self.assertEqual(122, atomsNumT)
