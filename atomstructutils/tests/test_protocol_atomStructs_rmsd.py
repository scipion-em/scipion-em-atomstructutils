# ***************************************************************************
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from atomstructutils.protocols import ProtRMSDAtomStructs
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols.protocol_import import ProtImportSetOfAtomStructs


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
      setupTestProject(cls)
      cls.ds = DataSet.getDataSet('model_building_tutorial')

    @classmethod
    def _importStructurePDBSet(cls):
        args = {'inputPdbData': ProtImportSetOfAtomStructs.IMPORT_FROM_FILES,
                'filesPath': cls.ds.getFile('PDBx_mmCIF'),
                'filesPattern': cls.ds.getFile('PDBx_mmCIF/1ake_*.pdb')
                }
        protImportPDBs = cls.newProtocol(ProtImportSetOfAtomStructs, **args)
        protImportPDBs.setObjLabel('import 1akes structures\n')
        cls.launchProtocol(protImportPDBs)

        return getattr(protImportPDBs, protImportPDBs._OUTNAME)

class TestRMSD(TestImportBase):

    def testRMSD(self):
        """"""
        inputASs = self._importStructurePDBSet()
        args = {'inputStructureSet': inputASs,
                'chains': '',
                'considerAtoms': 2,
                }

        protRMSDAll = self.newProtocol(ProtRMSDAtomStructs, **args)
        protRMSDAll.setObjLabel('RMSD for all atoms')
        self.launchProtocol(protRMSDAll)
