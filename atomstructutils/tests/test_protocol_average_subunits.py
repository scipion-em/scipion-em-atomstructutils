# ***************************************************************************
# * Authors:    Marta Martinez (mmmtnez@cnb.csic.es)
# *             Roberto Marabini (roberto@cnb.csic.es)
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

# protocol to test the operation on PDB files
import os
from atomstructutils.protocols import ProtAtomStrucOperate, ProtAverageSubunits
from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols.protocol_import import ProtImportPdb
from pwem.protocols.protocol_import.volumes import ProtImportVolumes
import pwem.protocols as emprot
from xmipp3.protocols import XmippProtExtractUnit
from xmipp3.constants import (XMIPP_I222r)


class TestImportBase(BaseTest):

    @classmethod
    def _importStructurePDB(cls, pdbID):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': pdbID
                }
        protImportPDB = cls.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import structure\n %s' % pdbID)
        cls.launchProtocol(protImportPDB)

        cls.assertTrue(cls,
            protImportPDB.outputPdb.getFileName().find(pdbID)!= -1,
                       "Error: Downloading atom struct id=%s" % pdbID)
        return protImportPDB.outputPdb

    @classmethod
    def _extractChain(cls, atomStruct):
        # extract all chain
        _dictOperations = ProtAtomStrucOperate.operationsDictInv
        args = {'pdbFileToBeRefined': atomStruct,
                'Operation': _dictOperations['extractAllChains']
                }
        protAtomStrucExtractChain = cls.newProtocol(ProtAtomStrucOperate, **args)
        protAtomStrucExtractChain.setObjLabel('split chains\n')

        cls.launchProtocol(protAtomStrucExtractChain)
        chainA = protAtomStrucExtractChain.outputPdb_chainA
        chainB = protAtomStrucExtractChain.outputPdb_chainB
        chainC = protAtomStrucExtractChain.outputPdb_chainC
        chainD = protAtomStrucExtractChain.outputPdb_chainD
        chainE = protAtomStrucExtractChain.outputPdb_chainE
        chainF = protAtomStrucExtractChain.outputPdb_chainF
        chainG = protAtomStrucExtractChain.outputPdb_chainG
        chainH = protAtomStrucExtractChain.outputPdb_chainH
        chainI = protAtomStrucExtractChain.outputPdb_chainI
        chainJ = protAtomStrucExtractChain.outputPdb_chainJ
        chainK = protAtomStrucExtractChain.outputPdb_chainK
        chainL = protAtomStrucExtractChain.outputPdb_chainL
        chainM = protAtomStrucExtractChain.outputPdb_chainM
        cls.assertTrue(cls, chainA.getFileName().find("Chain_A")!= -1,
                       "Error: extracting chains")
        return chainA, chainB, chainC, chainD, chainE, chainF,\
               chainG, chainH, chainI, chainJ, chainK, chainL,\
               chainM

    @classmethod
    def _joinChains(cls, chainList, title="join chains"):
        "Chain union"
        _dictOperations = ProtAtomStrucOperate.operationsDictInv
        args = {'pdbFileToBeRefined': chainList[0],
                'InputAtomStruct2': chainList[1:],
                'Operation': _dictOperations['addChain']
                }
        protAtomStrucExtractChain = cls.newProtocol(ProtAtomStrucOperate, **args)
        protAtomStrucExtractChain.setObjLabel(title)
        cls.launchProtocol(protAtomStrucExtractChain)
        chains = protAtomStrucExtractChain.outputPdb
        cls.assertTrue(cls, chains.getFileName().find("addChain")!= -1,
                       "Error: adding chains")
        return chains

    @classmethod
    def _importVolume(cls, emdbId):
        # Import 3D map from EMDB
        args = {'importFrom': ProtImportVolumes.IMPORT_FROM_EMDB,
                'emdbId': emdbId
                }

        protImportVol = cls.newProtocol(emprot.ProtImportVolumes, **args)
        protImportVol.setObjLabel('import vol %d,\n from EMDB' % emdbId)
        cls.launchProtocol(protImportVol)
        volume = protImportVol.outputVolume
        cls.assertTrue(cls, volume.getFileName().find("%d" % emdbId)!= -1,
                       "Error: adding chains")

        #extract unit cell
        args = {'inputVolumes': volume,
                'symmetryGroup': XMIPP_I222r,
                'symmetryOrder': 1,
                'innerRadius': 218,
                'outerRadius': 380,
                'expandFactor': .1,
                'offset': 0
                }

        protExtractUnit = cls.newProtocol(XmippProtExtractUnit, **args)
        protExtractUnit.setObjLabel('extract unit cell')
        cls.launchProtocol(protExtractUnit)
        unitCell = protExtractUnit.outputVolume
        cls.assertTrue(cls, unitCell.getFileName().find("output_volume")!= -1,
                       "Error: adding chains")
        return unitCell

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        # get atomic structure
        pdb = cls._importStructurePDB('6yba')
        # split atom struct in chains
        chainA, chainB, chainC, \
        chainD, chainE, chainF, \
        chainG, chainH, chainI, \
        chainJ, chainK, chainL, \
        cls.chainM = cls._extractChain(pdb)
        # recreate hexon1, hexon2, hexon3, hexon4
        cls.h1 = cls._joinChains([chainA, chainB, chainC], "h1")
        cls.h2 = cls._joinChains([chainD, chainE, chainF], "h2")
        cls.h3 = cls._joinChains([chainG, chainH, chainI], "h3")
        cls.h4 = cls._joinChains([chainJ, chainK, chainL], "h4")
        # import corresponding volume
        cls.unitcell = cls._importVolume(10768)

class TestAverageSubUnit(TestImportBase):

    def testAverageSubUnitDifferentNumberChains(self):
        """Call protocol with two atomic structs with different
        number of chains. This should produce an error"""
        # hexons 1 and 2
        h1 = self.h1; h2 = self.h2
        # pentom
        pentom = self.chainM
        # unitcell
        unitCell = self.unitcell

        # call new protocol with two atomStructs
        # with a different number of chains

        args = {'inputVolume': unitCell,
                'atomStructReference': h1,
                'otherAtomStructs': [pentom],
                'rangeStart': -1,
                'rangeEnd': -1,
                }

        protAverageSubunits = self.newProtocol(ProtAverageSubunits, **args)
        protAverageSubunits.setObjLabel('This protocol should FAIL\n'
                                        'invalid input parameter\n'
                                        ' diff num chains.')
        print("This test should raise the error message "
              "'Number of chains in reference atomic struct...'"
              " because the provided parameters are wrong")
        try:
            self.launchProtocol(protAverageSubunits)
        except Exception as e:
            print('Exception catched, this is OK '
                  'since we have provided wrong parameters')
            self.assertTrue(True)
            return
        self.assertTrue(False)

    def testAverageSubUnit(self):
        """Call protocol with four hexons"""
        # hexons 1, 2, 3, 4
        h1 = self.h1; h2 = self.h2
        h3 = self.h3; h4 = self.h4

        # unitcell
        unitCell = self.unitcell

        # call new protocol with four atomStructs
        # with a different number of chains

        args = {'inputVolume': unitCell,
                'atomStructReference': h1,
                'otherAtomStructs': [h2, h3, h4],
                'rangeStart': -1,
                'rangeEnd': -1,
                }

        protAverageSubunits = self.newProtocol(ProtAverageSubunits, **args)
        protAverageSubunits.setObjLabel('Four hexons')
        self.launchProtocol(protAverageSubunits)

        # let us check that there is an exit file
        result= protAverageSubunits.outVol.getFileName()
        self.assertTrue(os.path.exists(result))
