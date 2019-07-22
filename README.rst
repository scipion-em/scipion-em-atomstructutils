================================
Atom_struct_utils scipion plugin
================================

Utilities for handling PDB/mmcif atomic structure files:

Current plugin protocols:

* operate (option extract): extract a chain from an atom structure (pdb/cif file).
* operate (option add): perform union of several atomic structures.
* convert_sym: takes an atomStruct in any icosahedric symemtry and rotates it to another.

===================
Install this plugin
===================

You will need to use `2.0.0 <https://github.com/I2PC/scipion/releases/tag/v2.0>`_ version of Scipion to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion installp -p scipion-em-atom_struct_utils
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. Download repository: 

.. code-block::

            git clone https://github.com/scipion-em/scipion-em-atom_struct_utils.git

2. Install:

.. code-block::

            scipion installp -p path_to_scipion-em-atom_struct_utils --devel

- **Binary files** 

Atom_struct_utils plugin is a pure Python module, no binary files are required. 

- **Tests**

To check the installation, simply run the following Scipion test:
    * `scipion test --grep TestOperate --run`
    * `scipion test --grep TestConvertSymmetry --run`

========
Examples
========

See `Model Building Tutorial <https://github.com/I2PC/scipion/wiki/tutorials/tutorial_model_building_basic.pdf>`_



===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/atomstructutils_devel.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/atomstructutils_prod.svg
