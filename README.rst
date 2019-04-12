================================
Atom_struct_utils scipion plugin
================================

Utilities for handling PDB/mmcif atomic structure files:

So far it contains two utilities:

* extract a chain from an atom struct (pdb/cif file) 
* perform union of two atomic structs

implemented in a single protocol called `operate`


===================
Install this plugin
===================

You will need to use `2.0.0 <https://github.com/I2PC/scipion/releases/tag/v2.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

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


Atom_struct_utils plugin is a pure Python module, no binary files are required. 

- **TESTS**

To check the installation, simply run the following Scipion test: `scipion test  --grep TestOperate --run`


========
Examples
========

See `Model Building Tutorial <https://github.com/I2PC/scipion/wiki/tutorials/tutorial_model_building_basic.pdf>`_



===============
Buildbot status
===============

Status devel version: 

.. image:: http://arquimedes.cnb.csic.es:9980/badges/atomstructutils_devel.svg

Status production version: 

.. image:: http://arquimedes.cnb.csic.es:9980/badges/atomstructutils_prod.svg
