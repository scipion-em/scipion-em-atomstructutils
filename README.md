# Atom_struct_utils plugin

Utilities for handling PDB/mmcif atomic structure files:

So far it contains two utilities:

* extract a chain from an atom struct (pdb/cif file) 
* perform union of two atomic structs

implemented in a single protocol called `operate`

## Installation

You will need to use [2.0](https://github.com/I2PC/scipion/releases/tag/v2.0) version of Scipion to be able to run these protocols. To install the plugin, you have two options:

   a) Stable version
   ```
   scipion installp -p scipion-em-atom_struct_utils
   ```
   b) Developer's version
   * download repository 
   ```
    git clone https://github.com/scipion-em/scipion-em-atom_struct_utils.git
   ```
   * install 
   ```
    scipion installp -p path_to_scipion-em-atom_struct_utils --devel
   ```

Atom_struct_utils plugin is a pure Python module, no binary files are required. 

To check the installation, simply run the following Scipion test: `scipion test  --grep TestOperate --run`

## Examples
[See Model Building Tutorial](https://github.com/I2PC/scipion/wiki/tutorials/tutorial_model_building_basic.pdf)


## BuildBot Status
Status devel version: ![build status](http://arquimedes.cnb.csic.es:9980/badges/atomstructutils_devel.svg "Build status")

Status production version: ![build status](http://arquimedes.cnb.csic.es:9980/badges/atomstructutils_prod.svg "Build status")
