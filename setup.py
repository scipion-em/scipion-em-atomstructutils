"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='scipion-em-atomstructutils',
    version='1.0.0b',
    description='A Scipion plugin to manipulate atomic structure files (PDB/MMCIF)',
    long_description=long_description,
    url='https://github.com/scipion-em/scipion-em-atomstructutils',
    author='Roberto Marabini and Marta Martinez',
    author_email='scipion@cnb.csic.es',
    keywords='scipion pdb  scipion-1.2',
    packages=find_packages(),
    install_requires=[],
    package_data={
       'atom_struc_utils': ['tool.png'],
    }
)
