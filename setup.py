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
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='scipion-em-atomstructutils',
    version='3.0.1',
    description='A Scipion plugin to manipulate atomic structure files (PDB/MMCIF)',
    long_description=long_description,
    url='https://github.com/scipion-em/scipion-em-atomstructutils',
    author='Roberto Marabini and Marta Martinez',
    author_email='scipion@cnb.csic.es',
    keywords='scipion pdb  scipion-3',
    packages=find_packages(),
    install_requires=[requirements],
    include_package_data=True,
    package_data={
       'atom_struc_utils': ['tool.png'],
    },
    entry_points={
        'pyworkflow.plugin': 'atomstructutils = atomstructutils'
    }
)
