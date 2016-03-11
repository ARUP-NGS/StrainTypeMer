from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='StrainTypeMer',
    version='1.0',
    packages=['straintypemer', 'straintypemer/sub_commands',],
    package_dir={'straintypemer' : 'straintypemer' },
    package_data={'straintypemer' : ['data/*.jf',]},
    data_files = [('', ['straintypemer/data/dummy_A.jf'])],
    url='https://github.com/ARUP-NGS/StrainTypeMer',
    license='MIT',
    author='Keith E Simmon',
    author_email='ke.monk@gmail.com',
    description='kmer tool for strain typing',
    #packages = find_packages(),
    install_requires = ['numpy>=1.10', 'matplotlib>=1.5.0', 'biopython>=1.66'],
    include_package_data = True,
)

