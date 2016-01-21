from setuptools import setup, find_packages

setup(
    name='straintypemer',
    version='1.0b',
    packages=['straintypemer/sub_commands/',],
    url='https://github.com/ARUP-NGS/StrainTypeMer',
    license='MIT',
    author='Keith E Simmon',
    author_email='ke.monk@gmail.com',
    description='kmer tool for strain typing',
    #packages = find_packages(),
    install_requires = ['numpy>=1.10', 'matplotlib>=1.5.0', 'biopython>=1.66'],

    #package_data = { '/straintypemer/*/' : ['*.txt', '*.tfa'],}


)
