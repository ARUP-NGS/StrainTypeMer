from setuptools import setup

setup(
    name='StrainTypeMer',
    version='1.0b',
    packages=['straintypemer', 'straintypemer/sub_commands', ],
    package_dir={'straintypemer' : 'straintypemer', },
    package_data={'straintypemer' : ['data/*.fa','data/*.txt','data/*.jf', ]},
    data_files = [('', ['straintypemer/data/dummy_A.jf', 'straintypemer/data/ARmeta-genes.fa',
                        'straintypemer/data/categories.txt', 'straintypemer/data/AROtags.txt',
                        ])],
    url='https://github.com/ARUP-NGS/StrainTypeMer',
    license='MIT',
    author='Keith E Simmon',
    author_email='ke.monk@gmail.com',
    description='kmer tool for strain typing',
    #packages = find_packages(),
    install_requires = ['numpy>=1.10', 'biopython>=1.66', 'matplotlib>=1.5.1'],
    include_package_data = True,
)



