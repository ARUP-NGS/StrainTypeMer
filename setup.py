from setuptools import setup, find_packages
import os
import subprocess


here = os.path.abspath(os.path.dirname(__file__))
# jf_path = os.path.join(here, "straintypemer", "jellyfish" )
# ## add code to get and install jellyfish tarball
# tarball_path = "https://github.com/gmarcais/Jellyfish/releases/download/v2.2.5/jellyfish-2.2.5.tar.gz"
# subprocess.check_call(["wget", "-P", jf_path, tarball_path])
# subprocess.check_call(["tar", "xvf", os.path.join(jf_path, os.path.basename(tarball_path)), "--directory", jf_path])
# os.remove(os.path.join(jf_path, os.path.basename(tarball_path)))
# jf_version_path = os.path.join(jf_path, os.listdir(jf_path)[0])
# subprocess.check_call([ jf_version_path, "configure" ,"--prefix={0}".format(jf_path),
#                         "--enable-python-binding={0}".format(here)])
# subprocess.check_call(["make", "-C", jf_version_path ])


setup(
    name='StrainTypeMer',
    version='1.0',
    packages=['straintypemer', 'straintypemer/sub_commands', ],
    package_dir={'straintypemer' : 'straintypemer', },
    package_data={'straintypemer' : ['data/*.jf',]},
    data_files = [('', ['straintypemer/data/dummy_A.jf'])],
    url='https://github.com/ARUP-NGS/StrainTypeMer',
    license='MIT',
    author='Keith E Simmon',
    author_email='ke.monk@gmail.com',
    description='kmer tool for strain typing',
    #packages = find_packages(),
    install_requires = ['numpy>=1.10', 'matplotlib>=1.5.0', 'biopython>=1.66', 'scipy', 'bson'],
    include_package_data = True,
)



