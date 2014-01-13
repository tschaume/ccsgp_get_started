from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages
setup(
  name = 'ccsgp_get_started',
  version = '1.0',
  description = 'get started with ccsp',
  license = 'MIT',
  author = 'Patrick Huck',
  author_email = 'phuck@lbl.gov',
  url = 'https://github.com/tschaume/ccsgp_get_started',
  packages = find_packages(),
  install_requires = [
    'wsgiref==0.1.2', 'numpy==1.8.0', 'gnuplot-py==1.8'#, 'h5py==2.2.1'
  ]
)
