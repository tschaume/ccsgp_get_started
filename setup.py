from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

setup(
  name = 'ccsgp_get_started',
  version = '2.0',
  description = 'get started with ccsgp',
  license = 'MIT',
  author = 'Patrick Huck',
  author_email = 'patrick@the-huck.com',
  url = 'https://github.com/tschaume/ccsgp_get_started',
  download_url = 'https://github.com/tschaume/ccsgp_get_started/archive/v2.0.tar.gz',
  keywords = ['gnuplot', 'graph', 'plot', 'panel'],
  packages = find_packages(),
  install_requires = [
    #'h5py==2.2.1',
    'pint==0.5.1', 'pymodelfit==0.1.2',
    'uncertainties==2.4.4', 'gnuplot-py>=1.8', 'numpy>=1.8.0', 'wsgiref==0.1.2'
  ]
)
