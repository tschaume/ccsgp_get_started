import io, re
from setuptools import setup

package_name = 'ccsgp_get_started'
init_py = io.open('{}/__init__.py'.format(package_name)).read()
metadata = dict(re.findall("__([a-z]+)__ = '([^']+)'", init_py))
metadata['doc'] = re.findall('"""(.+)"""', init_py)[0]

setup(
    name = package_name,
    version = metadata['version'],
    description = metadata['doc'],
    author = metadata['author'],
    author_email = metadata['email'],
    url = metadata['url'],
    packages = [package_name],
    dependency_links = [
        'http://sourceforge.net/projects/gnuplot-py/files/Gnuplot-py/1.8/gnuplot-py-1.8.tar.gz/download'
    ],
    setup_requires = [ 'numpy==1.9.2' ],
    install_requires = [
        'numpy==1.9.2', 'Pint==0.5.1',
        'PyModelFit==0.1.2', 'uncertainties==2.4.4'
    ],
    license = open('LICENSE').read(),
    keywords = ['gnuplot', 'graph', 'plot', 'panel'],
)
