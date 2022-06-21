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
    packages = [
        package_name, '{}.ccsgp'.format(package_name),
        '{}.ccsgp.Gnuplot'.format(package_name),
        '{}.examples'.format(package_name)
    ],
    include_package_data=True,
    install_requires = [ 'numpy==1.22.0', 'Pint==0.5.1', 'uncertainties==2.4.4' ],
    license = 'MIT',
    keywords = ['gnuplot', 'graph', 'plot', 'panel'],
)
