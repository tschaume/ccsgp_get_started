"""get started with the ccsgp plotting library"""

__version__ = '2.1.0'
__url__ = 'https://github.com/tschaume/ccsgp_get_started'
__author__ = 'Patrick Huck'
__email__ = 'patrick@the-huck.com'

"""

If the documentation here doesn't contain screenshots it didn't build correctly.
In that case, please find a static manually uploaded version at
http://downloads.the-huck.com/ccsgp_get_started/

This package_ provides a sample setup to get started with ccsgp_. `ccsgp` is
initialized as a module and its usage demonstrated with dedicated functions in
the examples module. Helpful utility functions are also included to complement
the features of `ccsgp`. The following use cases are currently implemented or
will be in the future:

* keep input/output directory layouts according to the module's structure
* load data from text files with the correct column formats
* import data from ROOT objects via root-py
* and more ...

Please submit tickets on GitHub issues_.

.. _issues: https://github.com/tschaume/ccsgp_get_started/issues

Installation
^^^^^^^^^^^^

1. `ccsgp_get_started` requires gnuplot, git, virtualenv and optionally hdf5.
   Install these dependencies via::

     $ sudo port install gnuplot git-core py27-virtualenv [hdf5] # MacPorts
     $ sudo apt-get install gnuplot git python-virtualenv [libhdf5-dev] # Debian/Ubuntu

   `hdf5` is optional but if you'd like to save your image data to HDF5, you'll
   have to install it.

2. clone the `ccsgp_get_started` git repository::

     $ git clone https://github.com/tschaume/ccsgp_get_started.git --recursive
     $ cd ccsgp_get_started/

3. init the virtualenv_, activate it and install all requirements::

     $ virtualenv-2.7 env
     $ source env/bin/activate
     $ pip install -U numpy
     $ pip install -U -r requirements.txt --allow-external gnuplot-py --allow-unverified gnuplot-py

   Every time you start in a new terminal you have to activate the correct
   python environment by sourcing `env/bin/activate` again or instead use
   `env/bin/python` directly!

   The h5py package is currently omitted from the requirements. If you want to
   use it, uncomment the h5py requirement in ``requirements.txt`` and rerun ``$
   pip install -r requirements.txt``.

4. As member of the STAR collaboration you can also pull in the protected data::

     $ cd data
     $ git remote add dielec http://cgit.the-huck.com/dielectron_data_protected
     $ git checkout -b dielec
     $ git pull dielec master
       (enter STAR protected credentials)


.. _package: https://github.com/tschaume/ccsgp_get_started
.. _ccsgp: https://github.com/tschaume/ccsgp
.. _ccsgp_get_started: https://github.com/tschaume/ccsgp_get_started
.. _virtualenv: http://www.virtualenv.org/en/latest/virtualenv.html#usage
"""
