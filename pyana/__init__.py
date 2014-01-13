"""

This package_ provides a sample setup to get started with ccsgp_. `ccsgp` is
initialized as a module and its usage demonstrated with dedicated functions in
the examples module. Helpful utility functions are also included to complement
the features of `ccsgp`. The following use cases are currently implemented or
will be in the future:

* keep input/output directory layouts according to the module's structure
* load data from text files with the correct column formats
* import data from ROOT objects via root-py
* and more ...

Installation
^^^^^^^^^^^^

1. `ccsgp_get_started` requires gnuplot, git, virtualenv and optionally hdf5.
   Install these dependencies via MacPorts, for instance::

     $ sudo port install gnuplot git-core py27-virtualenv [hdf5]

   `hdf5` is optional but if you'd like to save your image data to HDF5, you'll
   have to install it.

2. clone the `ccsgp_get_started` git repository.

   Include the submodule with the test data by cloning the git repository
   recursively::

     $ git clone https://github.com/tschaume/ccsgp_get_started.git --recursive
     $ cd ccsgp_get_started/

   or alternatively, for source code only, clone and init the `ccsgp`
   submodule separately::

     $ git clone https://github.com/tschaume/ccsgp_get_started.git
     $ cd ccsgp_get_started/
     $ git submodule init pyana/ccsgp
     $ git submodule update

3. init the virtualenv_, activate it and install all requirements::

     $ virtualenv-2.7 env
     $ source env/bin/activate
     $ pip install -U numpy
     $ pip install -U -r requirements.txt

   - Every time you start in a new terminal you have to activate the correct
     python environment by sourcing `env/bin/activate` again or instead use
     `env/bin/python` directly!

   - The h5py package is currently omitted from the requirements. If you want to
     use it, uncomment the h5py requirement in ``requirements.txt`` and rerun ``$
     pip install -r requirements.txt``.

4. If you intend to run the examples (and you cloned recursively above), symlink
   *pyanaDir* to *ccsgp_get_started_data*::

     $ ln -s ccsgp_get_started_data pyanaDir

   `pyana.aux.utils.checkSymLink` checks for the *pyanaDir* symbolic link and
   the code won't run without it. Hence you need to generate a symlink called
   *pyanaDir* either to *ccsgp_get_started_data* or your own input/output
   directory (preferably separated from the code repository).

.. _package: https://github.com/tschaume/ccsgp_get_started
.. _ccsgp: https://github.com/tschaume/ccsgp
.. _ccsgp_get_started: https://github.com/tschaume/ccsgp_get_started
.. _virtualenv: http://www.virtualenv.org/en/latest/virtualenv.html#usage
"""
