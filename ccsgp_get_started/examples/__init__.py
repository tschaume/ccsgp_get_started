"""

The examples are based on a dataset of World Bank Indicators_.  You can use the
dataset yourself to play around [#]_. See the ``genExDat.sh`` script in the
same directory on how I extracted the data into the correct format for *ccsgp*.
To generate all example plots based on *ccsgp_get_started_data* you can run::

  $ python -m ccsgp_get_started
  
.. [#] ccsgp_get_started_data/input/examples/gp_datdir/{WorldBankIndicators.csv, genExDat.sh}

Alternatively, you can run a specific module, for instance::

  $ python -m ccsgp_get_started.examples.gp_datdir [--log] <country-initial> <#-most-populated>

and this way plot specific country initials. You can open all resulting pictures
via ``$ open examplesDir/examples/gp_datdir/*.pdf`` or use ``pdfnup`` to put
multiple plots on one page. To start on your own read the documentation below or
the source code and use one of the examples as a template.

.. _Indicators: http://www.tableausoftware.com/public/community/sample-data-sets#worldbank
"""
