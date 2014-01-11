
# TODO: use package __main__ to run all analyses
# with default arguments

import string
from aux.utils import checkSymLink
from examples.gp_datdir import gp_datdir

checkSymLink()
for l in list(string.ascii_uppercase):
  print gp_datdir(l)
