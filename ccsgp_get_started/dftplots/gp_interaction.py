import logging, argparse, os, sys
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_panel
from ..examples.utils import getWorkDirs
from ..ccsgp.utils import getOpts, colorscale
from ..ccsgp.config import default_colors
from string import ascii_lowercase

def gp_interaction():
  """example for plotting from a text file via numpy.loadtxt

  1. prepare input/output directories
  2. load the data into an OrderedDict() [adjust axes units]
  3. call ccsgp.make_panel with data from 2

  Below is an output image for country initial T and the 4 most populated
  countries for this initial (click to enlarge). Also see::

    $ python -m ccsgp_get_started.examples.gp_datdir -h

  for help on the command line options.

  .. image:: pics/T.png
     :width: 450 px

  .. image:: pics/U.png
     :width: 450 px

  :param initial: country initial
  :type initial: str
  :param topN: number of most populated countries to plot
  :type topN: int
  :ivar inDir: input directory according to package structure and initial
  :ivar outDir: output directory according to package structure
  :ivar data: OrderedDict with datasets to plot as separate keys
  :ivar file: data input file for specific country, format: [x y] OR [x y dx dy]
  :ivar country: country, filename stem of input file
  :ivar file_url: absolute url to input file
  :ivar nSets: number of datasets
  """
  # prepare color arrays for dataset
  my_color_set = [default_colors[i] for i in [0,2]]
  my_color_array = []
  for icol, color in enumerate(my_color_set):
      my_color_array.append(colorscale(color[-7:-1], 0.8))
      my_color_array.append(colorscale(color[-7:-1], 1.4))
  # prepare input/output directories
  inDir, outDir = getWorkDirs()
  dpt_dict = OrderedDict()
  for fold_idx,fold in enumerate(['totalE', 'Ex', 'Ec_LDA', 'Ec_vdW']):
    energy = os.path.splitext(fold)[0]
    key = '('+list(ascii_lowercase)[fold_idx]+')'
    # prepare data
    dpt_dict[key] = [[], [], []] #dpt
    for index,funct in enumerate(os.listdir(os.path.join(inDir, energy))):
      funct_name = os.path.splitext(funct)[0]
      file_url = os.path.join(inDir, energy, funct)
      data_import = np.loadtxt(open(file_url, 'rb')) #load data
      data_import[:,0] /= 10.
      dpt_dict[key][0].append(data_import)
      dpt_dict[key][1].append('with linespoints lw 4 pt 18 ps 1.5 lt 2 lc %s' %
                  my_color_array[index])
      dpt_dict[key][2].append(funct_name)
  #print dpt_dict
  logging.debug(dpt_dict) # shown if --log flag given on command line
  # generate plot using ccsgp.make_plot
  make_panel(
    dpt_dict = dpt_dict,
    name = os.path.join(outDir, 'interaction'),
    xr = [1.8, 4.53], yr = [-55, 20],
    key = ['bottom right', 'maxrows 2', 'width -1.1', 'nobox' ],
    ylabel = 'DFT binding energy (kJ/mol)',
    xlabel = 'Mg-O distance, r (nm)', rmargin = 0.98, size='7in,7.5in',
    lines = {'x=0': 'lc {} lw 4 lt 1'.format(default_colors[-10])},
    debug = False, layout = '2x2', key_subplot_id = 2
  )
  return 'done'

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  #parser.add_argument("initial", help="country initial = input subdir with txt files")
  #parser.add_argument("topN", help="number of most populated countries to plot")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_interaction()
