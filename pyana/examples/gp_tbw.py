import os, argparse, logging, math
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot, make_panel
from ..ccsgp.config import default_colors
import numpy as np

def gp_tbw():
    inDir, outDir = getWorkDirs()
    energies = ['19', '27', '39', '62']
    particles = OrderedDict()
    particles['{/Symbol \160}'] = ['pip', 'pim']
    particles['K'] = ['kp', 'km']
    particles['p'] = ['p', 'pbar']
    dtypes = ['tbw', 'data']
    data = OrderedDict()
    for particle,codes in particles.iteritems():
        if particle not in data:
            data[particle] = [[], [], []]
        for energy in energies:
            for code in codes:
                for dtype in dtypes:
                    filename = '_'.join([dtype, code, energy]) + '.dat'
                    data_import = np.loadtxt(open(
                        os.path.join(inDir, filename), 'rb'
                    ))
                    data[particle][0].append(data_import)
                    data[particle][1].append('lt 1')
                    data[particle][2].append(filename)
    print data['p']
    make_panel(
        dpt_dict = data,
        name = os.path.join(outDir, 'tbw'),
        yr = [1e-6,100], xr = [0,2.2], ylog = True,
        xlabel = 'p_{T} (GeV/c)', ylabel = 'dN/2pipTdpTdy',
        layout = '3x1', size = '4.5in,12in', tmargin = 0.9,
        key = ['nobox', 'width -3', 'at screen 1.0,1.0', 'maxrows 2'],
    )


if __name__ == '__main__':
    checkSymLink()
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", help="show log output", action="store_true")
    args = parser.parse_args()
    loglevel = 'DEBUG' if args.log else 'WARNING'
    logging.basicConfig(
        format='%(message)s', level=getattr(logging, loglevel)
    )
    gp_tbw()
