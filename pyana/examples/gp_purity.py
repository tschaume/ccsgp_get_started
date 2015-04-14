import os, argparse, logging, math
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_panel
from ..ccsgp.config import default_colors
import numpy as np

def gp_purity():
    inDir, outDir = getWorkDirs()
    energies = ['19', '27', '39', '62']
    particles = ['electrons', 'positrons']
    contams = ['{/Symbol \160}', 'K', 'p', '{/Symbol \160}{/Symbol \160}']
    dtypes = ['purity'] + ['contam{}'.format(i+1) for i in range(len(contams))]
    data = OrderedDict()
    for eidx,energy in enumerate(energies):
        ekey = ' '.join([getEnergy4Key(str(energy)), 'GeV'])
        data[ekey] = [[], [], []]
        for pidx,particle in enumerate(particles):
            for didx,dtype in enumerate(dtypes):
                infile = os.path.join(inDir, '{}_{}_{}.dat'.format(
                    particle, energy, dtype
                ))
                data_import = np.loadtxt(open(infile,'r'))
                data_import[:,2] = 0
                #data_import = data_import[data_import[:,1]>0.001]
                propsstr = 'with linespoints lt {} lw 3 pt 18 lc {}' if didx == 0 \
                        else 'with lines lt {} lw 3 lc {}'
                props = propsstr.format(pidx+1, default_colors[didx])
                data[ekey][0].append(data_import)
                data[ekey][1].append(props)
                title = contams[didx-1] if didx > 0 else 'e'
                title += '^{%s}' % ('+' if particle == 'positrons' else '-')
                if title == 'p^{-}': title = '@^{\261}p'
                if title == 'p^{+}': title = 'p'
                data[ekey][2].append(title)
    make_panel(
        dpt_dict = data, name = os.path.join(outDir, 'purity'),
        yr = [-0.03,1.03], xr = [0.15,1.85],
        xlabel = 'momentum, p (GeV/c)',
        ylabel = 'fraction of candidate sample (%)',
        layout = '2x2', size = '5in,7in', key = ['nobox', 'at graph 0.4,0.75'],
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
    gp_purity()
