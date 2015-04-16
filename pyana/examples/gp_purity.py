import os, argparse, logging, math, glob
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_panel
from ..ccsgp.config import default_colors
import numpy as np

particles = ['electrons', 'positrons']
contams = ['{/Symbol \160}', 'K', 'p', '{/Symbol \160}{/Symbol \160}']

def gp_purity():
    inDir, outDir = getWorkDirs()
    energies = ['19', '27', '39', '62']
    dtypes = ['purity'] + ['contam{}'.format(i+1) for i in range(len(contams))]
    data = OrderedDict()
    sampfr = {0.4: [], 0.7: [], 1: []}
    for eidx,energy in enumerate(energies):
        ekey = ' '.join([getEnergy4Key(str(energy)), 'GeV'])
        data[ekey] = [[], [], []]
        for pidx,particle in enumerate(particles):
            infile = os.path.join(inDir, '{}_{}_sampfr.dat'.format(particle, energy))
            data_import = np.loadtxt(open(infile,'r'))
            for x in sampfr:
                sampfr[x].append(sum(data_import[data_import[:,0]<x][:,1]))
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
    sampfr_tot = {}
    for x in sampfr:
        sampfr_tot[x] = sum(sampfr[x])/len(sampfr[x])
    make_panel(
        dpt_dict = data, name = os.path.join(outDir, 'purity'),
        yr = [-0.03,1.03], xr = [0.12,1.87],
        xlabel = 'momentum, p (GeV/c)',
        ylabel = 'fraction of candidate sample',
        layout = '2x2', size = '5in,7in',
        key = ['nobox', 'at graph 0.4,0.72'], lines = dict(
            ('y={}'.format(x), 'lc {} lt 3 lw 3'.format(default_colors[-8]))
            for x in sampfr
        ), key_subplot_id = 2, gpcalls = [
            'label %d "{/Helvetica=18 %.0f%%}" at %f,0.5 rotate center textcolor %s' % (
                i+5, sampfr_tot[x]*100., x-0.04, default_colors[-8]
            ) for i,x in enumerate(sampfr_tot)
        ]
    )

def gp_nsigmael():
    inDir, outDir = getWorkDirs()
    energy, mom_ranges = '39', []
    for path in glob.glob(os.path.join(
        inDir, '{}_{}_e_*'.format(particles[0], energy)
    )):
        mom_ranges.append(path.split('_')[-1][:-4])
    data = OrderedDict()
    for ridx,r in enumerate(mom_ranges):
        if ridx > 5: break # skip last momentum bin (cover 99%)
        data[r] = [[], [], []]
        for pidx,particle in enumerate(particles):
            infile = os.path.join(inDir, '{}_{}_data_{}.dat'.format(particle, energy, r))
            data_import = np.loadtxt(open(infile,'r'))
            data_import[:,2] = 0
            props = 'lw 3 pt {} lc {}'.format(4+pidx*2, default_colors[-1])
            data[r][0].append(data_import)
            data[r][1].append(props)
            data[r][2].append(particle)
    make_panel(
        dpt_dict = data, name = os.path.join(outDir, 'nsigmael_overview'),
        yr = [1e-7,0.2], xr = [-10,12], ylog = True, lmargin = 0.06,
        xlabel = 'nsigma_{el}', layout = '2x3', size = '8.3in,7.5in',
        key = ['nobox'], gpcalls = ['bars small']
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
    #gp_purity()
    gp_nsigmael()
