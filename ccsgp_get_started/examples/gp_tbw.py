import os, argparse, logging, math
from .utils import getWorkDirs, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot, make_panel
from ..ccsgp.config import default_colors
import numpy as np
from scipy.interpolate import interp1d
import uncertainties.unumpy as unp

def gp_tbw(shift = False):
    inDir, outDir = getWorkDirs()
    energies = ['19', '27', '39', '62']
    scale = { '19': 0.5, '27': 2, '39': 4, '62': 15 }
    particles = OrderedDict()
    particles['{/Symbol \160}'] = ['pip', 'pim']
    particles['K {/Symbol \264}3'] = ['kp', 'km']
    particles['p {/Symbol \264}5'] = ['p', 'pbar']
    dtypes = ['tbw', 'data']
    data = OrderedDict()
    for particle,codes in particles.iteritems():
        if particle not in data:
            data[particle] = [[], [], []]
        for eidx,energy in enumerate(energies):
            for cidx,code in enumerate(codes):
                if energy == '62' and code in ['p', 'pbar']: continue
                for dtype in dtypes:
                    filename = '_'.join([dtype, code, energy]) + '.dat'
                    data_import = np.loadtxt(open(
                        os.path.join(inDir, filename), 'rb'
                    ))
                    if shift: data_import[:,(1,3,4)] *= scale[energy]
                    props = None
                    if dtype == 'data':
                        data_import[:,3] = np.sqrt(
                            data_import[:,3]*data_import[:,3] +
                            data_import[:,4]*data_import[:,4]
                        )
                        data_import[:,(2,4)] = 0.
                        props = 'lt 1 lw 2 pt 18 lc %s' % default_colors[eidx]
                    elif dtype == 'tbw':
                        props = 'with lines lt %d lw 2 lc %s' % (cidx+1, default_colors[eidx])
                    if particle.startswith('K'):
                        data_import[:,(1,3)] *= 3.
                    if particle.startswith('p'):
                        data_import[:,(1,3)] *= 5.
                    data[particle][0].append(data_import)
                    data[particle][1].append(props)
                    data[particle][2].append('')
    pseudo = np.array([[-1,-1,0,0,0]])
    for cidx in range(2):
        data[particles.keys()[0]][0].append(pseudo)
        data[particles.keys()[0]][1].append('with lines lt %d lw 2 lc %s' % (cidx+1, default_colors[-1]))
        data[particles.keys()[0]][2].append('charge -1' if cidx else 'charge +1')
    for eidx,energy in enumerate(energies):
        data[particles.keys()[0]][0].append(pseudo)
        data[particles.keys()[0]][1].append('lt 1 lw 2 pt 18 lc %s ps 1.3' % default_colors[eidx])
        data[particles.keys()[0]][2].append(' '.join([ getEnergy4Key(energy), 'GeV', (
            ('{/Symbol \264}%g' % scale[energy]) if shift else ''
        )]))
    make_panel(
        dpt_dict = data,
        name = os.path.join(outDir, ('tbw_shift' if shift else 'tbw')),
        yr = [1e-2,(2e4 if shift else 500)], xr = [-0.05,2.05], ylog = True,
        xlabel = 'transverse momentum, p_{T} (GeV/c)',
        ylabel = 'd^{2}N/2{/Symbol \160}p_{T}dp_{T}dy',
        layout = '3x1', size = '3.5in,9in',
        key = ['nobox', 'samplen 0.7', 'at graph %f,1.0' % (1.1 if shift else 1.05)],
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", help="show log output", action="store_true")
    args = parser.parse_args()
    loglevel = 'DEBUG' if args.log else 'WARNING'
    logging.basicConfig(
        format='%(message)s', level=getattr(logging, loglevel)
    )
    gp_tbw()
    gp_tbw(shift = True)
