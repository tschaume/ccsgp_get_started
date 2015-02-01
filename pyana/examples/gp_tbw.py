import os, argparse, logging, math
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot, make_panel
from ..ccsgp.config import default_colors
import numpy as np
from scipy.interpolate import interp1d
import uncertainties.unumpy as unp

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
                    props = None
                    if dtype == 'data':
                        data_import[:,3] = np.sqrt(
                            data_import[:,3]*data_import[:,3] +
                            data_import[:,4]*data_import[:,4]
                        )
                        data_import[:,(2,4)] = 0.
                        props = 'lt 1 lw 2 pt 18'
                    elif dtype == 'tbw':
                        props = 'with lines lt 1 lw 2'
                        if energy == '19':
                            # scale 19 tbw to data
                            filename19 = '_'.join(['data', code, energy]) + '.dat'
                            data19 = np.loadtxt(open(
                                os.path.join(inDir, filename19), 'rb'
                            ))
                            data19[:,3] = np.sqrt(
                                data19[:,3]*data19[:,3] + data19[:,4]*data19[:,4]
                            )
                            data19[:,(2,4)] = 0.
                            tbwF = interp1d(data_import[:,0], data_import[:,1])
                            tbwY = np.array([ tbwF(x) for x in data19[:,0] if x < 2 ])
                            data19Y = unp.uarray(data19[:,1], data19[:,3])
                            ratios = data19Y[:len(tbwY)] / tbwY
                            normfactor = np.average(
                                unp.nominal_values(ratios), weights=1./unp.std_devs(ratios)
                            )
                            print filename, normfactor
                            data_import[:,1] *= normfactor
                    data[particle][0].append(data_import)
                    data[particle][1].append(props)
                    data[particle][2].append(filename)
    make_panel(
        dpt_dict = data,
        name = os.path.join(outDir, 'tbw'),
        yr = [8e-3,300], xr = [-0.05,2.05], ylog = True,
        xlabel = 'p_{T} (GeV/c)', ylabel = 'd^{2}N/2{/Symbol \160}p_{T}dp_{T}dy',
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
