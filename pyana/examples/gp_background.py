import logging, argparse, re, os
import numpy as np
from .utils import getWorkDirs, checkSymLink, getEnergy4Key
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.config import default_colors

def gp_background():
    """ plot background methods and S/B vs energy """
    inDir, outDir = getWorkDirs()
    graph_data = []
    energy, REBIN = None, None
    for infile in [ 'eps', 'ngm_corr', 'epm_normed', 'rgm', 'ac' ]:
        file_url = os.path.realpath(os.path.join(inDir, infile) + '.dat')
        if energy is None: energy = file_url.split('/')[-3]
        if os.path.isdir(file_url): continue
        data_import = np.loadtxt(open(file_url, 'rb'))
        if REBIN is None: REBIN = int(data_import[-1][2]*2*1000) # MeV
        # reset syst. uncertainties to stat. for filledcurves
        data_import[:,4] = data_import[:,3]
        data_import[:,3] = 0
        data_import[:,2] = 0 # no x errors
        if infile == 'ngm_corr':
            data_import = data_import[data_import[:,0] <= 0.9]
        if infile == 'epm_normed':
            data_import = data_import[data_import[:,0] > 0.9]
        graph_data.append(data_import)
    # unsubtracted background
    make_plot(
        name = '%s/methods%s' % (outDir, energy),
        xr = [0,3.5], yr = [0.9,2e5], ylog = True,
        data = graph_data[:3],
        properties = [
            'with filledcurves lt 1 lw 3 lc %s pt 0' % default_colors[i]
            for i in xrange(3)
        ],
        titles = [
            'SE_{+-}', 'SE@^{acc-corr\'ed}_{/Symbol \\261\\261}',
            'ME@^{norm\'ed}_{+-}'
        ],
        xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
        ylabel = 'counts / %d MeV/c^{2}' % REBIN,
        labels = {
            '%s GeV' % energy: (2.2, 5e3),
            'NR: M_{ee} > 0.9 GeV/c^{2}': (2.2, 2e3),
            'A_{/Symbol \\261} = 0.026813(25)': (2.2, 0.7e3)
        },
    )
    # background ratio and acc.corr.
    make_plot(
        name = '%s/ratios%s' % (outDir, energy),
        xr = [0,1.6], yr = [0.95,1.2],
        data = graph_data[3:],
        properties = [
            'with filledcurves lt 1 lw 3 lc %s pt 0' % default_colors[i]
            for i in xrange(2)
        ],
        titles = [
            'SE_{/Symbol \\261\\261} / ME_{/Symbol \\261\\261}',
            'f_{acc} = ME_{+-} / ME_{/Symbol \\261\\261}'
        ],
        xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
        ylabel = '', key = [ 'width -2' ],
        labels = { '%s GeV' % energy: (0.4, 0.97) }
    )
    # signal-to-background ratio in rho/omega region vs. energy
    graph_data_sn = []
    for infile in os.listdir(os.path.join(inDir, 'sn')):
        energy = re.compile('\d+').search(infile).group()
        file_url = os.path.join(inDir, 'sn', infile)
        data_import = np.loadtxt(open(file_url, 'rb'))
        mask = (data_import[:,0] > 0.3) & (data_import[:,0] < 0.75)
        data_import = data_import[mask]
        weights = 1./data_import[:,3]
        sn = np.average(data_import[:,1], weights = weights)
        sn_err = np.average((data_import[:,1]-sn)**2, weights = weights)
        graph_data_sn.append([float(getEnergy4Key(energy)), sn, 0, sn_err, 0])
    graph_data_sn = np.array(graph_data_sn)
    make_plot(
        name = '%s/SNvsEnergy' % (outDir), xr = [15,210],
        yr = [1e-3, .11], xlog = True, ylog = True,
        data = [ np.array([[15,0.1,0,0,0],[210,0.1,0,0,0]]), graph_data_sn ],
        properties = [
            'with lines lt 2 lw 4 lc 0',
            'lt 1 lw 3 lc %s pt 18 ps 2' % default_colors[0]
        ],
        titles = ['']*2,
        xlabel = '{/Symbol \326}s_{NN} (GeV)',
        ylabel = 'S/B for 0.3 < M_{ee} < 0.75 GeV/c^{2}',
        lmargin = 0.1, gpcalls = [
            'nokey', 'format x "%g"',
            'xtics (20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
        ],
        labels = { 'p+p': (100, 0.09) }
    )
    return 'done'

if __name__ == '__main__':
    checkSymLink()
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", help="show log output", action="store_true")
    args = parser.parse_args()
    loglevel = 'DEBUG' if args.log else 'WARNING'
    logging.basicConfig(
        format='%(message)s', level=getattr(logging, loglevel)
    )
    print gp_background()
