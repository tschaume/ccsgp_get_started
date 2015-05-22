import string, logging, argparse, os
from itertools import permutations, chain
from examples.gp_datdir import gp_datdir
from examples.gp_lcltpt import gp_lcltpt
from examples.gp_panel import gp_panel
from examples.gp_stack import gp_stack
from examples.gp_rdiff import gp_rdiff
from examples.gp_ptspec import gp_ptspec
#from examples.gp_rapp import gp_rapp # TODO: only produced for Nu Xu
#from examples.gp_xfac import gp_xfac # TODO: produce dynamically in gp_rdiff

versions = [
    'QM12', 'QM12Incl39', 'QM12InclMed', 'QM12Latest200',
    'Latest19200_PatrickQM12', 'LatestBingchuJoeyPatrickJieYi',
    'LatestPatrickJieYi'
]
ltf = [True, False]
ltf_perms = [zip(x,ltf) for x in permutations(ltf, len(ltf))]
ltf_combs = list(chain(*ltf_perms))

def getBaseDir(plot):
    # TODO: de-hardcode
    return 'data/examples/%s/input' % plot

def inBaseDirExists(plot): return os.path.exists(getBaseDir(plot))

def inDirExists(plot, version):
    return os.path.exists('%s/%s' % (getBaseDir(plot), version))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", help="show log output", action="store_true")
    args = parser.parse_args()
    loglevel = 'DEBUG' if args.log else 'WARNING'
    logging.basicConfig(
        format='%(message)s', level=getattr(logging, loglevel)
    )
    #for l in list(string.ascii_uppercase): gp_datdir(l, 4)
    #gp_lcltpt()
    for v in versions:
        if inDirExists('gp_panel', v):
            logging.info('gp_panel: %s' % v)
            gp_panel(v, None)
        if inDirExists('gp_stack', v):
            for med,fit in ltf_combs:
                logging.info('gp_stack: %s %d %d' % (v, med, fit))
                gp_stack(v, None, med, fit)
        if inDirExists('gp_rdiff', v):
            for nomed,noxerr in ltf_combs:
                logging.info('gp_rdiff: %s %d %d' % (v, nomed, noxerr))
                gp_rdiff(v, nomed, noxerr, False) # absolute
                gp_rdiff(v, nomed, noxerr, True) # relative
    if inBaseDirExists('gp_ptspec'):
        logging.info('gp_ptspec')
        gp_ptspec()
