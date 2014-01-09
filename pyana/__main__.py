import logging, argparse
from examples.gp_from_txt import gp_from_txt

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("subname", help="input subdir with txt files")
parser.add_argument("--log", help="show log output", action="store_true")
args = parser.parse_args()
loglevel = 'DEBUG' if args.log else 'WARNING'
logging.basicConfig(
  format='%(message)s', level=getattr(logging, loglevel)
)

# run analysis
logging.debug( gp_from_txt(args.subname) )
