import sys, os, itertools, inspect, logging

def checkSymLink():
  """check for symbolic link to input directory"""
  link_name = __name__.split('.')[0] + 'Dir'
  if not os.path.islink(link_name):
    logging.critical('create symlink %s to continue!' % link_name)
    sys.exit(1)

def getWorkDirs():
  """get input/output dirs (same input/output layout as for package)"""
  # get caller module
  caller_fullurl = inspect.stack()[1][1]
  caller_relurl = os.path.relpath(caller_fullurl)
  caller_modurl = os.path.splitext(caller_relurl)[0]
  # split caller_url & append 'Dir' to package name
  dirs = caller_modurl.split('/')
  dirs[0] += 'Dir'
  # get, check and create outdir
  outDir = os.path.join(*(dirs + ['output']))
  if not os.path.exists(outDir): os.makedirs(outDir)
  # get and check indir
  dirs.append('input')
  inDir = os.path.join(*dirs)
  if not os.path.exists(inDir):
    logging.critical('create input dir %s to continue!' % inDir)
    sys.exit(1)
  return inDir, outDir
