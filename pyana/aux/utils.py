import sys, os, itertools

def getOpts(i):
  opts = 'lt 1 lw 4 ps 2 '
  if not i%2: opts += 'lc %d pt 18' % int(i/2)
  else: opts += 'lc 0'
  return opts

def checkSymLink():
  link_name = __name__.split('.')[0] + 'Dir'
  if not os.path.islink(link_name):
    logging.critical('create symlink %s to continue!' % link_name)
    sys.exit(1)

def getWorkDir(module, isInput = False):
  # use same output layout as for package
  # split module & append 'Dir' to package name
  dirs = module.split('.')
  dirs[0] += 'Dir'
  if isInput: dirs[0] += '/input'
  pth = '/'.join(dirs)
  if not isInput and not os.path.exists(pth): os.makedirs(pth)
  return pth

def zip_flat(a, b):
  return list(itertools.chain.from_iterable(zip(a, b)))

