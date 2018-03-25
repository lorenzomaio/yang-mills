#!/usr/bin/env python3

import sys
import os

try:
  with open("statanalysis/jackknife.py") as f: 
    pass
except:
  print("Get the package 'statanalysis' from")
  print("https://github.com/claudio-bonati/statanalysis")
  print("and put it in this directory.")
  print()
  print("If you have 'git' installed you can use")
  print("git clone https://github.com/claudio-bonati/statanalysis")
  print()
  sys.exit(1)

sys.path.append(os.getcwd()+"/statanalysis")

import numpy as np
import os.path
import statanalysis.jackknife as jack

# main
if __name__=="__main__":

  try:
    infile =sys.argv[1]
  except:
    print("USE: %s file_name block_size\n" % sys.argv[0])
    sys.exit(1)

  try:
    blocksize =int(sys.argv[2])
  except:
    print("USE: %s file_name block_size\n" % sys.argv[0])
    sys.exit(1)

  if not os.path.isfile(infile):
    print("ERROR: file %s does not exists\n" % infile)
    sys.exit(1)

  # functions to be evaluated
  def id(x):
    return x

  def modulo(x):
    return abs(x)

  # data acquisition
  indata=np.loadtxt(infile, skiprows=1, dtype=np.float)
  data=np.transpose(indata)     #column ordered

  # primary observables
  for element in data:
    ris, err = jack.jackknife_for_primary(id, element, blocksize)
    print(ris, err, end=' ')
    ris, err = jack.jackknife_for_primary(modulo, element, blocksize)
    print(ris, err)
