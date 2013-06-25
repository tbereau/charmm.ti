#!/usr/bin/env python

# Generate, manage, and analyze TI calculations in CHARMM.

# Tristan BEREAU. March 2013

import sys, os
import numpy as np
import math
import argparse
import time
from briz import charmm

# Parse command-line options
parser = argparse.ArgumentParser(description=
  'Generate, manage, and analyze TI calculations in CHARMM',
  epilog='Tristan BEREAU (2013)')
parser.add_argument('--ti', dest='ti', type=str, required=True, 
                  choices=['vdw', 'pc', 'pcsg', 'mtp'],
                  help='TI flag: either "vdw" or "elec"')
parser.add_argument('--chm', dest='charmm', default='charmm', type=str, 
                  help='CHARMM executable')
parser.add_argument('--tps', dest='tps', type=str, 
                  help='Topology file for solute', required=True)
parser.add_argument('--top', dest='top', type=str, action='append', 
                  help='Other topology file(s)')
parser.add_argument('--par', dest='par', type=str, action='append',
                  help='Parameter file(s) for solute', required=True)
parser.add_argument('--slu', dest='slu', type=str, 
                  help='PDB of the solute molecule', required=True)
parser.add_argument('--slv', dest='slv', type=str, 
                  help='PDB of the solvent (ignore for gas phase)')
parser.add_argument('--lpun', dest='lpun', type=str,
                  help='LPUN file of the solute molecule')
parser.add_argument('--lmb', dest='lmb', type=float, nargs=3,
                  default=[0.0, 0.1, 1.0],
                  help='lambda spacing for TI: l_min, l_space, l_max. ' +
                       'Default: 0.0 0.1 1.0')
parser.add_argument('--nst', dest='nsteps', type=int,
                  default=20000,
                  help='Number of MD integration steps for each window')
parser.add_argument('--neq', dest='nequil', type=int,
                  default=10000,
                  help='Number of MD equilibration steps for each window')
parser.add_argument('--rem', dest='remote', type=str,
                  help='Name of remote cluster')
parser.add_argument('--num', dest='numproc', type=int, default=1,
                  help='Number of procs for each job submitted to remote')
parser.add_argument('--sub', dest='submit', action='store_true',
                  help='run one cycle through collection of lambda values and exit.')
parser.add_argument('--dir', dest='subdir', type=str, default="",
                  help='set remote directory to collect files from')
parser.add_argument('--back', dest='backward', action='store_true',
                  help='run backward TI calculation (default: off -> forward)')

args = parser.parse_args()


if args.ti == 'mtp':
  if args.lpun is None:
    print "Error: 'mtp' option requires *.lpun file."
    exit(1)
if args.lmb[0] < 0.0 or args.lmb[2] > 1.0 or args.lmb[1] < 0.0 or \
  args.lmb[1] > args.lmb[2]-args.lmb[0]:
  print "Error: Lambda spacing values out of bounds."
  print "If running a backward simulation, use corresponding forward bounds with option --back."
  exit(1)
if args.nsteps < 0 or args.nequil < 0 or args.nsteps < args.nequil:
  print "Error: Number of equilibration/integration steps out of bounds"
  exit(1)

print "# Parameters: nsteps=%d; nequil=%d" % (args.nsteps, args.nequil)

# Reference folder for remote
bucketDir = "bucket"

################ FUNCTIONS ################
def die(arg):
  print "Error: "+ arg + " file hasn't been set or doesn't exist."
  exit(1)

def getTrjFile(lambdai, lambdan, lambdaf):
  return "tmp.perform.ti.charmm.script.%.6f.%.6f.%.6f.%s.trj" % (lambdai, 
    lambdan, lambdaf, args.ti)

def generateCHARMMScript(lambdai, lambdan, lambdaf, nstep, nequil, sim=True,
  analyzeLambdai=True):
  # 'sim' flag: run either MD simulation or extract energy from trajectory.
  # 'analyzeLambdai': when extracting energies, use lambda_i as coupling
  # (otherwise we use lambda_f).
  trjFile = getTrjFile(lambdai, lambdan, lambdaf)
  topSnippet = ''
  if args.ti in ['pc','mtp'] or (args.ti == 'pcsg' and sim == True):
    topSnippet = 'READ RTF CARD        NAME -\n %s\n' % \
      getScaleTopFile('%.6f' % lambdai)
  elif (args.ti == 'pcsg' and sim == False):
    lambdaEnergy = lambdaf
    if analyzeLambdai:
      lambdaEnergy = lambdai
    topSnippet = 'READ RTF CARD        NAME -\n %s\n' % \
      getScaleTopFile('%.6f' % lambdaEnergy)
  elif args.ti == 'vdw':
    topSnippet = 'READ RTF CARD        NAME -\n %s\n' % \
      getScaleTopFile('%.6f' % 0.00)
  if args.top:
    for i in range(len(args.top)):
      topFile = args.top[i]
      if args.remote:
        if topFile.find("/") != -1:
          topFile = topFile[topFile.rfind("/")+1:]
      topSnippet += 'READ RTF CARD APPEND NAME -\n %s\n' % topFile
  parSnippet = 'READ PARAM CARD        NAME -\n %s\n' % args.par[0]
  for i in range(1,len(args.par)):
    parFile = args.par[i]
    if args.remote:
      if parFile.find("/") != -1:
        parFile = parFile[parFile.rfind("/")+1:]
    parSnippet += 'READ PARAM CARD APPEND NAME -\n %s\n' % parFile
  rscaSnippet = ''
  noPSSPSnippet = 'PSSP -'
  solventSnippet1 = ''
  solventSnippet2 = ''
  pressureSnippet = '-'
  if args.slv:
    pressureSnippet = 'PCONST PINTERNAL PREFERENCE 1.0 -'
    solventSnippet1 = \
'''OPEN UNIT 10 READ CARD NAME %s
READ SEQUENCE PDB UNIT 10
GENERATE WAT NOANgle NODIhedral
REWIND UNIT 10
READ COOR PDB UNIT 10
CLOSE UNIT 10
''' % args.slv
    solventSnippet2 = \
'''SHAKE FAST WATER SELECT SEGI WAT END
'''
    if sim == True:
      solventSnippet2 += \
'''COOR STAT
CALC BOXX = ?XMAX - ?XMIN
CALC BOXY = ?YMAX - ?YMIN
CALC BOXZ = ?ZMAX - ?ZMIN
CRYSTAL DEFI ORTH @BOXX @BOXY @BOXZ 90. 90. 90. 
CRYSTAL BUILD nope 0
IMAGE BYRES XCEN 0.0 YCEN 0.0 ZCEN 0.0 SELE ALL END

NBONDS ATOM EWALD PMEWALD KAPPA 0.32  -
  FFTX 32 FFTY 32 FFTZ 32 ORDER 6 -
  CUTNB 12.0  CTOFNB 11.0 CTONnb 10.0

CONS HMCM FORCE 5.0 WEIGH REFX 0. REFY 0. REFZ 0. -
  SELECT SEGI SOLU END
'''

  dcdSnippet = 'NPRINT 1000 NSAVC -1 -'
  if args.ti == 'vdw':
    rscaSnippet = 'PERT SELE SEGI SOLU END \nSCALAR RSCA SET 0. SELE SEGI SOLU END'
    noPSSPSnippet = 'PSSP -'
  elif args.ti == 'pc':
    addTopFiles = ''
    if args.top:
      for i in range(len(args.top)):
        addTopFiles += 'READ RTF CARD APPEND NAME -\n %s\n' % args.top[i]
    rscaSnippet = \
'''PERT SELE SEGI SOLU END
DELETE ATOM SELE ALL END
READ RTF CARD NAME -
 %s
%s
OPEN UNIT 10 READ CARD NAME %s
READ SEQUENCE PDB UNIT 10
GENERATE SOLU
REWIND UNIT 10
READ COOR PDB UNIT 10
CLOSE UNIT 10

%s
''' % (getScaleTopFile('%.6f' % lambdaf), addTopFiles, args.slu, solventSnippet1)
    noPSSPSnippet = "NOPSSP"
  elif args.ti == 'pcsg':
    if sim == True:
      rscaSnippet = "OPEN WRITE UNIT 50 NAME %s" % trjFile
      dcdSnippet = 'NPRINT 100 NSAVC 100 IUNCRD 50'
    else:
      rscaSnippet = ''
      dcdSnippet = ''
      pressureSnippet = '-'
    noPSSPSnippet = "NOPSSP -"
  elif args.ti == 'mtp':
    if sim == True:
      rscaSnippet = \
'''OPEN WRITE UNIT 50 NAME %s

OPEN UNIT 40 CARD READ NAME %s
MTP MTPUNIT 40
CLOSE UNIT 40
''' % (trjFile, getScaleLpunFile('%.6f' % lambdai))
      dcdSnippet = 'NPRINT 100 NSAVC 100 IUNCRD 50'
    else:
      lambdaEnergy = lambdaf
      if analyzeLambdai:
        lambdaEnergy = lambdai
      rscaSnippet = \
'''OPEN UNIT 40 CARD READ NAME %s
MTP MTPUNIT 40
CLOSE UNIT 40
''' % getScaleLpunFile('%.6f' % lambdaEnergy)
      dcdSnippet = ''
      pressureSnippet = '-'
    noPSSPSnippet = "NOPSSP -"
  procSnippet = \
'''DYNA LEAP STRT NSTEP %s TIMESTEP 0.001 -
  NTRFRQ 100 -
  IPRFRQ 0 INBFRQ -1 IMGFRQ 250 -
  LSTART %s  LAMBDA %s  LSTOP %s  PSTART  %s -
  PSTOP  %s  PSLOW LINCR %s %s -''' % (str(nstep), str('%.6f' % lambdai), 
    str('%.6f' % lambdan), str('%.6f' % lambdaf), str(nequil), 
    str(nstep), str('%.6f' % (lambdaf-lambdai)), noPSSPSnippet)
  if args.slv:
    procSnippet = procSnippet + '''
  IHTFRQ 0 IEQFRQ 0 -
  TSTRUCT 298.0 FINALT 298.0 FIRSTT 298.0 -
  CPT PCONst PREF 1.0 PGAMMA 20.0 PMASs 500 HOOVER -
  REFT 298.9 TMASS 5000.0 -'''
  if args.ti in ['pcsg','mtp'] and sim == False:
    procSnippet = \
'''OPEN READ UNIT 50 NAME %s

CONS FIX SELE .NOT. SEGI SOLU END
TRAJECTORY FIRSTUNIT 50 NUNIT 6 SKIP 1
SET NTRAJ 0
LABEL SNAP
TRAJ READ
SET TIME ?TIME
ENERGY CUTNB 99
INCR NTRAJ BY 1
IF @NTRAJ .LT. %s GOTO SNAP
CLOSE UNIT 50
''' % (trjFile, str(int(nstep/100.)))
    dcdSnippet = ''
    pressureSnippet = ''
  returnScript = \
'''* TI calculation for file %s
* Tristan BEREAU (2013)
*

SET PRNLEV 2

%s
%s
OPEN UNIT 10 READ CARD NAME %s
READ SEQUENCE PDB UNIT 10
GENERATE SOLU
REWIND UNIT 10
READ COOR PDB UNIT 10
CLOSE UNIT 10

%s
FAST OFF
%s

ENERGY NBXMOD 5 ATOM CDIEL EPS 1.0 SHIFt VATOM VDISTANCE -
  VSWItch CUTNb 12.0 CTOFnb 11. CTONnb 10. E14Fac 1.0 

%s

%s
  %s
  %s

STOP
'''
  return returnScript % (args.slu, topSnippet, parSnippet, 
    args.slu, solventSnippet1, rscaSnippet, solventSnippet2, 
    procSnippet, pressureSnippet, dcdSnippet)

def saveToFile(data, filename):
  try:
    pipe = open(filename,'w')
    pipe.write(data)
    pipe.close()
  except IOError, e:
    raise "I/O Error",e
  if args.remote:
    rmtChm.putFile(filename)
  return

def getSubName(simCounter):
  return "ti" + str(rmtChm.subdir[3:]) + str("%05d" % simCounter)

def getSubSubDir(lambdai, lambdan, lambdaf):
  return "%.6f.%.6f.%.6f.%s" % (lambdai, lambdan, lambdaf, args.ti)

def runCHARMMScript(inpFile, outFile, simCounter, 
  dependID=0, description="", noMPI=False):
  # Run/submit CHARMM script. Returns -1*(exit code) if we run it locally or job ID
  # if we submit it.
  retvalue = 0
  if args.remote:
    # Check whether we've already submitted this script
    if inpFile in lambdaVals['submitted']:
      return 0, simCounter
    # Generate submission script (part 1: calls charmm)
    chmSubScript = rmtChm.generateCharmmJob(inpFile, outFile, 
      otherInpDir=bucketDir, email=False)
    chmSubFile   = inpFile[:inpFile.rfind(".")] + ".sub"
    saveToFile(chmSubScript, chmSubFile)
    # Generate submission script (part 2: invokes qsub)
    subName      = getSubName(simCounter)
    numProc      = args.numproc
    if noMPI:
      numProc      = 1
    retvalue = rmtChm.submitJob(subName, numProc, 
      rmtChm.remdir + "/" + rmtChm.subsubdir + "/" + chmSubFile, dependID)
    lambdaVals['submitted'].append(inpFile)
    lambdaVals['jobID'][inpFile] = retvalue
    simCounter += 1
  else:
    sysCommand = args.charmm + " < " + inpFile + " > " + outFile
    retvalue = -1*abs(int(os.system(sysCommand)))
    checkAbnormalTermination(outFile)
  return retvalue, simCounter

def checkAbnormalTermination(charmmScript):
  try:
    f = open(charmmScript,'r')
    s = f.readlines()
    f.close()
  except IOError, e:
    raise "I/O Error",e
  normalFlag = False
  energyChangeFlag = False
  for i in range(len(s)):
    if "NORMAL TERMINATION BY NORMAL STOP" in s[i]:
      normalFlag = True
    if "TOTAL ENERGY CHANGE EXCEEDED" in s[i]:
      energyChangeFlag = True
  if normalFlag is False:
    # If we get a "TOTAL ENERGY CHANGE EXCEEDED" error, then let it go. We'll
    # subdivide the lambda interval later.
    if energyChangeFlag is True:
      return
    print "Error. CHARMM script did not end properly. "
    print "Please check CHARMM output file:",charmmScript
    exit(1)
  return

def extractWindowEnergy(outFile):
  # Read TI output file. Extract energy for the current window.
  try:
    f = open(outFile,'r')
    s = f.readlines()
    f.close()
  except IOError, e:
    raise "I/O Error",e
  for i in range(len(s)):
    line = s[i].split()
    if 'PERTRES>' in line and 'EPRTOT=' in line:
      return float(line[7])
  print "Error. Can't find PERTRES statement in output file",outFile
  exit(1)


def checkTIConvergence(outFile):
  # Read TI output file. Identify DIFFLC and make sure it's below 0.5.
  # Otherwise we divide the lambda space in two.
  try:
    f = open(outFile,'r')
    s = f.readlines()
    f.close()
  except IOError, e:
    raise "I/O Error",e
  for i in range(len(s)):
    line = s[i].split()
    if 'DIFFLC' in line and len(line) == 3:
      if float(line[2]) < 0.5:
        return True
      else:
        return False
    if 'DIFFLC' in line and len(line) == 2:
      # DIFFLC is so large that it touches the equal sign:
      #   DIFFLC =*********
      return False
  print "Error. Can't find DIFFLC statement in output file",outFile
  exit(1)

def subdivideLambdaWindow(index):
  # Index corresponds to the lambdavals array. It will delete that entry and
  # replace with a finer grid.
  initVal = lambdaVals['initial'][index]
  finalVal = lambdaVals['final'][index]
  delVal = finalVal - initVal
  if abs(delVal) < 1e-6:
    print "Can't further refine lambda spacing."
    exit(1)
  lambdaVals['initial'].insert(index+1,initVal + delVal/2.)
  lambdaVals['final'].insert(index, initVal + delVal/2.)
  lambdaVals['done'][index] = False
  lambdaVals['done'].insert(index+1,False)
  lambdaVals['energy'][index] = 99999.9
  lambdaVals['energy'].insert(index, 99999.9)
  return

def scaleChargesInTop(lambda_current):
  # Scale charges and store new topology file
  try:
    f = open(args.tps,'r')
    s = f.readlines()
    f.close()
  except IOError, e:
    raise "I/O Error",e
  newTopFile = getScaleTopFile(lambda_current)
  f = open(newTopFile,'w')
  for i in range(len(s)):
    words = s[i].split()
    if len(words) > 1 and words[0] == 'ATOM':
      charge = str(float(words[3]) * math.sqrt(float(lambda_current)))
      s[i] = s[i][:s[i].find(words[3])] + charge +"\n"
    f.write(s[i])
  f.close()
  if args.remote:
    rmtChm.putFile(newTopFile)
  return

def scaleMTPInLpun(lambda_current):
  # Scale all multipole moments and store new lpun file
  try:
    f = open(args.lpun,'r')
    s = f.readlines()
    f.close()
  except IOError, e:
    raise "I/O Error",e
  newTopFile = getScaleLpunFile(lambda_current)
  f = open(newTopFile,'w')
  for i in range(len(s)):
    if s[i][0] == '#' or 'Rank' in s[i] or 'LRA:' in s[i]:
      f.write(s[i])
    else:
      words = s[i].split()
      line = ''
      for word in words:
        word = str(float(word) * math.sqrt(float(lambda_current)))
        line += word + " "
      line += "\n"
      f.write(line)
  f.close()
  if args.remote:
    rmtChm.putFile(newTopFile)
  return
 

def getScaleTopFile(lambda_current):
  return "tmp.perform.ti." + args.tps[:args.tps.rfind('.')] + \
            '.' + str(lambda_current) + \
            args.tps[args.tps.rfind('.'):]

def getScaleLpunFile(lambda_current):
  return "tmp.perform.ti." + args.lpun[:args.lpun.rfind('.')] + \
            '.' + str(lambda_current) + \
            args.lpun[args.lpun.rfind('.'):]

def extractEnergyDiff(trajFile1, trajFile2):
  trajEne = []
  for i in [1,2]:
    curFile    = str(eval('trajFile' + str(i)))
    curTrajEne = []
    try:
      f = open(curFile,'r')
      s = f.readlines()
      f.close()
    except IOError, e:
      raise "I/O Error",e  
    for i in range(len(s)):
      if "ENER EXTERN>" in s[i]:
        curTrajEne.append(float(s[i].split()[3]))
        if math.isnan(curTrajEne[-1]):
          print "Error. MTP energy is NaN."
          exit(1)
    if len(curTrajEne) == 0:
      print "Error. No energy recorded in",curFile
      exit(1)  
    trajEne.append(curTrajEne)
  if len(trajEne[0]) != len(trajEne[1]):
    print "Error. Trajectories have different numbers of energy data points."
    exit(1)
  diff = []
  difft = []
  for i in range(len(trajEne[0])):
    diff.append(trajEne[1][i] - trajEne[0][i])
    difft.append(diff[i] - diff[0])
  avg = sum(diff)/len(diff)
  avgt = sum(difft)/len(difft)
  std = math.sqrt(sum([i**2-avgt**2 for i in difft]))/len(difft)
  # Save statistics to file
  retName = trajFile1[:trajFile1.rfind(".1.")] + '.stat.dat'
  f = open(retName,'w')
  f.write("DIFFLC = " + str(std) + "\n")
  f.write("PERTRES> . . . . . EPRTOT= " + str(avg) + "\n")
  f.close()
  return retName


def runLambdaInterval(index, nstep, nequil, simCounter):
  # Run TI window given by index of the lambdaVals dictionary. Routine returns
  # True if the run was analyzed, False if we're still waiting for incoming data.
  lambda_i = lambdaVals['initial'][index]
  lambda_f = lambdaVals['final'][index]
  delVal   = lambda_f - lambda_i
  lambda_n = delVal/2. + lambda_i
  if delVal > 0. and delVal < 1.:
    if lambda_i == 0.0:
      lambda_n = 0.0
    if lambda_f == 1.0:
      lambda_n = 1.0
  elif delVal < 0. and delVal > -1.:
    if lambda_i == 1.0:
      lambda_n = 1.0
    if lambda_f == 0.0:
      lambda_n = 0.0
  # Only for args.remote -- subdirectory of the submitted simulation
  if args.remote:
    rmtChm.setSubSubDir(getSubSubDir(lambda_i, lambda_n, lambda_f))
  if args.ti in ['pc', 'pcsg', 'mtp']:
    # Scale charges
    scaleChargesInTop('%.6f' % lambda_i)
    scaleChargesInTop('%.6f' % lambda_f)
  if args.ti == 'mtp':
    # Scale MTPs
    scaleMTPInLpun('%.6f' % lambda_i)
    scaleMTPInLpun('%.6f' % lambda_f)
  print "# lambda: (%.6f - %.6f)..." % (lambda_i, lambda_f),
  sys.stdout.flush()
  if args.remote:
    print ""
  inpFiles  = []
  outFiles  = []
  trjFile   = getTrjFile(lambda_i, lambda_n, lambda_f)
  inpScript = generateCHARMMScript(lambda_i, lambda_n, lambda_f, nstep, nequil)
  inpFile   = 'tmp.perform.ti.charmm.script.'+ str('%.6f' % lambda_i)+ '.' + \
              str("%.6f" % lambda_n) + '.' + str('%.6f' % lambda_f) + '.' + \
              str(args.ti) + '.inp'
  outFile   = 'tmp.perform.ti.charmm.script.'+ str('%.6f' % lambda_i)+ '.' + \
              str("%.6f" % lambda_n) + '.' + str('%.6f' % lambda_f) + '.' + \
              str(args.ti) + '.out'
  inpFiles.append(inpFile)
  outFiles.append(outFile)
  if args.remote:
    # Leave the routine if the simulation and analysis is complete
    if outFile in lambdaVals['completed']:
      return False, simCounter
    # Don't submit new run if the file is stored locally
    if rmtChm.localFileExists(outFile) and \
      (args.ti in ['pc','vdw'] or (args.ti in ['pcsg','mtp'] and \
        rmtChm.trjfileConsistent(trjFile) ) ):
      if inpFile not in lambdaVals['submitted']:
        lambdaVals['submitted'].append(inpFile)
      jobReturn=0
    elif rmtChm.remoteFileExists(outFile) and \
      (args.ti in ['pc','vdw'] or (args.ti in ['pcsg','mtp'] and \
        rmtChm.trjfileConsistent(trjFile) ) ):
      if inpFile not in lambdaVals['submitted']:
        lambdaVals['submitted'].append(inpFile)
      rmtChm.consistentAndGet(outFile)
      if args.ti in ['pcsg','mtp']:
        rmtChm.getFile(trjFile)
      jobReturn=0
    else:
      if inpFile in lambdaVals['jobID']:
        if not rmtChm.jobIsRunning(lambdaVals['jobID'][inpFile]):
          print "delete",lambdaVals['jobID'][inpFile]," ",inpFile
          if inpFile in lambdaVals['jobID']:
            del lambdaVals['jobID'][inpFile]
          if inpFile in lambdaVals['submitted']:
            lambdaVals['submitted'].remove(inpFile)        
  saveToFile(inpScript, inpFile)
  jobReturn, simCounter = runCHARMMScript(inpFile, outFile, simCounter)
  if jobReturn < 0:
    subdivideLambdaWindow(index)
    return False, simCounter
  if args.ti in ['pcsg','mtp']:
    # Loop over lambda_i and lambda_f to extract energies
    inpReadFiles = []
    outReadFiles = []
    for i in [1,2]:
      lambdaCurIsI = True
      if i == 2:
        lambdaCurIsI = False
      # Run analysis script
      inpReadScript = generateCHARMMScript(lambda_i, lambda_n, lambda_f, nstep, 
                      nequil, sim=False, analyzeLambdai=lambdaCurIsI)
      inpReadFile   = 'tmp.perform.ti.charmm.script.'+ str('%.6f' % lambda_i)+ \
                      '.' + str("%.6f" % lambda_n) + '.' + str('%.6f' % lambda_f) \
                      + '.' + args.ti + '.read.' + str(i) + '.inp'
      outReadFile   = 'tmp.perform.ti.charmm.script.'+ str('%.6f' % lambda_i)+ \
                      '.' + str("%.6f" % lambda_n) + '.' + str('%.6f' % lambda_f) \
                      + '.' + args.ti + '.read.' + str(i) + '.out'
      inpReadFiles.append(inpReadFile)
      inpFiles.append(inpReadFile)
      outReadFiles.append(outReadFile)
      outFiles.append(outReadFile)
      saveToFile(inpReadScript, inpReadFile)
      if args.remote:
        # Don't submit new run if the file is stored locally
        if rmtChm.localFileExists(outReadFile):
          if inpReadFile not in lambdaVals['submitted']:  
            lambdaVals['submitted'].append(inpReadFile)
        else:
          # Check if the run exists on the remote server
          if rmtChm.consistentAndGet(outReadFile):
            # No need for a new simulation
            if inpReadFile not in lambdaVals['submitted']:  
              lambdaVals['submitted'].append(inpReadFile)
          else:
            # We'll need a simulation. Don't make it depend on previous
            # simulation in case the trj file is present.
            if rmtChm.remoteFileExists(trjFile):
              # Check that it's consistent: download the file and run a
              # catdcd.
              rmtChm.getFile(trjFile)
              if rmtChm.trjfileConsistent(trjFile):
                jobReturn = 0
            else:
              # Check integrity of trjfile
              if rmtChm.trjfileConsistent(trjFile):
                rmtChm.putFile(trjFile)
                jobReturn = 0
      jobTmpReturn, simCounter = runCHARMMScript(inpReadFile, outReadFile, 
        simCounter, dependID=jobReturn, noMPI=True)
    for i in [1,2]:
      if args.remote:
        allFilesFinished = True
        for j in range(len(inpFiles)):
          # Look for files locally (e.g., from a previous run)
          if rmtChm.localFileExists(outFiles[j]) == False:
            allFilesFinished = False
            # Don't analyze anything unless the files are ready
            if rmtChm.remoteFileExists(outFiles[j]) == False:
              if inpFiles[j] in lambdaVals['jobID']:
                if not rmtChm.jobIsRunning(lambdaVals['jobID'][inpFiles[j]]):
                  del lambdaVals['jobID'][inpFiles[j]]
                  lambdaVals['submitted'].remove(inpFiles[j])
            else:
              if rmtChm.consistentAndGet(outFiles[j]) is False:
                lambdaVals['submitted'].remove(inpFiles[j])
                rmtChm.delFile(outFiles[j])
              else:
                if outFiles[j] == outFile:
                  rmtChm.getFile(trjFile)
        if allFilesFinished == False:
          return False, simCounter
        lambdaVals['completed'].append(outFile)  
    tempFile = extractEnergyDiff(outReadFiles[0], outReadFiles[1])
    # outFile is now tempFile--it contains the free energy result.
    outFile = tempFile
  else:
    # vdW
    if args.remote:
      # First look for file locally, then remotely
      if rmtChm.localFileExists(outFile) == False:
        if rmtChm.remoteFileExists(outFile) == False:
          if inpFile in lambdaVals['jobID']:
            if not rmtChm.jobIsRunning(lambdaVals['jobID'][inpFile]):
              del lambdaVals['jobID'][inpFile]
              lambdaVals['submitted'].remove(inpFile)
          return False, simCounter
        else:
          if rmtChm.consistentAndGet(outFile) is False:
            lambdaVals['submitted'].remove(inpFile)
            rmtChm.delFile(outFile)
            return False, simCounter
      lambdaVals['completed'].append(outFile)
  if args.remote:
    rmtChm.delRemoteSubSubDir()      
  if checkTIConvergence(outFile):
    lambdaVals['done'][index] = True
    lambdaVals['energy'][index] = extractWindowEnergy(outFile)
    if args.ti == 'vdw':
      # TI for vdW is set up the wrong way: from lambda=1 to 0.
      lambdaVals['energy'][index] *= -1.0
    print "(dE: %9.5f) succeeded" % lambdaVals['energy'][index]
    return True, simCounter
  else:
    print "failed. Subdividing."
    subdivideLambdaWindow(index)
    return True, simCounter


def allTIDone():
  # Check whether all TI simulations have been performed successfully. Return
  # True if that's the case. Otherwise, return also the index of the first
  # simulation that still needs to be done.
  for i in range(len(lambdaVals['done'])): 
    if lambdaVals['done'][i] is False:
      return False
  # All done, return True
  return True

############################################

if args.remote:
  # Connect to the server
  rmtChm = charmm.RunCharmmRemotely(server=args.remote, 
    subdir=args.subdir)
  # Copy some files to subfolder 'bucket'
  rmtChm.setSubSubDir(bucketDir)
  rmtChm.putFile(args.tps)
  if args.top:
    for myFile in args.top:
      rmtChm.putFile(myFile)
  for myFile in args.par:
    rmtChm.putFile(myFile)
  rmtChm.putFile(args.slu)
  if args.slv:
    rmtChm.putFile(args.slv)
  if args.lpun:
    rmtChm.putFile(args.lpun)
else:
  if args.submit:
    print "Error. Can't submit without remote."
    exit(1)

# Backward simulation?
if args.backward:
  # Swap forward bounds args.lmb[0] and [2]. Invert sign of args.lmb[1].
  tmp = args.lmb[2]
  args.lmb[2]  = args.lmb[0]
  args.lmb[0]  = tmp
  args.lmb[1] *= -1.

# Initialize lambda windows
lambdaVals = {}
lambdaVals['initial'] = [args.lmb[0]]
lambdaVals['final'] = []
lambdaVals['done'] = [False]
# arbitrarily large energy
lambdaVals['energy'] = [99999.9]
lmb_cur = args.lmb[0] + args.lmb[1]
while (not args.backward and lmb_cur <= args.lmb[2] - args.lmb[1]) \
  or (args.backward and lmb_cur >= args.lmb[2] - args.lmb[1]):
  lambdaVals['initial'].append(lmb_cur)
  lambdaVals['final'].append(lmb_cur)
  lambdaVals['done'].append(False)
  # arbitrarily large energy
  lambdaVals['energy'].append(99999.9)
  lmb_cur += args.lmb[1]
lambdaVals['final'].append(args.lmb[2])
if args.remote:
  # Are jobs submitted?
  lambdaVals['submitted'] = []
  # Are jobs completed?
  lambdaVals['completed'] = []
  # ID of running jobs -- dictionary
  lambdaVals['jobID'] = {}

# MTP or vdW -> no charges in top file
if args.ti in ['mtp','vdw']:
  scaleChargesInTop('%.6f' % 0.00)

simCounter = 0
# Main loop. Calculates all lambda windows.
allDone = allTIDone()
while allDone is False:
  # Don't wait if we've gotten new analysis back.
  newAnalysis = False
  for simIndex in range(len(lambdaVals['done'])):
    if lambdaVals['done'][simIndex] is False:
      status, simCounter = runLambdaInterval(index=simIndex, nstep=args.nsteps, 
        nequil=args.nequil, simCounter=simCounter)
      if status == True:
        newAnalysis = True
  if args.submit ==  True:
    print "# Submitted all simulations."
    print "# Remote directory: %s. " % (rmtChm.getDir())
    print "# jobID: %s" % getSubName(0)[:5]
    print "# Exiting"
    exit(0)
  if newAnalysis == False:
    # Sleep for a while
    time.sleep(60)
  allDone = allTIDone()

# Print results
totalEnergy = 0.0
print "# lambda_i  lambda_f    deltaG"
for i in range(len(lambdaVals['initial'])):
  print "  %.6f - %.6f: %9.5f" % (lambdaVals['initial'][i], lambdaVals['final'][i],
    lambdaVals['energy'][i])
  totalEnergy += lambdaVals['energy'][i]
print "##############################"
print "# %.6f - %.6f: %9.5f kcal/mol" % (lambdaVals['initial'][0], 
  lambdaVals['final'][-1], totalEnergy)

# Remove remote directory
if args.remote:
  rmtChm.delRemoteSubDir()

print "# Normal termination"
