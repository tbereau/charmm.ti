#!/bin/bash

# Initialize variables
scriptDir=${0%/*}
pyScriptName=$scriptDir/perform.ti.py
charmm="charmm"
numproc=1
PARFILES=()
topfile=
TOPFILES=()
slv=solvent.shifted.pdb

function show_help 
{
  echo "Usage: $0 [-c charmm] [-n numProc] [-p file.par]"
}

function die
{
  echo $1
  exit 1
}

OPTIND=1
while getopts "h?:c:n:p:t:q:" opt; do
  case "$opt" in
    h|\?)
      show_help
      exit 0
      ;;
    c)
      charmm=$OPTARG
      ;;
    n)
      numproc=$OPTARG
      [ $numproc -lt 1 ] && die "Error in number of CPUs"
      ;;
    p)
      PARFILES+=("--par $OPTARG")
      ;;
    t)
      topfile=$OPTARG
      ;;
    q)
      TOPFILES+=("--top $OPTARG")
      ;;
  esac
done
shift $((OPTIND-1)) # Shift off the options

echo ${PARFILES[@]}
[ -z $topfile ] && die "Missing topology file for solute"
exit

# Vacuum simulations
# ~/hamm/scripts/perform.ti.py --chm ~/soft/charmm.c36b1.mtpl.w.forces --tps acetone.top --top ~/soft/c36b1.mtpl.w.forces/toppar/top_all27_prot_na.rtf --par hamm.opt.par --par ~/soft/c36b1.mtpl.w.forces/toppar/par_all27_prot_na.prm --ti pcsg --slu acetone.pdb --lpun acetone.lpun --nst 20000 --neq 10000 --rem verdi --lmb 0.0 0.02 1.0 --sub > ti.pcsg.2ns.dat| tee

# Water simulations
for sim in vdw pcsg mtp; do
  echo "Submitting $sim"
  # Submit jobs
  $pyScriptName \
    # CHARMM executable
    --chm $charmm \
    # Topology file of solute
    --tps $topfile \
    # Add list of top files
    ${TOPFILES[@]} \
    # Add list of par files
    ${PARFILES[@]} \
    --ti $sim \
    --slu $1.pdb \
    --lpun $1.lpun \
    --nst 20000 \
    --neq 10000 \
    --rem verdi \
    --lmb 0.0 0.02 1.0 \
    --slv $slv \
    --sub \
    --num $numproc > ti.$sim.2ns.dat | tee
done



