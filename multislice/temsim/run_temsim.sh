#!/bin/bash
# run a temsim job from python
# ./run_temsim.sh -d DECK -s PROGNAME [-h] [-v] 

TEMSIM=/home/ronan/Documents/CCP4/src/multislice/kirkland/temsim

usage() {
    cat <<EOF
USAGE:
	./run_temsim.sh -d DECK -s PROGNAME -h -v
EXAMPLE:
        ./run_temsim.sh dat/test/Si110_autoslic.in autoslic  -v 
ARGUMENTS:
    -d : DECK is the full path to the input deck
    -s : PROGNAME = atompot,mulslic,autoslic
    -v : verbose
    -h : help
EOF
    exit 1;
}

## get arguments 
while getopts "hvs:d:" o; do
  case "${o}" in
      h) usage;;
      v) v=1;;
      d) deck=${OPTARG};;
      s) bin=${OPTARG};;
      *) usage;;
  esac
done

################## call temsim
cmd="cat $deck | $TEMSIM/$bin"
if [ $v ];then printf "$cmd\n";fi
#$cmd
