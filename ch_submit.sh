#!/bin/sh
#$ -S /bin/sh
#$ -N xxxx
#$ -cwd
#$ -o logrun.log
#$ -e logerr.log

/home/catalys/mcce251/bin/mcce
# my local version implements 'extra energy titration'
# with and without YMC

# Previously used Xuyu's version:
##/home/xzhu/mcce_version/Old/mcce2.4.4_extra/mcce
