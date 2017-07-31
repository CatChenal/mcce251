#!/bin/bash

if [ $(pwd) != /home/catalys/mcce251 ]; then
  echo "Not in /home/catalys/mcce251. Check if correct before compiling local mcce."
  exit 0
fi
cd /lib; make; 
cd ../; make
