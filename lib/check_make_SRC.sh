#!/bin/bash

outfile="dir_SRC_files"
if [[ -f $outfile ]];then
  /bin/rm $outfile
fi
touch $outfile

ls -l *.c|awk '{print $NF}'>c_files

for doc in $(cat c_files)
do
  grep $doc Makefile &> /dev/null
  if [ $? == 0 ]; then
    printf "%s \tIncluded\n" $doc >> $outfile
  else
     printf "%s \tNot included\n" $doc >> $outfile
  fi
done
echo
echo "Output (check dir files in Makefile): "$outfile
grep Not $outfile

outfile="SRC_dir_files"
if [[ -f $outfile ]];then
  /bin/rm $outfile
fi
touch $outfile

sed -n '/SRC/,/OBJ/p; /OBJ/q' Makefile | sed 's/SRC//; s/=//; s/\\/ /; s/\t/ /; /OBJ/d; s/  */\n/g' | sed 's/^$/@/g; /@/d'|sort > make_SRC_files

for src in $(cat make_SRC_files)
do
  grep $src c_files &> /dev/null
  if [ $? == 0 ]; then
    printf "%s \tIncluded\n" $src >> $outfile
  else
     printf "%s \tNot included\n" $src >> $outfile
  fi
done
echo
echo "Output (check SRC files in dir .c files): "$outfile
grep Not $outfile
echo

/bin/rm c_files make_SRC_files
