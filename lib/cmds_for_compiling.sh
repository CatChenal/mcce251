#!/bin/bash


# 1) type [+enter]:

make -f Makefile &> make.log.txt

# 2) type [+enter]:

make -f my.Makefile >> make.log.txt

# new mcce exec file created

#Clean up error-free lines in log:
sed -i '/^gcc/d' make.log.txt

/bin/mv *.o o_files/.
