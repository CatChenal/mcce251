The unformated pdb file test tool

This program reads an unformatted pdb file and print out the entries. It is
for testing the connection between mcce and delphi programs.

Syntax:
   testfmt unpdb

   unpdb: unformated pdb filename (such as fort.13)

Remarks:
   It has been found that delphi uses different unformated pdb file on different 
   platforms. 
   On Linux, one entry contains 5 4-byte fields that are x, y, x, c, r.

Author:
   Junjun Mao (jmao@sci.ccny.cuny.edu)
