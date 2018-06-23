# Introduction

Contact Number Calculator 2.0 by Yiming Tang @ Fudan.

I hate manuals so I'm not going to talk about how this program is constructed. Nevertheless, it is useful to tell you a little about how to make full use of this program to calculate contact number between groups.

Firstly, you should use command line to select three sets of groups: (1) -all: this should always contain all atoms in your trajectory file. (2) -reference and (3) -select: these two sets of groups support multivalue-selection. The contacts are calculated between each group-pair between these two sets reference1-select1, reference1-select2, ..., reference2-select1, etc.

Secondly, you should tell the program what output you want. Three kinds of output files can be processed: (1) -map: the frame-averaged contact number/probability map in XPM file which can be further transfered into eps using gmx xpm2ps command. (3) -verbose: the frame-averaged contact number/probability between each group-pair as a function of time.

If -probability flag is specified, a Contact Probability Map is calculated instead of contact number. In Contact Probability Map, two residues are said to be in contact if there exist at least one contact between them.

NOTICE: As this program contains two group-selections that support multivalue, you should always use conmmand line option of select instead of using user-interface.

For further information contact Yiming Tang (tym@tymworld.com) but there's significant chance that I ignore you. Bite me!

# Features V2.0
- This program is totally rewritten to use neighborhood searching algorithm to speed up calculations.

# Features V2.1
- This program now supports "-probability" flag which computes a contact-probability map instead of contact-number map.


