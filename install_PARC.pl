#!/usr/bin/perl
# install.pl #PARC 2015
#open
#
use strict;
`gfortran -c rot.f`;
`gfortran -c diag.f`;
`g++ -c PARC2015.cpp`;
`g++ -o parc2015 PARC2015.o rot.o diag.o -lgfortran`;
`g++ -c PARC_surf.cpp`;
`g++ -o parc_surfE PARC_surf.o rot.o diag.o -lgfortran`;
`g++ -c PARC_nat.cpp`;
`g++ -o parc_natE PARC_nat.o rot.o diag.o -lgfortran`;
`g++ -c parc2015trj.cpp`;
`g++ -o parc2015trj parc2015trj.o rot.o diag.o -lgfortran`;
`mkdir ene`;
`mkdir results`;
`mkdir rec`;
`mkdir trj`;
`mkdir natE`;
`mkdir surfE`;
`mkdir SS`;
print "Install Completed!\n";
