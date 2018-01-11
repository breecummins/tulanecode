#!/bin/sh

gfortran -c density.f -o dens_f
gfortran scvt_ini.f randgen.f dens_f -o scvtini
./scvtini
cp scvt_mc.dat scvt_s.dat
gfortran scvt_opt.f process.f svtgen.f dens_f -o scvtopt
./scvtopt
cp scvt_lloyd.dat nodes.dat
gfortran draw_diag.f svtgen.f -o draw
./draw
open voronoi.eps
open deltri.eps