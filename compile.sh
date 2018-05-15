#!/bin/bash
rm *.o *.mod *.exe
gfortran -c module_pet.f90  -I/home/pvk03/local//netcdf-4.5.0/include -L/home/pvk03/local/netcdf-4.5.0/lib/ -lnetcdf -lnetcdff
gfortran -c penman_et.f90  -I/home/pvk03/local//netcdf-4.5.0/include -L/home/pvk03/local/netcdf-4.5.0/lib/ -lnetcdf -lnetcdff
gfortran -o pet.exe module_pet.o penman_et.o  -I/home/pvk03/local//netcdf-4.5.0/include -L/home/pvk03/local/netcdf-4.5.0/lib/ -lnetcdf -lnetcdff
exit 0
