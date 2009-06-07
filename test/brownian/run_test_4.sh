#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

for f in out/brownian_mc_????_00000001.nc ; do
    f1=${f/_00000001.nc/}
    f2=${f1/_mc_/_mc_size_mass_}.txt
    echo "../../extract_state_aero_size_mass 1e-10 1e-4 220 ${f1}_ $f2"
    ../../extract_state_aero_size_mass 1e-10 1e-4 220 ${f1}_ $f2
done
echo "../../numeric_average out/brownian_mc_size_mass_average.txt out/brownian_mc_size_mass_????.txt"
../../numeric_average out/brownian_mc_size_mass_average.txt out/brownian_mc_size_mass_????.txt

echo "../../extract_summary_aero_size_mass out/brownian_sect_0001.nc out/brownian_sect_size_mass.txt"
../../extract_summary_aero_size_mass out/brownian_sect_0001.nc out/brownian_sect_size_mass.txt

echo "../../numeric_diff out/brownian_mc_size_mass_average.txt out/brownian_sect_size_mass.txt 0 0.7 0 0 2 0"
../../numeric_diff out/brownian_mc_size_mass_average.txt out/brownian_sect_size_mass.txt 0 0.7 0 0 2 0
exit $?