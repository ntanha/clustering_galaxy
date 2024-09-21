#!/bin/bash

sim_type='PE_PI_SN noPE_PI_SN PE_noPI_SN PE_PI_SN_3Myr'

age_bin='10Myrs 100Myrs 5Myrs 40Myrs'

for sim in $sim_type
do
        python fit_smooth.py $sim 5Myrs
done