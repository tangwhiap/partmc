run_type sect                   # sectional code
output_file out/sedi_exp_sect_out.d # name of output file
num_conc 1d9                    # particle concentration (#/m^3)
kernel sedi                     # coagulation kernel

t_max 600                       # total simulation time (s)
del_t 1                         # timestep (s)
t_output 60                     # output interval (0 disables) (s)
t_progress 60                   # progress printing interval (0 disables) (s)

n_bin 220                       # number of bins
v_min 1d-24                     # volume of smallest bin (m^3)
scal 4                          # scale factor (integer)

temp_profile temp_constant_15C.dat # temperature profile file
RH 0.999                        # initial relative humidity (1)
pressure 1d5                    # initial pressure (Pa)
rho_a 1.25                      # initial air density (kg/m^3)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

aerosol_data aerosol_data_water.dat  # file containing aerosol data

dist_type exp                   # type of initial distribution
dist_mean_vol 4.1886d-15        # mean diameter (m)
