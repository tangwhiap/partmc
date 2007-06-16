run_type mc                     # Monte Carlo run
output_file out/golovin_mc_out.d # name of output file
state_prefix out/golovin_mc_state # prefix of state files
n_loop 10                       # number of Monte Carlo loops
num_conc 1d9                    # particle concentration (#/m^3)
kernel golovin                  # coagulation kernel

t_max 600                       # total simulation time (s)
del_t 1                         # timestep (s)
t_output 60                     # output interval (0 disables) (s)
t_state 0                       # state output interval (0 disables) (s)
t_progress 60                   # progress printing interval (0 disables) (s)

temp_profile temp_constant_15C.dat # temperature profile file
RH 0.999                        # initial relative humidity (1)
pressure 1d5                    # initial pressure (Pa)
rho_a 1.25                      # initial air density (kg/m^3)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

n_bin 160                       # number of bins
v_min 1d-24                     # volume of smallest bin (m^3)
scal 3                          # scale factor (integer)

gas_data gas_data_simple.dat    # file containing gas data
gas_init gas_init_simple.dat    # initial gas concentrations
gas_emissions gas_none.dat      # gas emissions file
gas_background gas_none.dat     # background gas concentrations file

aerosol_data aerosol_data_water.dat # file containing aerosol data
aerosol_init golovin_mc_init.dat # aerosol initial condition file
aerosol_emissions aerosol_none.dat # aerosol emissions file
aerosol_background aerosol_none.dat # aerosol background file

rand_init 17                    # random initialization (0 to auto-generate)
do_coagulation yes              # whether to do coagulation (yes/no)
allow_double yes                # whether to allow doubling (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic no                    # whether to do MOSAIC (yes/no)
do_restart no                   # whether to restart from stored state (yes/no)
restart_name XXXX.d             # filename to restart from
