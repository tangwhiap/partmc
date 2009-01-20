! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Read the data generated by the emission testcase and output the time
! history of the different particle sources.

program test_emission_process

  use pmc_output_state_netcdf
  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_state
  use pmc_gas_data
  use pmc_gas_state
  use pmc_env_state

  character(len=1000) :: filename
  type(bin_grid_t) :: bin_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(gas_data_t) :: gas_data
  type(gas_state_t) :: gas_state
  type(env_state_t) :: env_state
  integer :: index
  real*8 :: time
  real*8 :: del_t
  integer :: i_loop

  filename = "out/emission_mc_state_0001_00000001.nc"
  call input_state_netcdf(filename, bin_grid, aero_data, &
       aero_state, gas_data, gas_state, env_state, index, time, &
       del_t, i_loop)

  write(*,*) 'index ', index
  write(*,*) 'time ', time
  write(*,*) 'del_t ', del_t
  write(*,*) 'i_loop ', i_loop

end program test_emission_process