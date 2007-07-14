! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process output files to produce number of small particles and volume
! of big particle.

program sedi_bidisperse_state_to_count

  use mod_bin_grid
  use mod_environ
  use mod_aero_data
  use mod_output_state
  use mod_aero_state
  use mod_gas_data
  use mod_gas_state
  use mod_aero_dist

  integer, parameter :: out_unit = 33   ! output unit number
  integer, parameter :: state_unit = 34 ! state file unit number
  character(len=*), parameter :: out_name = "out/sedi_bidisperse_mc_counts.d"
  character(len=*), parameter :: state_prefix &
       = "out/sedi_bidisperse_mc_state_0001_"
  integer, parameter :: n_time = 600    ! number of state files
  integer, parameter :: time_inc = 10   ! increment for state files
  
  character(len=1000) :: state_name     ! name of state file to read
  type(bin_grid_t) :: bin_grid          ! bin_grid
  type(aero_data_t) :: aero_data        ! aerosol data
  type(aero_state_t) :: aero_state      ! aerosol state
  type(gas_data_t) :: gas_data          ! gas data
  type(gas_state_t) :: gas_state        ! gas state
  type(environ) :: env                  ! environment state
  real*8 :: time                        ! current time (s)

  integer i, j, k, n_small
  real*8 v_big

  open(unit=out_unit, file=out_name)

  do i = 0,n_time,time_inc
     write(state_name, "(a,i8.8,a)") state_prefix, i, ".d"
     call inout_read_state(state_name, bin_grid, aero_data, &
          aero_state, gas_data, gas_state, env, time)

     ! if there is only one particle, assume it is big
     n_small = 0
     if (total_particles(aero_state) == 1) then
        n_small = 0
        j = 1
     else
        do j = 1,bin_grid%n_bin
           if (aero_state%n(j) > 0) then
              n_small = aero_state%n(j)
              exit
           end if
        end do
     end if

     v_big = 0d0
     do j = (j + 1),bin_grid%n_bin
        do k = 1,aero_state%n(j)
           v_big = v_big + particle_volume(aero_state%v(j)%p(k,:))
        end do
     end do

     write(*,'(a8,a14,a14)') &
          't', 'n_small', 'v_big'
     write(*,'(f8.1,e14.5,e14.5)') &
          time, dble(n_small) / aero_state%comp_vol / bin_grid%dlnr, &
          v_big / aero_state%comp_vol / bin_grid%dlnr
     write(out_unit,'(e20.10,e20.10,e20.10)') &
          time, dble(n_small) / aero_state%comp_vol / bin_grid%dlnr, &
          v_big / aero_state%comp_vol / bin_grid%dlnr
  end do

  close(out_unit)

end program sedi_bidisperse_state_to_count