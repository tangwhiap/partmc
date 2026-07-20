! Copyright (C) 2024 University of Illinois at Urbana-Champaign
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
! The pmc_ice_nucleation module.

module pmc_ice_nucleation
  use pmc_aero_state
  use pmc_env_state
  use pmc_aero_data
  use pmc_util
  use pmc_aero_particle
  use pmc_constants
  use pmc_rand
  use pmc_ice_nucleation_data

  implicit none

  !> Type code for an undefined or invalid immersion freezing scheme.
  integer, parameter :: IMMERSION_FREEZING_SCHEME_INVALID = 0
  !> Type code for constant ice nucleation rate (J_het) immersion freezing
  !> scheme.
  integer, parameter :: IMMERSION_FREEZING_SCHEME_CONST = 1
  !> Type code for the singular (INAS) immersion freezing scheme.
  integer, parameter :: IMMERSION_FREEZING_SCHEME_SINGULAR = 2
  !> Type code for the ABIFM immersion freezing scheme.
  integer, parameter :: IMMERSION_FREEZING_SCHEME_ABIFM = 3

  !> Type code for an undefined or invalid particle-water criterion.
  integer, parameter :: FREEZING_WATER_CRITERION_INVALID = 0
  !> Require a minimum water mass fraction.
  integer, parameter :: FREEZING_WATER_CRITERION_MASS = 1
  !> Require a minimum water-to-dry-particle volume ratio.
  integer, parameter :: FREEZING_WATER_CRITERION_VOLUME_RATIO = 2

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Main subroutine for immersion freezing simulation.
  subroutine ice_nucleation_immersion_freezing(aero_state, aero_data, &
       env_state, del_t, immersion_freezing_scheme_type, &
       freezing_rate, do_freezing_naive, INAS_a, INAS_b, ins_data, &
       water_criterion_type, water_criterion_value)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    !> Immersion freezing scheme type.
    integer, intent(in) :: immersion_freezing_scheme_type
    !> Freezing rate (only used for the constant rate scheme).
    real(kind=dp), intent(in) :: freezing_rate
    !> Whether to use the naive algorithm for time-dependent scheme.
    !> (If false, use the binned tau-leaping algorithm.)
    logical, intent(in) :: do_freezing_naive
    !> Slope parameter for the INAS parameterization (singular scheme only).
    real(kind=dp), intent(in) :: INAS_a
    !> Intercept parameter for the INAS parameterization (singular scheme only).
    real(kind=dp), intent(in) :: INAS_b
    !> Species- and scheme-specific ice nucleation data.
    type(ice_nucleation_data_t), intent(in) :: ins_data
    !> Particle-water criterion used by immersion freezing.
    integer, intent(in) :: water_criterion_type
    !> Threshold associated with the particle-water criterion.
    real(kind=dp), intent(in) :: water_criterion_value

    ! Call the immersion freezing subroutine according to the immersion
    ! freezing scheme.
    if (env_state%temp <= const%water_freeze_temp) then
       if ((immersion_freezing_scheme_type == IMMERSION_FREEZING_SCHEME_ABIFM) &
            .OR. (immersion_freezing_scheme_type == &
            IMMERSION_FREEZING_SCHEME_CONST)) then
          if (do_freezing_naive) then
             call ice_nucleation_immersion_freezing_time_dependent_naive( &
                  aero_state, aero_data, env_state, del_t, &
                  immersion_freezing_scheme_type, freezing_rate, ins_data, &
                  water_criterion_type, water_criterion_value)
          else
             call ice_nucleation_immersion_freezing_time_dependent( &
                  aero_state, aero_data, env_state, del_t, &
                  immersion_freezing_scheme_type, freezing_rate, ins_data, &
                  water_criterion_type, water_criterion_value)
          end if
       else if (immersion_freezing_scheme_type == &
            IMMERSION_FREEZING_SCHEME_SINGULAR) then
          call ice_nucleation_singular_initialize(aero_state, aero_data, &
               INAS_a, INAS_b)
          call ice_nucleation_immersion_freezing_singular(aero_state, &
               aero_data, env_state, water_criterion_type, &
               water_criterion_value)
       else
          call assert_msg(121370299, .false., &
               'Invalid immersion freezing scheme type')
       end if
    end if

  end subroutine ice_nucleation_immersion_freezing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Main subroutine for homogeneous freezing simulation.
  subroutine ice_nucleation_homogeneous_freezing(aero_state, aero_data, &
       env_state, del_t, water_criterion_type, water_criterion_value)

    type(aero_state_t), intent(inout) :: aero_state
    type(aero_data_t), intent(in) :: aero_data
    type(env_state_t), intent(inout) :: env_state
    real(kind=dp), intent(in) :: del_t
    integer, intent(in) :: water_criterion_type
    real(kind=dp), intent(in) :: water_criterion_value

    if (env_state%temp <= const%water_homo_freeze_temp) then
       call ice_nucleation_homogeneous_freezing_time_dependent_naive( &
            aero_state, aero_data, env_state, del_t, water_criterion_type, &
            water_criterion_value)
    end if

  end subroutine ice_nucleation_homogeneous_freezing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialization for the singular scheme, sampling the freezing
  !> temperature for each particle.
  subroutine ice_nucleation_singular_initialize(aero_state, aero_data, &
       INAS_a, INAS_b)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Slope parameter for the INAS parameterization.
    real(kind=dp), intent(in) :: INAS_a
    !> Intercept parameter for the INAS parameterization.
    real(kind=dp), intent(in) :: INAS_b

    integer :: i_part
    real(kind=dp) :: p, S, T0, temp
    real(kind=dp) :: aerosol_diameter

    T0 = const%water_freeze_temp
    do i_part = 1,aero_state_n_part(aero_state)
       if (aero_state%apa%particle(i_part)%imf_temperature == 0d0) then
          aerosol_diameter = aero_particle_dry_diameter( &
               aero_state%apa%particle(i_part), aero_data)
          S = const%pi * aerosol_diameter**2
          p = pmc_random()
          temp = (log(1d0 - p) + exp(-S * exp(-INAS_a * T0 + INAS_b))) / (-S)
          aero_state%apa%particle(i_part)%imf_temperature = T0 + (log(temp) &
             - INAS_b) / INAS_a
       end if
    end do

  end subroutine ice_nucleation_singular_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Simulation for singular scheme, deciding whether to freeze for each
  !> particle. Run in each time step.
  subroutine ice_nucleation_immersion_freezing_singular(aero_state, &
       aero_data, env_state, water_criterion_type, water_criterion_value)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    integer, intent(in) :: water_criterion_type
    real(kind=dp), intent(in) :: water_criterion_value
    integer :: i_part

    do i_part = 1,aero_state_n_part(aero_state)
       if (aero_state%apa%particle(i_part)%frozen) then
          cycle
       end if
       if (.not. freezing_water_criteria( &
            aero_state%apa%particle(i_part), aero_data, &
            water_criterion_type, water_criterion_value)) cycle
       if (env_state%temp <= &
            aero_state%apa%particle(i_part)%imf_temperature) then
          aero_state%apa%particle(i_part)%frozen = .true.
          aero_state%apa%particle(i_part)%den_ice = &
               const%reference_ice_density
          aero_state%apa%particle(i_part)%ice_shape_phi = 1d0
       end if
    end do

  end subroutine ice_nucleation_immersion_freezing_singular

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Simulation for time-dependent scheme (e.g., ABIFM, constant rate),
  !> deciding whether to freeze for each particle. Run in each time step.
  !> This subroutine applies the binned-tau leaping algorithm for speeding up.
  subroutine ice_nucleation_immersion_freezing_time_dependent(aero_state, &
       aero_data, env_state, del_t, immersion_freezing_scheme_type, &
       freezing_rate, ins_data, water_criterion_type, &
       water_criterion_value)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    !> Freezing rate (only used for the constant rate scheme).
    real(kind=dp), intent(in) :: freezing_rate
    !> Immersion freezing scheme type.
    integer, intent(in) :: immersion_freezing_scheme_type
    type(ice_nucleation_data_t), intent(in) :: ins_data
    integer, intent(in) :: water_criterion_type
    real(kind=dp), intent(in) :: water_criterion_value

    integer :: i_part, i_bin, i_class, n_bins, n_class
    real(kind=dp) :: a_w_ice, pis, pvs
    real(kind=dp) :: p_freeze

    real(kind=dp) :: p_freeze_max, radius_max, diameter_max

    integer :: k_th, n_parts_in_bin
    real(kind=dp) :: rand
    real(kind=dp) :: j_het_max
    integer :: rand_geo

    call aero_state_sort(aero_state, aero_data)
    pvs = env_state_saturated_vapor_pressure_wrt_water(env_state%temp)
    pis = env_state_saturated_vapor_pressure_wrt_ice(env_state%temp)
    a_w_ice = pis / pvs
    if (immersion_freezing_scheme_type == IMMERSION_FREEZING_SCHEME_ABIFM) then
       j_het_max = ice_nucleation_data_max_abifm_rate(ins_data, a_w_ice)
       if (j_het_max <= 0d0) return
    end if

    n_bins = aero_sorted_n_bin(aero_state%aero_sorted)
    n_class = aero_sorted_n_class(aero_state%aero_sorted)

    if (immersion_freezing_scheme_type == IMMERSION_FREEZING_SCHEME_CONST) then
       p_freeze_max = 1d0 - exp(-freezing_rate * del_t)
    else
       p_freeze_max = const%nan
    end if

    loop_bins: do i_bin = 1,n_bins
       loop_classes: do i_class = 1,n_class
          n_parts_in_bin = integer_varray_n_entry(&
               aero_state%aero_sorted%size_class%inverse(i_bin, i_class))
          radius_max = aero_state%aero_sorted%bin_grid%edges(i_bin + 1)
          diameter_max = radius_max * 2
          if (immersion_freezing_scheme_type == &
               IMMERSION_FREEZING_SCHEME_ABIFM) then
             p_freeze_max = ABIFM_Pfrz_max(diameter_max, j_het_max, del_t)
          end if

          if (p_freeze_max <= 0d0) cycle loop_classes

          k_th = n_parts_in_bin + 1
          loop_chosen_particles: do while(.true.)
             rand_geo = rand_geometric(p_freeze_max)
             k_th = k_th - rand_geo
             if (k_th <= 0) then
                EXIT loop_chosen_particles
             end if
             i_part = aero_state%aero_sorted%size_class &
                  %inverse(i_bin, i_class)%entry(k_th)
             if (aero_state%apa%particle(i_part)%frozen) then
                cycle
             end if
             if (.not. freezing_water_criteria( &
                  aero_state%apa%particle(i_part), aero_data, &
                  water_criterion_type, water_criterion_value)) cycle
             if (immersion_freezing_scheme_type == &
                  IMMERSION_FREEZING_SCHEME_ABIFM) then
                p_freeze = ABIFM_Pfrz_particle( &
                     aero_state%apa%particle(i_part), aero_data, ins_data, &
                     a_w_ice, del_t)
                call warn_assert_msg(301184565, p_freeze <= p_freeze_max,&
                     "p_freeze > p_freeze_max.")
                rand = pmc_random()
                if (rand < p_freeze / p_freeze_max) then
                   aero_state%apa%particle(i_part)%frozen = .true.
                   aero_state%apa%particle(i_part)%den_ice = &
                        const%reference_ice_density
                   aero_state%apa%particle(i_part)%ice_shape_phi = 1d0

                end if
             else
                aero_state%apa%particle(i_part)%frozen = .true.
                aero_state%apa%particle(i_part)%den_ice = &
                     const%reference_ice_density
                aero_state%apa%particle(i_part)%ice_shape_phi = 1d0
             end if

          end do loop_chosen_particles
       end do loop_classes
    end do loop_bins

  end subroutine ice_nucleation_immersion_freezing_time_dependent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Simulation for time-dependent scheme (e.g., ABIFM, constant rate),
  !> deciding whether to freeze for each particle. Run in each time step.
  !> This subroutine applies the naive algorithm that checks each particle.
  subroutine ice_nucleation_immersion_freezing_time_dependent_naive( &
       aero_state, aero_data, env_state, del_t, &
       immersion_freezing_scheme_type, freezing_rate, ins_data, &
       water_criterion_type, water_criterion_value)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    !> Type of the immersion freezing scheme
    integer, intent(in) :: immersion_freezing_scheme_type
    !> Freezing rate (only used for the constant rate scheme).
    real(kind=dp), intent(in) :: freezing_rate
    type(ice_nucleation_data_t), intent(in) :: ins_data
    integer, intent(in) :: water_criterion_type
    real(kind=dp), intent(in) :: water_criterion_value

    integer :: i_part
    real(kind=dp) :: a_w_ice, pis, pvs
    real(kind=dp) :: p_freeze
    real(kind=dp) :: rand

    pvs = env_state_saturated_vapor_pressure_wrt_water(env_state%temp)
    pis = env_state_saturated_vapor_pressure_wrt_ice(env_state%temp)
    a_w_ice = pis / pvs

    if (immersion_freezing_scheme_type == IMMERSION_FREEZING_SCHEME_CONST) then
       p_freeze = 1d0 - exp(-freezing_rate * del_t)
    else
       p_freeze = const%nan
    end if

    do i_part = 1,aero_state_n_part(aero_state)
       if (aero_state%apa%particle(i_part)%frozen) cycle
       if (.not. freezing_water_criteria( &
            aero_state%apa%particle(i_part), aero_data, &
            water_criterion_type, water_criterion_value)) cycle

       if (immersion_freezing_scheme_type == &
            IMMERSION_FREEZING_SCHEME_ABIFM) then
          p_freeze = ABIFM_Pfrz_particle(aero_state%apa%particle(i_part), &
               aero_data, ins_data, a_w_ice, del_t)
       end if

       if (p_freeze <= 0d0) cycle
       rand = pmc_random()
       if (rand < p_freeze) then
          aero_state%apa%particle(i_part)%frozen = .true.
          aero_state%apa%particle(i_part)%den_ice = &
               const%reference_ice_density
          aero_state%apa%particle(i_part)%ice_shape_phi = 1d0
       end if
    end do

  end subroutine ice_nucleation_immersion_freezing_time_dependent_naive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Simulates homogeneous freezing with the Koop et al. (2000) rate.
  subroutine ice_nucleation_homogeneous_freezing_time_dependent_naive( &
       aero_state, aero_data, env_state, del_t, water_criterion_type, &
       water_criterion_value)

    type(aero_state_t), intent(inout) :: aero_state
    type(aero_data_t), intent(in) :: aero_data
    type(env_state_t), intent(inout) :: env_state
    real(kind=dp), intent(in) :: del_t
    integer, intent(in) :: water_criterion_type
    real(kind=dp), intent(in) :: water_criterion_value

    integer :: i_part
    real(kind=dp) :: a_w_ice, pis, pvs, p_freeze, rand

    pvs = env_state_saturated_vapor_pressure_wrt_water(env_state%temp)
    pis = env_state_saturated_vapor_pressure_wrt_ice(env_state%temp)
    a_w_ice = pis / pvs

    do i_part = 1,aero_state_n_part(aero_state)
       if (aero_state%apa%particle(i_part)%frozen) cycle
       if (.not. freezing_water_criteria( &
            aero_state%apa%particle(i_part), aero_data, &
            water_criterion_type, water_criterion_value)) cycle
       p_freeze = Homo_Koop_Pfrz_particle( &
            aero_state%apa%particle(i_part), aero_data, a_w_ice, del_t)
       if (p_freeze <= 0d0) cycle
       rand = pmc_random()
       if (rand < p_freeze) then
          aero_state%apa%particle(i_part)%frozen = .true.
          aero_state%apa%particle(i_part)%den_ice = &
               const%reference_ice_density
          aero_state%apa%particle(i_part)%ice_shape_phi = 1d0
       end if
    end do

  end subroutine ice_nucleation_homogeneous_freezing_time_dependent_naive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Simulates melting: if the environmental temperature is above the freezing
  !> temperature of water, all particles are set to be unfrozen.
  subroutine ice_nucleation_melting(aero_state, aero_data, env_state)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state

    integer :: i_part

    if (env_state%temp > const%water_freeze_temp) then
       do i_part = 1,aero_state_n_part(aero_state)
          aero_state%apa%particle(i_part)%frozen = .false.
          aero_state%apa%particle(i_part)%den_ice = const%nan
          aero_state%apa%particle(i_part)%ice_shape_phi = const%nan
       end do
    end if

  end subroutine ice_nucleation_melting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates homogeneous freezing probability using Koop et al. (2000).
  real(kind=dp) function Homo_Koop_Pfrz_particle(aero_particle, aero_data, &
       a_w_ice, del_t)

    type(aero_particle_t), intent(in) :: aero_particle
    type(aero_data_t), intent(in) :: aero_data
    real(kind=dp), intent(in) :: a_w_ice, del_t

    real(kind=dp) :: wet_volume, dry_volume, water_volume
    real(kind=dp) :: delta_aw, log10_j_hom, j_hom

    wet_volume = aero_particle_volume(aero_particle)
    dry_volume = aero_particle_dry_volume(aero_particle, aero_data)
    water_volume = wet_volume - dry_volume
    delta_aw = 1d0 - a_w_ice
    log10_j_hom = -906.7d0 + 8502.0d0 * delta_aw &
         - 26924.0d0 * delta_aw**2 + 29180.0d0 * delta_aw**3
    j_hom = 10d0**log10_j_hom * 1d6
    Homo_Koop_Pfrz_particle = 1d0 - exp(-water_volume * j_hom * del_t)

  end function Homo_Koop_Pfrz_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Whether a particle meets a configured freezing water-content criterion.
  logical function freezing_water_criteria(aero_particle, aero_data, &
       criterion_type, criterion_value)

    type(aero_particle_t), intent(in) :: aero_particle
    type(aero_data_t), intent(in) :: aero_data
    integer, intent(in) :: criterion_type
    real(kind=dp), intent(in) :: criterion_value

    real(kind=dp) :: total_mass, water_mass
    real(kind=dp) :: wet_volume, dry_volume, water_volume

    select case (criterion_type)
    case (FREEZING_WATER_CRITERION_MASS)
       total_mass = sum(aero_particle%vol * aero_data%density)
       water_mass = aero_particle%vol(aero_data%i_water) * &
            aero_data%density(aero_data%i_water)
       if (total_mass > 0d0) then
          freezing_water_criteria = &
               water_mass / total_mass >= criterion_value
       else
          freezing_water_criteria = .false.
       end if
    case (FREEZING_WATER_CRITERION_VOLUME_RATIO)
       wet_volume = aero_particle_volume(aero_particle)
       dry_volume = aero_particle_dry_volume(aero_particle, aero_data)
       water_volume = max(0d0, wet_volume - dry_volume)
       freezing_water_criteria = water_volume >= criterion_value * dry_volume
    case default
       call assert_msg(496442384, .false., &
            'invalid freezing water criterion type')
       freezing_water_criteria = .false.
    end select

  end function freezing_water_criteria

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculating the freezing probability for the particle (i_part) using ABIFM
  !> method (Knopf et al.,2013)
  real(kind=dp) function ABIFM_Pfrz_particle(aero_particle, aero_data, &
       ins_data, a_w_ice, del_t)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Species- and scheme-specific ice nucleation data.
    type(ice_nucleation_data_t), intent(in) :: ins_data
    !> The water activity w.r.t. ice.
    real(kind=dp), intent(in) :: a_w_ice
    !> Time interval.
    real(kind=dp), intent(in) :: del_t

    real(kind=dp) :: core_volume, core_surface_area, species_surface_area
    real(kind=dp) :: abifm_m, abifm_c
    real(kind=dp) :: j_het, j_het_x_area
    integer :: i_entry, i_spec

    core_volume = 0d0
    do i_spec = 1,aero_data_n_spec(aero_data)
       if (i_spec == aero_data%i_water) cycle
       if (aero_data%is_soluble(i_spec) == 0) then
          core_volume = core_volume + aero_particle%vol(i_spec)
       end if
    end do

    if (core_volume <= 0d0) then
       ABIFM_Pfrz_particle = 0d0
       return
    end if
    core_surface_area = (36d0 * const%pi * core_volume**2)**(1d0 / 3d0)

    j_het_x_area = 0d0
    do i_entry = 1,ins_data%n_entry
       if (ins_data%scheme_type(i_entry) /= INS_SCHEME_ABIFM) cycle
       i_spec = ins_data%aero_spec(i_entry)
       if (aero_particle%vol(i_spec) <= 0d0) cycle
       abifm_m = ins_data%parameter(1, i_entry)
       abifm_c = ins_data%parameter(2, i_entry)
       j_het = 10d0**(abifm_m * (1d0 - a_w_ice) + abifm_c) * 10000d0
       species_surface_area = core_surface_area * &
            aero_particle%vol(i_spec) / core_volume
       j_het_x_area = j_het_x_area + j_het * species_surface_area
    end do

    ABIFM_Pfrz_particle = 1d0 - exp(-j_het_x_area * del_t)

  end function ABIFM_Pfrz_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculating the maximum freezing probability for particles in
  !> one bin using ABIFM method (Knopf et al.,2013). Only used by
  !> the binned-tau leaping algorithm.
  real(kind=dp) function ABIFM_Pfrz_max(diameter_max, j_het_max, &
       del_t)

    !> Maximum diameter.
    real(kind=dp), intent(in) :: diameter_max
    !> Time interval.
    real(kind=dp), intent(in) :: del_t
    !> Maximum J_het among all species.
    real(kind=dp), intent(in) :: j_het_max

    real(kind=dp) :: immersed_surface_area

    immersed_surface_area = const%pi * diameter_max**2
    ABIFM_Pfrz_max = 1d0 - exp(-j_het_max * immersed_surface_area * del_t)

  end function ABIFM_Pfrz_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for immersion freezing scheme.
  subroutine spec_file_read_immersion_freezing_scheme_type(file, &
       immersion_freezing_scheme_type)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Immersion freezing scheme type.
    integer, intent(out) :: immersion_freezing_scheme_type

    character(len=SPEC_LINE_MAX_VAR_LEN) :: imf_scheme

    !> \page input_format_imf_scheme Input File Format: Immersion Freezing Scheme
    !!
    !! The immersion freezing scheme is specified by the parameter:
    !!   - \b immersion_freezing_scheme (string): the type of
    !!   immersion freezing scheme must be one of:
    !!   \c const for freezing with a constant freezing rate;
    !!   \c singular for the INAS scheme based on singular perspective;
    !!   \c ABIFM for the ABIFM scheme based on CNT perspective.
    !!
    !! If \c immersion_freezing_scheme is \c singular, the parameters INAS_a
    !! and INAS_b need to be provided. 
    !! If \c immersion_freezing_scheme is \c const, the parameter freezing_rate
    !! needs to be provided.
    !! If \c immersion_freezing_scheme is \c ABIFM or \c const, a logical
    !! variable \c do_freezing_naive needs to be provided.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'immersion_freezing_scheme', imf_scheme)
    if (trim(imf_scheme) == 'const') then
       immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_CONST
    elseif (trim(imf_scheme) == 'singular') then
       immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_SINGULAR
    elseif (trim(imf_scheme) == 'ABIFM') then
       immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_ABIFM
    else
       call spec_file_die_msg(920761229, file, &
        "Unknown immersion freezing scheme: " // trim(imf_scheme))
    end if

  end subroutine spec_file_read_immersion_freezing_scheme_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read and validate a particle-water criterion and its threshold.
  subroutine spec_file_read_freezing_water_criterion(file, criterion_name, &
       mass_threshold_name, volume_threshold_name, criterion_type, &
       criterion_value)

    type(spec_file_t), intent(inout) :: file
    character(len=*), intent(in) :: criterion_name
    character(len=*), intent(in) :: mass_threshold_name
    character(len=*), intent(in) :: volume_threshold_name
    integer, intent(out) :: criterion_type
    real(kind=dp), intent(out) :: criterion_value

    character(len=SPEC_LINE_MAX_VAR_LEN) :: criterion

    call spec_file_read_string(file, criterion_name, criterion)
    select case (trim(criterion))
    case ('water_mass_threshold')
       criterion_type = FREEZING_WATER_CRITERION_MASS
       call spec_file_read_real(file, mass_threshold_name, criterion_value)
       if (criterion_value < 0d0 .or. criterion_value > 1d0) then
          call spec_file_die_msg(815744920, file, trim(mass_threshold_name) // &
               ' must be between 0 and 1')
       end if
    case ('water_volume_ratio_threshold')
       criterion_type = FREEZING_WATER_CRITERION_VOLUME_RATIO
       call spec_file_read_real(file, volume_threshold_name, criterion_value)
       if (criterion_value < 0d0) then
          call spec_file_die_msg(650300792, file, &
               trim(volume_threshold_name) // ' must be non-negative')
       end if
    case default
       call spec_file_die_msg(177143483, file, &
            'unknown freezing water criterion: ' // trim(criterion))
    end select

  end subroutine spec_file_read_freezing_water_criterion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_ice_nucleation
