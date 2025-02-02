! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
! Copyright (C) 2009 Joseph Ching
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_condense module.

!> Water condensation onto aerosol particles.
!!
!! The model here assumes that the temperature \f$ T \f$ and pressure
!! \f$ p \f$ are prescribed as functions of time, while water content
!! per-particle and relative humidity are to be calculated by
!! integrating their rates of change.
!!
!! The state of the system is defined by the per-particle wet
!! diameters \f$ D_i \f$ and the relative humidity \f$ H \f$. The
!! state vector stores these in the order \f$ (D_1,\ldots,D_n,H)
!! \f$. The time-derivative of the state vector and the Jacobian
!! (derivative of the time-derivative with repsect to the state) all
!! conform to this ordering.
!!
!! The SUNDIALS ODE solver is used to compute the system evolution
!! using an implicit method. The system Jacobian is explicitly
!! inverted as its structure is very simple.
!!
!! All equations used in this file are written in detail in the file
!! \c doc/condense.tex.
module pmc_condense

  use pmc_aero_state
  use pmc_env_state
  use pmc_aero_data
  use pmc_util
  use pmc_aero_particle
  use pmc_constants
#ifdef PMC_USE_SUNDIALS
  use iso_c_binding
#endif

  !> Whether to numerically test the Jacobian-solve function during
  !> execution (for debugging only).
  logical, parameter :: CONDENSE_DO_TEST_JAC_SOLVE = .false.
  !> Whether to print call-counts for helper routines during execution
  !> (for debugging only).
  logical, parameter :: CONDENSE_DO_TEST_COUNTS = .false.

  !> Result code indicating successful completion.
  integer, parameter :: PMC_CONDENSE_SOLVER_SUCCESS        = 0
  !> Result code indicating failure to allocate \c y vector.
  integer, parameter :: PMC_CONDENSE_SOLVER_INIT_Y         = 1
  !> Result code indicating failure to allocate \c abstol vector.
  integer, parameter :: PMC_CONDENSE_SOLVER_INIT_ABSTOL    = 2
  !> Result code indicating failure to create the solver.
  integer, parameter :: PMC_CONDENSE_SOLVER_INIT_CVODE_MEM = 3
  !> Result code indicating failure to initialize the solver.
  integer, parameter :: PMC_CONDENSE_SOLVER_INIT_CVODE     = 4
  !> Result code indicating failure to set tolerances.
  integer, parameter :: PMC_CONDENSE_SOLVER_SVTOL          = 5
  !> Result code indicating failure to set maximum steps.
  integer, parameter :: PMC_CONDENSE_SOLVER_SET_MAX_STEPS  = 6
  !> Result code indicating failure of the solver.
  integer, parameter :: PMC_CONDENSE_SOLVER_FAIL           = 7

  integer, parameter :: CONDENSE_ICE_DEP_DENSITY_SCHEME_CHENLAMB = 1
  integer, parameter :: CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_CTL = 2
  integer, parameter :: CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_LOCAL = 3
  integer, parameter :: CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_FIX = 4
  
  integer, parameter :: CONDENSE_GAMMA_DIMENSION = 60
  integer, parameter :: CONDENSE_POKRIFKA_SI_DIMENSION = 100
  integer, parameter :: CONDENSE_POKRIFKA_DEPDEN_DIMENSION = 100

  !> Internal-use structure for storing the inputs for the
  !> rate-calculation function.
  type condense_rates_inputs_t
     !> Temperature (K).
     real(kind=dp) :: T
     !> Rate of change of temperature (K s^{-1}).
     real(kind=dp) :: Tdot
     !> Relative humidity (1).
     real(kind=dp) :: H
     !> Pressure (Pa).
     real(kind=dp) :: p
     !> Rate of change of pressure (Pa s^{-1}).
     real(kind=dp) :: pdot
     !> Computational volume (m^3).
     real(kind=dp) :: V_comp
     !> Particle diameter (m).
     real(kind=dp) :: D
     !> Particle dry diameter (m).
     real(kind=dp) :: D_dry
     !> Kappa parameter (1).
     real(kind=dp) :: kappa
  end type condense_rates_inputs_t

  !> Internal-use structure for storing the outputs from the
  !> rate-calculation function.
  type condense_rates_outputs_t
     !> Change rate of diameter (m s^{-1}).
     real(kind=dp) :: Ddot
     !> Change rate of relative humidity due to this particle (s^{-1}).
     real(kind=dp) :: Hdot_i
     !> Change rate of relative humidity due to environment changes (s^{-1}).
     real(kind=dp) :: Hdot_env
     !> Sensitivity of \c Ddot to input \c D (m s^{-1} m^{-1}).
     real(kind=dp) :: dDdot_dD
     !> Sensitivity of \c Ddot to input \c H (m s^{-1}).
     real(kind=dp) :: dDdot_dH
     !> Sensitivity of \c Hdot_i to input \c D (s^{-1} m^{-1}).
     real(kind=dp) :: dHdoti_dD
     !> Sensitivity of \c Hdot_i to input \c D (s^{-1}).
     real(kind=dp) :: dHdoti_dH
     !> Sensitivity of \c Hdot_env to input \c D (s^{-1} m^{-1}).
     real(kind=dp) :: dHdotenv_dD
     !> Sensitivity of \c Hdot_env to input \c D (s^{-1}).
     real(kind=dp) :: dHdotenv_dH
  end type condense_rates_outputs_t

  !> Internal-use structure for storing the inputs for the
  !> rate-calculation function for ice.
  type iceGrowth_rates_inputs_t
     !> Temperature (K).
     real(kind=dp) :: T
     !> Rate of change of temperature (K s^{-1}).
     real(kind=dp) :: Tdot
     !> Relative humidity (1).
     real(kind=dp) :: H
     !> Pressure (Pa).
     real(kind=dp) :: p
     !> Rate of change of pressure (Pa s^{-1}).
     real(kind=dp) :: pdot
     !> Computational volume (m^3).
     real(kind=dp) :: V_comp
     !> Particle diameter (m).
     real(kind=dp) :: D
     !> Particle dry diameter (m).
     real(kind=dp) :: D_dry
     !> Kappa parameter (1).
     !real(kind=dp) :: kappa

     !> Ice density
     real(kind=dp) :: den_ice
     real(kind=dp) :: den_dep
     real(kind=dp) :: den_sub
     !> Ice shape parameter
     real(kind=dp) :: ice_shape_phi
     real(kind=dp) :: fvi
     real(kind=dp) :: fti
     real(kind=dp) :: dfv_dR
     real(kind=dp) :: dft_dR

     !real(kind=dp) :: a_very_very_strange_bug
     !integer :: AA

  end type iceGrowth_rates_inputs_t

  !> Internal-use structure for storing the outputs from the
  !> rate-calculation function for ice.
  type iceGrowth_rates_outputs_t
     !> Change rate of diameter (m s^{-1}).
     real(kind=dp) :: Ddot
     !> Change rate of relative humidity due to this particle (s^{-1}).
     real(kind=dp) :: Hdot_i
     !> Change rate of relative humidity due to environment changes (s^{-1}).
     real(kind=dp) :: Hdot_env
     !> Sensitivity of \c Ddot to input \c D (m s^{-1} m^{-1}).
     real(kind=dp) :: dDdot_dD
     !> Sensitivity of \c Ddot to input \c H (m s^{-1}).
     real(kind=dp) :: dDdot_dH
     !> Sensitivity of \c Hdot_i to input \c D (s^{-1} m^{-1}).
     real(kind=dp) :: dHdoti_dD
     !> Sensitivity of \c Hdot_i to input \c D (s^{-1}).
     real(kind=dp) :: dHdoti_dH
     !> Sensitivity of \c Hdot_env to input \c D (s^{-1} m^{-1}).
     real(kind=dp) :: dHdotenv_dD
     !> Sensitivity of \c Hdot_env to input \c D (s^{-1}).
     real(kind=dp) :: dHdotenv_dH

     !real(kind=dp) :: AA
     !real(kind=dp) :: BB
  end type iceGrowth_rates_outputs_t

  !> Internal-use variable for storing the aerosol data during calls
  !> to the ODE solver.
  type(aero_data_t) :: condense_saved_aero_data
  !> Internal-use variable for storing the initial environment state
  !> during calls to the ODE solver.
  type(env_state_t) :: condense_saved_env_state_initial
  !> Internal-use variable for storing the rate of change of the
  !> temperature during calls to the ODE solver.
  real(kind=dp) :: condense_saved_Tdot
  !> Internal-use variable for storing the rate of change of the
  !> pressure during calls to the ODE solver.
  real(kind=dp) :: condense_saved_pdot
  !> Internal-use variable for storing the per-particle kappa values
  !> during calls to the ODE solver.
  real(kind=dp), allocatable :: condense_saved_kappa(:)
  !> Internal-use variable for storing the per-particle dry diameters
  !> during calls to the ODE solver.
  real(kind=dp), allocatable :: condense_saved_D_dry(:)
  !> Internal-use variable for storing the per-particle number
  !> concentrations during calls to the ODE solver.
  real(kind=dp), allocatable :: condense_saved_num_conc(:)

  !> Internal-use variable for storing the per-particle frozen state
  !> concentrations during calls to the ODE solver.
  !> TangWenhan
  integer, parameter :: condense_Rice_avg = 10d-6
  logical, allocatable :: condense_saved_frozen(:)
  real(kind=dp), allocatable :: condense_saved_den_ice(:)
  real(kind=dp), allocatable :: condense_saved_ice_shape_phi(:)
  real(kind=dp), allocatable :: condense_saved_ice_density_dep(:)
  real(kind=dp), allocatable :: condense_saved_ice_density_sub(:)
  real(kind=dp), allocatable :: condense_saved_fv(:)
  real(kind=dp), allocatable :: condense_saved_ft(:)
  real(kind=dp), allocatable :: condense_saved_fva(:)
  real(kind=dp), allocatable :: condense_saved_fvc(:)
  real(kind=dp), allocatable :: condense_saved_dfv_dR(:)
  real(kind=dp), allocatable :: condense_saved_dft_dR(:)
  !real(kind=dp), allocatable :: condense_saved_pokrifka_depden_dist_si(:)
  !real(kind=dp), allocatable :: condense_saved_pokrifka_depden_dist_den(:)
  !real(kind=dp), allocatable :: condense_saved_pokrifka_depden_cdf(:, :)
  !real(kind=dp), allocatable :: condense_saved_pokrifka_depden_part(:)
  !real(kind=dp), allocatable :: condense_saved_pokrifka_depden_func_part(:, :)

  !integer :: CONDENSE_POKRIFKA_SI_DIMENSION = 0
  !integer :: CONDENSE_POKRIFKA_DEPDEN_DIMENSION = 0

  logical :: condense_saved_do_ice_shape
  logical :: condense_saved_do_ice_density
  logical :: condense_saved_do_ice_ventilation
  real(kind=dp) :: condense_saved_gammaFindTrip(CONDENSE_GAMMA_DIMENSION)
  real(kind=dp) :: condense_saved_pokrifka_depden_dist_si( &
       CONDENSE_POKRIFKA_SI_DIMENSION)
  real(kind=dp) :: condense_saved_pokrifka_depden_dist_den( &
       CONDENSE_POKRIFKA_DEPDEN_DIMENSION)
  real(kind=dp) :: condense_saved_pokrifka_depden_cdf( &
       CONDENSE_POKRIFKA_SI_DIMENSION, CONDENSE_POKRIFKA_DEPDEN_DIMENSION)
  !real(kind=dp) :: condense_saved_pokrifka_depden_part( &
  !     CONDENSE_POKRIFKA_DEPDEN_DIMENSION)
  !real(kind=dp) :: condense_saved_pokrifka_depden_func_part( &
  !     CONDENSE_POKRIFKA_SI_DIMENSION, CONDENSE_POKRIFKA_DEPDEN_DIMENSION)
  logical :: condense_gamma_accessed = .false.
  logical :: condense_pokrifka_dendist_accessed = .false.
  !logical :: condense_pokrifka_depden_part_allocated = .false.
  logical :: condense_pokrifka_depden_func_part_allocated = .false.
  real(kind=dp) :: condense_saved_ice_supersat_density

  !> Internal-use variable for counting calls to the vector field
  !> subroutine.
  integer, save :: condense_count_vf
  !> Internal-use variable for counting calls to the Jacobian-solving
  !> subroutine.
  integer, save :: condense_count_solve

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do condensation to all the particles for a given time interval,
  !> including updating the environment to account for the lost
  !> water vapor.
  subroutine condense_particles(aero_state, aero_data, env_state_initial, &
       env_state_final, del_t, do_ice_shape, do_ice_density, &
       do_ice_ventilation, ice_dep_density_scheme_type)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state at the start of the timestep.
    type(env_state_t), intent(in) :: env_state_initial
    !> Environment state at the end of the timestep. The rel_humid
    !> value will be ignored and overwritten with the
    !> condensation-computed value.
    type(env_state_t), intent(inout) :: env_state_final
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    logical, intent(in) :: do_ice_shape
    logical, intent(in) :: do_ice_density
    logical, intent(in) :: do_ice_ventilation
    integer, intent(in) :: ice_dep_density_scheme_type
    real(kind=dp) :: P0

    integer :: i_part, n_eqn, i_eqn
    real(kind=dp) :: state(aero_state_n_part(aero_state) + 1)
    real(kind=dp) :: init_time, final_time
    real(kind=dp) :: abs_tol_vector(aero_state_n_part(aero_state) + 1)
    real(kind=dp) :: num_conc
    real(kind=dp) :: reweight_num_conc(aero_state_n_part(aero_state))
    real(kind=dp) :: water_vol_conc_initial, water_vol_conc_final
    real(kind=dp) :: vapor_vol_conc_initial, vapor_vol_conc_final
    real(kind=dp) :: d_water_vol_conc, d_vapor_vol_conc
    real(kind=dp) :: V_comp_ratio, water_rel_error
    !> TangWenhan
    real(kind=dp) :: ice_D3
    real(kind=dp), allocatable :: particle_volume_initial(:)
    real(kind=dp), allocatable :: particle_volume_final(:)
#ifdef PMC_USE_SUNDIALS
    real(kind=c_double), target :: state_f(aero_state_n_part(aero_state) + 1)
    real(kind=c_double), target :: abstol_f(aero_state_n_part(aero_state) + 1)
    type(c_ptr) :: state_f_p, abstol_f_p
    integer(kind=c_int) :: n_eqn_f, solver_stat
    real(kind=c_double) :: reltol_f, t_initial_f, t_final_f
#endif

#ifdef PMC_USE_SUNDIALS
#ifndef DOXYGEN_SKIP_DOC
    interface
       integer(kind=c_int) function condense_solver(neq, x_f, abstol_f, &
            reltol_f, t_initial_f, t_final_f) bind(c)
         use iso_c_binding
         integer(kind=c_int), value :: neq
         type(c_ptr), value :: x_f
         type(c_ptr), value :: abstol_f
         real(kind=c_double), value :: reltol_f
         real(kind=c_double), value :: t_initial_f
         real(kind=c_double), value :: t_final_f
       end function condense_solver
    end interface
#endif
#endif
    !print*, "start condensation"
    ! TangWenhan
    if (do_ice_shape .OR. do_ice_density .OR. do_ice_ventilation) then
        if (.NOT. condense_gamma_accessed) then
            !print*, "READ from gammaFindTrip"
            OPEN(60, file = 'gammaFindTrip')
            READ(60,*) condense_saved_gammaFindTrip
            CLOSE(60)
            condense_gamma_accessed = .true.
        end if
    end if
    if (do_ice_density .AND. ( &
         (ice_dep_density_scheme_type .eq. &
         CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_CTL) .OR. &
         (ice_dep_density_scheme_type .eq. &
         CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_LOCAL) .OR. &
         (ice_dep_density_scheme_type .eq. &
         CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_FIX) )) then
       call condense_access_pokrifka_depden_dist()
       !if ((ice_dep_density_scheme_type .eq. &
       !     CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_CTL)) then
            !.and. (.not. condense_pokrifka_depden_part_allocated)) then
       !   call condense_ice_depden_pokrifka_ctl_initialization(aero_state)
       !end if
    end if
    !print*, do_ice_shape, do_ice_density, do_ice_ventilation
    ! initial water concentration in the aerosol particles
    water_vol_conc_initial = 0d0
    do i_part = 1,aero_state_n_part(aero_state)
       num_conc = aero_weight_array_num_conc(aero_state%awa, &
            aero_state%apa%particle(i_part), aero_data)
       water_vol_conc_initial = water_vol_conc_initial &
            + aero_state%apa%particle(i_part)%vol(aero_data%i_water) * num_conc
    end do

    ! save data for use within the timestepper
    condense_saved_aero_data = aero_data
    condense_saved_env_state_initial = env_state_initial
    condense_saved_Tdot &
         = (env_state_final%temp - env_state_initial%temp) / del_t
    condense_saved_pdot &
         = (env_state_final%pressure - env_state_initial%pressure) / del_t

    !> TangWenhan
    condense_saved_do_ice_shape = do_ice_shape
    condense_saved_do_ice_density = do_ice_density
    condense_saved_do_ice_ventilation = do_ice_ventilation

    ! construct initial state vector from aero_state and env_state
    allocate(condense_saved_kappa(aero_state_n_part(aero_state)))
    allocate(condense_saved_D_dry(aero_state_n_part(aero_state)))
    allocate(condense_saved_num_conc(aero_state_n_part(aero_state)))
    !> TangWenhan
    allocate(condense_saved_frozen(aero_state_n_part(aero_state)))
    allocate(condense_saved_den_ice(aero_state_n_part(aero_state)))
    allocate(condense_saved_ice_density_dep(aero_state_n_part(aero_state)))
    allocate(condense_saved_ice_density_sub(aero_state_n_part(aero_state)))
    allocate(condense_saved_ice_shape_phi(aero_state_n_part(aero_state)))
    allocate(condense_saved_fv(aero_state_n_part(aero_state)))
    allocate(condense_saved_ft(aero_state_n_part(aero_state)))
    allocate(condense_saved_fva(aero_state_n_part(aero_state)))
    allocate(condense_saved_fvc(aero_state_n_part(aero_state)))
    allocate(condense_saved_dfv_dR(aero_state_n_part(aero_state)))
    allocate(condense_saved_dft_dR(aero_state_n_part(aero_state)))
    allocate(particle_volume_initial(aero_state_n_part(aero_state)))
    allocate(particle_volume_final(aero_state_n_part(aero_state)))

    !> TangWenhan
    call condense_ice_supersat_density(env_state_initial)

    ! work backwards for consistency with the later number
    ! concentration adjustment, which has specific ordering
    ! requirements
    do i_part = aero_state_n_part(aero_state),1,-1
       condense_saved_kappa(i_part) &
            = aero_particle_solute_kappa(aero_state%apa%particle(i_part), &
            aero_data)
       condense_saved_D_dry(i_part) = aero_data_vol2diam(aero_data, &
            aero_particle_solute_volume(aero_state%apa%particle(i_part), &
            aero_data))
       condense_saved_num_conc(i_part) &
            = aero_weight_array_num_conc(aero_state%awa, &
            aero_state%apa%particle(i_part), aero_data)

       !> TangWenhan
       condense_saved_frozen(i_part) = aero_state%apa%particle(i_part)%frozen
       condense_saved_den_ice(i_part) = aero_state%apa%particle(i_part)%den_ice
       condense_saved_ice_shape_phi(i_part) = &
           aero_state%apa%particle(i_part)%ice_shape_phi

       condense_saved_fv(i_part) = 1d0
       condense_saved_ft(i_part) = 1d0
       condense_saved_fva(i_part) = 1d0
       condense_saved_fvc(i_part) = 1d0

       state(i_part) = aero_particle_diameter(&
            aero_state%apa%particle(i_part), aero_data)

       !if (i_part .eq. 1) then
       !    print*, "state(1) = ", state(i_part)
       !end if
       !> TangWenhan
       !print*, i_part, aero_state%apa%particle(i_part)%frozen
       if (aero_state%apa%particle(i_part)%frozen) then
            
            ice_D3 = state(i_part) **3d0 + 6d0/const%pi * &
                aero_state%apa%particle(i_part)%vol(aero_data%i_water) * &
                    (const%water_density / aero_state%apa%particle(i_part)%den_ice - 1d0)
            !if (i_part .eq. 1) then
            !    print*, "state(1)=", state(i_part), aero_state%apa%particle(i_part)%vol(aero_data%i_water) * &
            !        (const%water_density), &
            !        aero_state%apa%particle(i_part)%den_ice, &
            !        aero_state%apa%particle(i_part)%frozen
            !end if
            state(i_part) = ice_D3 ** (1d0/3d0)
            
            particle_volume_initial(i_part) = const%pi / 6d0 * ice_D3
            if (do_ice_density) then
                call condense_ice_density_dep_sub( &
                    particle_volume_initial(i_part), &
                    aero_state%apa%particle(i_part)%den_ice, &
                    env_state_initial, &
                    condense_saved_ice_density_dep(i_part), &
                    condense_saved_ice_density_sub(i_part), &
                    condense_saved_do_ice_shape, &
                    ice_dep_density_scheme_type, &
                    aero_state%apa%particle(i_part))
                                
            end if
            !if (i_part .eq. 1143) then
                !print*, "start: ", state(i_part), i_part
            !    continue
            !end if
       end if

       abs_tol_vector(i_part) = max(1d-30, &
            1d-8 * (state(i_part) - condense_saved_D_dry(i_part)))
    end do

    if (do_ice_ventilation) then
        call condense_ice_ventilation(state(1:aero_state_n_part(aero_state)) / 2d0, &
                aero_state_n_part(aero_state), env_state_final)
    end if
    !print*, condense_saved_fv(1), condense_saved_ft(1), condense_saved_fva(1), &
    !            condense_saved_fvc(1)


    state(aero_state_n_part(aero_state) + 1) = env_state_initial%rel_humid
    abs_tol_vector(aero_state_n_part(aero_state) + 1) = 1d-10


    !print*, "init = ", state(1), state(aero_state_n_part(aero_state) + 1)

    !print*, "delta rou"
    !print*, condense_saved_ice_supersat_density
    !print*, "dep"
    !print*, condense_saved_ice_density_dep(1)
    !print*, "sub"
    !print*, condense_saved_ice_density_sub(1)

    !print*, "R", state(1) * 1e6, "delr", condense_saved_ice_supersat_density, "dep", condense_saved_ice_density_dep(1),&
    !   "sub", condense_saved_ice_density_sub(1), "den_ice", condense_saved_den_ice(1)

    


#ifdef PMC_USE_SUNDIALS
    ! call SUNDIALS solver
    n_eqn = aero_state_n_part(aero_state) + 1
    n_eqn_f = int(n_eqn, kind=c_int)
    reltol_f = real(1d-8, kind=c_double)
    t_initial_f = real(0, kind=c_double)
    t_final_f = real(del_t, kind=c_double)
    do i_eqn = 1,n_eqn
       state_f(i_eqn) = real(state(i_eqn), kind=c_double)
       abstol_f(i_eqn) = real(abs_tol_vector(i_eqn), kind=c_double)
    end do
    state_f_p = c_loc(state_f)
    abstol_f_p = c_loc(abstol_f)
    condense_count_vf = 0
    condense_count_solve = 0
    solver_stat = condense_solver(n_eqn_f, state_f_p, abstol_f_p, reltol_f, &
         t_initial_f, t_final_f)
    call condense_check_solve(solver_stat)
    if (CONDENSE_DO_TEST_COUNTS) then
       write(0,*) 'condense_count_vf ', condense_count_vf
       write(0,*) 'condense_count_solve ', condense_count_solve
    end if
    do i_eqn = 1,n_eqn
       state(i_eqn) = real(state_f(i_eqn), kind=dp)
    end do
#endif


    ! unpack result state vector into env_state_final
    env_state_final%rel_humid = state(aero_state_n_part(aero_state) + 1)
    !print*, "final rel_humid = ", env_state_final%rel_humid

    ! unpack result state vector into aero_state, compute the final
    ! water volume concentration, and adjust particle number to
    ! account for number concentration changes
    water_vol_conc_final = 0d0
    call aero_state_num_conc_for_reweight(aero_state, aero_data, &
         reweight_num_conc)
    do i_part = 1,aero_state_n_part(aero_state)
       num_conc = aero_weight_array_num_conc(aero_state%awa, &
            aero_state%apa%particle(i_part), aero_data)

       ! translate output back to particle
       aero_state%apa%particle(i_part)%vol(aero_data%i_water) &
            = aero_data_diam2vol(aero_data, state(i_part)) &
            - aero_particle_solute_volume(aero_state%apa%particle(i_part), &
            aero_data)


       !if (i_part .eq. 1143) then
           !print*, "end: ", state(i_part)
       !    continue
       !end if
       ! ensure volumes stay positive
       aero_state%apa%particle(i_part)%vol(aero_data%i_water) = max(0d0, &
            aero_state%apa%particle(i_part)%vol(aero_data%i_water))

       particle_volume_final(i_part) = max(&
               aero_particle_solute_volume(aero_state%apa%particle(i_part), aero_data),&
                aero_data_diam2vol(aero_data, state(i_part)))
       !print*, particle_volume_initial(i_part), particle_volume_final(i_part),&
       !    aero_state%apa%particle(i_part)%vol(aero_data%i_water)

       ! convert ice real volume to mass equivalent volume wrt water
       ! TangWenhan
       if (aero_state%apa%particle(i_part)%frozen) then
            !print*, "frozen2"
            if (do_ice_density) then
                aero_state%apa%particle(i_part)%den_ice = &
                        condense_update_ice_density(&
                            particle_volume_initial(i_part), &
                            particle_volume_final(i_part), &
                            aero_state%apa%particle(i_part)%den_ice, &
                            condense_saved_ice_density_dep(i_part), &
                            condense_saved_ice_density_sub(i_part) )
            end if
            aero_state%apa%particle(i_part)%vol(aero_data%i_water) = &
                aero_state%apa%particle(i_part)%vol(aero_data%i_water) * &
                    aero_state%apa%particle(i_part)%den_ice / &
                        const%water_density
       end if
       ! add up total water volume, using old number concentrations
       water_vol_conc_final = water_vol_conc_final &
            + aero_state%apa%particle(i_part)%vol(aero_data%i_water) * num_conc
    end do
    ! TangWenhan
    if (do_ice_shape) then
        do i_part = 1,aero_state_n_part(aero_state)
            if (aero_state%apa%particle(i_part)%frozen) then
                aero_state%apa%particle(i_part)%ice_shape_phi = &
                    condense_update_ice_shape_phi(&
                        aero_state%apa%particle(i_part)%ice_shape_phi, &
                        particle_volume_initial(i_part), &
                        particle_volume_final(i_part), &
                        env_state_initial, &
                        condense_saved_fva(i_part), &
                        condense_saved_fvc(i_part))
            end if
        end do
    end if
    !> TangWenhan
    !if (do_ice_density) then
    !    do i_part = 1,aero_state_n_part(aero_state)
    !        if (aero_state%apa%particle(i_part)%frozen) then
    !            aero_state%apa%particle(i_part)%den_ice = &
    !                condense_update_ice_density(&
    !                    particle_volume_initial(i_part), &
    !                    particle_volume_final(i_part), &
    !                    aero_state%apa%particle(i_part)%den_ice, &
    !                    condense_saved_ice_density_dep(i_part), &
    !                    condense_saved_ice_density_sub(i_part) )
    !        end if
    !    end do
    !end if
    ! adjust particles to account for weight changes
    call aero_state_reweight(aero_state, aero_data, reweight_num_conc)

    ! Check that water removed from particles equals water added to
    ! vapor. Note that water concentration is not conserved (due to
    ! computational volume changes), and we need to consider particle
    ! weightings correctly.
    V_comp_ratio = env_state_final%temp * env_state_initial%pressure &
         / (env_state_initial%temp * env_state_final%pressure)
    !print*, water_vol_conc_final * const%water_density
    P0 = env_state_saturated_vapor_pressure_water(env_state_initial%temp)
    vapor_vol_conc_initial = aero_data%molec_weight(aero_data%i_water) &
         / (const%univ_gas_const * env_state_initial%temp) &
         * env_state_sat_vapor_pressure(env_state_initial) &
         !* P0 &
         * env_state_initial%rel_humid &
         / aero_particle_water_density(aero_data)
    P0 = env_state_saturated_vapor_pressure_water(env_state_final%temp)
    vapor_vol_conc_final = aero_data%molec_weight(aero_data%i_water) &
         / (const%univ_gas_const * env_state_final%temp) &
         * env_state_sat_vapor_pressure(env_state_final) &
         !* P0 &
         * env_state_final%rel_humid &
         * V_comp_ratio / aero_particle_water_density(aero_data)

    !print*, water_vol_conc_final * const%water_density + vapor_vol_conc_final * const%water_density

    d_vapor_vol_conc = vapor_vol_conc_final - vapor_vol_conc_initial
    d_water_vol_conc = water_vol_conc_final - water_vol_conc_initial
    water_rel_error = (d_vapor_vol_conc + d_water_vol_conc) &
         / (vapor_vol_conc_final + water_vol_conc_final)
    !print*, water_rel_error
    call warn_assert_msg(477865387, abs(water_rel_error) < 1d-6, &
         "condensation water imbalance too high: " &
         // trim(real_to_string(water_rel_error)))

    

    deallocate(condense_saved_kappa)
    deallocate(condense_saved_D_dry)
    deallocate(condense_saved_num_conc)
    deallocate(condense_saved_frozen)
    deallocate(condense_saved_den_ice)
    deallocate(condense_saved_ice_density_dep)
    deallocate(condense_saved_ice_density_sub)
    deallocate(condense_saved_ice_shape_phi)
    deallocate(condense_saved_fv)
    deallocate(condense_saved_ft)
    deallocate(condense_saved_fva)
    deallocate(condense_saved_fvc)
    deallocate(condense_saved_dfv_dR)
    deallocate(condense_saved_dft_dR)
    deallocate(particle_volume_initial)
    deallocate(particle_volume_final)
    !if (condense_pokrifka_dendist_accessed) then
    !   deallocate(condense_saved_pokrifka_depden_dist_si)
    !   deallocate(condense_saved_pokrifka_depden_dist_den)
    !   deallocate(condense_saved_pokrifka_depden_cdf)
    !end if
    !if (condense_pokrifka_depden_part_allocated) then
    !   deallocate(condense_saved_pokrifka_depden_part)
    !end if
    !if (condense_pokrifka_depden_func_part_allocated) then
    !   deallocate(condense_saved_pokrifka_depden_func_part)
    !end if

    !print*, "end condensation"
  end subroutine condense_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_SUNDIALS
  !> Check the return code from the condense_solver() function.
  subroutine condense_check_solve(value)

    !> Return code to check.
    integer(kind=c_int), intent(in) :: value

    if (value == PMC_CONDENSE_SOLVER_SUCCESS) then
       return
    elseif (value == PMC_CONDENSE_SOLVER_INIT_Y) then
       call die_msg(123749472, "condense_solver: " &
            // "failed to allocate y vector")
    elseif (value == PMC_CONDENSE_SOLVER_INIT_ABSTOL) then
       call die_msg(563665949, "condense_solver: " &
            // "failed to allocate abstol vector")
    elseif (value == PMC_CONDENSE_SOLVER_INIT_CVODE_MEM) then
       call die_msg(700541443, "condense_solver: " &
            // "failed to create the solver")
    elseif (value == PMC_CONDENSE_SOLVER_INIT_CVODE) then
       call die_msg(297559183, "condense_solver: " &
            // "failure to initialize the solver")
    elseif (value == PMC_CONDENSE_SOLVER_SVTOL) then
       call die_msg(848342417, "condense_solver: " &
            // "failed to set tolerances")
    elseif (value == PMC_CONDENSE_SOLVER_SET_MAX_STEPS) then
       call die_msg(275591501, "condense_solver: " &
            // "failed to set maximum steps")
    elseif (value == PMC_CONDENSE_SOLVER_FAIL) then
       call die_msg(862254233, "condense_solver: solver failed")
    else
       call die_msg(635697577, "condense_solver: unknown return code: " &
            // trim(integer_to_string(value)))
    end if

  end subroutine condense_check_solve
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the rate of change of particle diameter and relative
  !> humidity for a single particle, together with the derivatives of
  !> the rates with respect to the input variables.
  subroutine condense_rates(inputs, outputs)

    !> Inputs to rates.
    type(condense_rates_inputs_t), intent(in) :: inputs
    !> Outputs rates.
    type(condense_rates_outputs_t), intent(out) :: outputs

    real(kind=dp) :: rho_w, M_w, P_0, dP0_dT_div_P0, rho_air, k_a, D_v, U
    real(kind=dp) :: V, W, X, Y, Z, k_ap, dkap_dD, D_vp, dDvp_dD
    real(kind=dp) :: a_w, daw_dD, delta_star, h, dh_ddelta, dh_dD
    real(kind=dp) :: dh_dH, ddeltastar_dD, ddeltastar_dH
    integer :: newton_step

    !print*, "ee", inputs%H, inputs%V_comp, inputs%D, inputs%D_dry, inputs%kappa
    rho_w = const%water_density
    M_w = const%water_molec_weight
    P_0 = const%water_eq_vap_press &
         * 10d0**(7.45d0 * (inputs%T - const%water_freeze_temp) &
         / (inputs%T - 38d0))
    dP0_dT_div_P0 = 7.45d0 * log(10d0) * (const%water_freeze_temp - 38d0) &
         / (inputs%T - 38d0)**2
    rho_air = const%air_molec_weight * inputs%p &
         / (const%univ_gas_const * inputs%T)

    k_a = 1d-3 * (4.39d0 + 0.071d0 * inputs%T)
    D_v = 0.211d-4 / (inputs%p / const%air_std_press) &
         * (inputs%T / 273d0)**1.94d0
    U = const%water_latent_heat * rho_w / (4d0 * inputs%T)
    V = 4d0 * M_w * P_0 / (rho_w * const%univ_gas_const * inputs%T)
    W = const%water_latent_heat * M_w / (const%univ_gas_const * inputs%T)
    X = 4d0 * M_w * const%water_surf_eng &
         / (const%univ_gas_const * inputs%T * rho_w)
    Y = 2d0 * k_a / (const%accom_coeff * rho_air &
         * const%air_spec_heat) &
         * sqrt(2d0 * const%pi * const%air_molec_weight &
         / (const%univ_gas_const * inputs%T))
    Z = 2d0 * D_v / const%accom_coeff * sqrt(2d0 * const%pi * M_w &
         / (const%univ_gas_const * inputs%T))

    outputs%Hdot_env = - dP0_dT_div_P0 * inputs%Tdot * inputs%H &
         + inputs%H * inputs%pdot / inputs%p
    outputs%dHdotenv_dD = 0d0
    outputs%dHdotenv_dH = - dP0_dT_div_P0 * inputs%Tdot &
         + inputs%pdot / inputs%p

    if (inputs%D <= inputs%D_dry) then
       k_ap = k_a / (1d0 + Y / inputs%D_dry)
       dkap_dD = 0d0
       D_vp = D_v / (1d0 + Z / inputs%D_dry)
       dDvp_dD = 0d0
       a_w = 0d0
       daw_dD = 0d0

       delta_star = U * V * D_vp * inputs%H / k_ap

       outputs%Ddot = k_ap * delta_star / (U * inputs%D_dry)
       outputs%Hdot_i = - 2d0 * const%pi / (V * inputs%V_comp) &
            * inputs%D_dry**2 * outputs%Ddot

       dh_ddelta = k_ap
       dh_dD = 0d0
       dh_dH = - U * V * D_vp

       ddeltastar_dD = - dh_dD / dh_ddelta
       ddeltastar_dH = - dh_dH / dh_ddelta

       outputs%dDdot_dD = 0d0
       outputs%dDdot_dH = k_ap / (U * inputs%D_dry) * ddeltastar_dH
       outputs%dHdoti_dD = - 2d0 * const%pi / (V * inputs%V_comp) &
            * inputs%D_dry**2 * outputs%dDdot_dD
       outputs%dHdoti_dH = - 2d0 * const%pi / (V * inputs%V_comp) &
            * inputs%D_dry**2 * outputs%dDdot_dH

       return
    end if

    k_ap = k_a / (1d0 + Y / inputs%D)
    dkap_dD = k_a * Y / (inputs%D + Y)**2
    D_vp = D_v / (1d0 + Z / inputs%D)
    dDvp_dD = D_v * Z / (inputs%D + Z)**2
    a_w = (inputs%D**3 - inputs%D_dry**3) &
         / (inputs%D**3 + (inputs%kappa - 1d0) * inputs%D_dry**3)
    daw_dD = 3d0 * inputs%D**2 * inputs%kappa * inputs%D_dry**3 &
         / (inputs%D**3 + (inputs%kappa - 1d0) * inputs%D_dry**3)**2

    delta_star = 0d0
    h = 0d0
    dh_ddelta = 1d0
    do newton_step = 1,5
       ! update delta_star first so when the newton loop ends we have
       ! h and dh_ddelta evaluated at the final delta_star value
       delta_star = delta_star - h / dh_ddelta
       h = k_ap * delta_star - U * V * D_vp &
            * (inputs%H - a_w / (1d0 + delta_star) &
            * exp(W * delta_star / (1d0 + delta_star) &
            + (X / inputs%D) / (1d0 + delta_star)))
       dh_ddelta = &
            k_ap - U * V * D_vp * a_w / (1d0 + delta_star)**2 &
            * (1d0 - W / (1d0 + delta_star) &
            + (X / inputs%D) / (1d0 + delta_star)) &
            * exp(W * delta_star / (1d0 + delta_star) &
            + (X / inputs%D) / (1d0 + delta_star))
    end do
    !call warn_assert_msg(387362320, &
    !     abs(h) < 1d3 * epsilon(1d0) * abs(U * V * D_vp * inputs%H), &
    !     "condensation newton loop did not satisfy convergence tolerance")

    !print*, abs(h) - 1d3 * epsilon(1d0) * abs(U * V * D_vp * inputs%H)
    !if (.not. (abs(h) < 1d3 * epsilon(1d0) * abs(U * V * D_vp * inputs%H))) then
       ! print*, abs(h)
       !print*, 1d3 * epsilon(1d0) * abs(U * V * D_vp * inputs%H)
        !print*, inputs%H, inputs%V_comp, inputs%D, inputs%D_dry, inputs%kappa
    !end if
    outputs%Ddot = k_ap * delta_star / (U * inputs%D)
    outputs%Hdot_i = - 2d0 * const%pi / (V * inputs%V_comp) &
         * inputs%D**2 * outputs%Ddot

    dh_dD = dkap_dD * delta_star &
         - U * V * dDvp_dD * inputs%H + U * V &
         * (a_w * dDvp_dD + D_vp * daw_dD &
         - D_vp * a_w * (X / inputs%D**2) / (1d0 + delta_star)) &
         * (1d0 / (1d0 + delta_star)) &
         * exp((W * delta_star) / (1d0 + delta_star) &
         + (X / inputs%D) / (1d0 + delta_star))
    dh_dH = - U * V * D_vp

    ddeltastar_dD = - dh_dD / dh_ddelta
    ddeltastar_dH = - dh_dH / dh_ddelta

    outputs%dDdot_dD = dkap_dD * delta_star / (U * inputs%D) &
         + k_ap * ddeltastar_dD / (U * inputs%D) &
         - k_ap * delta_star / (U * inputs%D**2)
    outputs%dDdot_dH = k_ap / (U * inputs%D) * ddeltastar_dH
    outputs%dHdoti_dD = - 2d0 * const%pi / (V * inputs%V_comp) &
         * (2d0 * inputs%D * outputs%Ddot + inputs%D**2 * outputs%dDdot_dD)
    outputs%dHdoti_dH = - 2d0 * const%pi / (V * inputs%V_comp) &
         * inputs%D**2 * outputs%dDdot_dH

  end subroutine condense_rates

  !> TangWenhan
  !> Compute the rate of change of ice particle diameter and relative
  !> humidity for a single ice particle, together with the derivatives of
  !> the rates with respect to the input variables.
  subroutine iceGrowth_rates(inputs, outputs)

    implicit none
    !> Inputs to rates.
    type(iceGrowth_rates_inputs_t), intent(in) :: inputs
    !> Outputs rates.
    type(iceGrowth_rates_outputs_t), intent(out) :: outputs

    real(kind = dp) :: R, Rdot, dRdot_dR, dRdot_dH, dHdoti_dR, Hdot_i, dHdoti_dH
    real(kind = dp) :: k_a, D_v, dP0_dT_div_P0, P0, P0_ice, G, rho_air, Cap, phi
    real(kind = dp) :: IWCi_dot, MRii_dot, Ls, Rd, Rv, d_Cap_d_R, dG_dR
    real(kind = dp) :: fvi, fti
    real(kind = dp) :: Lv = 2.5d6, Cpv = 1850d0, Cpw = 4200d0
    real(kind = dp) :: spec_den
    !air_spec_heat

    !print*, "Here"

    dP0_dT_div_P0 = 7.45d0 * log(10d0) * (const%water_freeze_temp - 38d0) &
        / (inputs%T - 38d0)**2
    outputs%Hdot_env = - dP0_dT_div_P0 * inputs%Tdot * inputs%H &
        + inputs%H * inputs%pdot / inputs%p
    outputs%dHdotenv_dD = 0d0
    outputs%dHdotenv_dH = - dP0_dT_div_P0 * inputs%Tdot &
        + inputs%pdot / inputs%p
    !print*, "inner = ", outputs%dHdotenv_dH
    

    !if (inputs%D - inputs%D_dry .lt. 1e-10) then
    if (inputs%D .lt. 1e-10) then
        outputs%Ddot = 0d0
        outputs%Hdot_i = 0d0
        outputs%dDdot_dD = 0d0
        outputs%dDdot_dH = 0d0
        outputs%dHdoti_dD = 0d0
        outputs%dHdoti_dH = 0d0
    else
        R = inputs%D / 2d0
        phi = inputs%ice_shape_phi
        fvi = inputs%fvi
        fti = inputs%fti
        !k_a = 1d-3 * (4.39d0 + 0.071d0 * inputs%T)
        k_a = 5.69 + 0.019 * (inputs%T - 273.15) * 418.684 * 1d-5
        D_v = 0.211d-4 / (inputs%p / const%air_std_press) &
            * (inputs%T / 273d0)**1.94d0
        Rd = const%univ_gas_const / const%air_molec_weight
        Rv = const%univ_gas_const / const%water_molec_weight
        rho_air = inputs%p / (Rd * inputs%T)
        P0 = env_state_saturated_vapor_pressure_water_2(inputs%T)
        P0_ice = env_state_saturated_vapor_pressure_ice(inputs%T)
        Ls = Lv + (Cpv - Cpw) * (inputs%T - const%water_freeze_temp)
        G = (Ls / (Rv * inputs%T) - 1d0) * Ls / (k_a * fti * inputs%T) + &
            Rv * inputs%T / (D_v * fvi * P0_ice)
        !print*, Ls, Lv, D_v, k_a
        if (condense_saved_do_ice_shape) then
            Cap = condense_iceGrowth_capacitance(R, phi)
            d_Cap_d_R = condense_dC_dR(phi)
        else
            Cap = R
            d_Cap_d_R = 1d0
        end if
        if (condense_saved_do_ice_density) then
            if (condense_saved_ice_supersat_density .ge. 0d0) then
                spec_den = inputs%den_dep
            else
                spec_den = inputs%den_sub
            end if
        else
            spec_den = inputs%den_ice
        end if
        !if (condense_saved_do_ice_ventilation .and. .False.) then
        if (condense_saved_do_ice_ventilation) then
            dG_dR =  -((Ls / (Rv * inputs%T) - 1d0) * Ls / (K_a * inputs%T) / &
                    (fti**2) * inputs%dft_dR + Rv * inputs%T / &
                    (D_v * P0_ice) / (fvi**2) * inputs%dfv_dR)
        else
            dG_dR = 0d0
        end if
        !print*, dG_dR
        !Rdot = (inputs%H * P0 / P0_ice - 1d0) * Cap / (R**2 * inputs%den_ice * G)

        Rdot = (inputs%H * P0 / P0_ice - 1d0) * Cap / (R**2 * spec_den * G)
        !!!!!! Temporary
        !Rdot = 0.148645043 * Cap / (R**2 * spec_den * G)


        MRii_dot = 4d0 * const%pi * R**2 * spec_den * Rdot / &
            inputs%V_comp / rho_air
        !Hdot_i = -inputs%H * (Lv * Ls / (Rv * inputs%T**2 * &
        !        const%air_spec_heat) + Rv * inputs%p / &
        !        (P0_ice * inputs%H * Rd)) * MRii_dot
        Hdot_i = -(Rv * inputs%p / (P0 *  Rd)) * MRii_dot

        !dRdot_dR = -(inputs%H*P0/P0_ice - 1d0) * 2d0 * Cap / (inputs%den_ice * G * R**3)
        dRdot_dR = -(inputs%H*P0/P0_ice - 1d0) / (spec_den * G) * (2d0 / &
                R**3 * Cap - 1 / R**2 * d_Cap_d_R) - &
                (inputs%H*P0/P0_ice - 1d0) * Cap / (R**2 * spec_den * G**2) * dG_dR
        !print*, dRdot_dR
        dRdot_dH = P0 / P0_ice * Cap / (R**2 * spec_den * G)
        !dHdoti_dR = 0d0
        dHdoti_dR = -(Rv * inputs%p) / (P0_ice * Rd) / rho_air * 4 * const%pi * &
            spec_den / inputs%V_comp * ( 2 * R * Rdot + R**2 * dRdot_dR)
        !dHdoti_dH = -4*const%pi*R**2 *inputs%den_ice/inputs%V_comp/rho_air *&
        !    (Lv*Ls/(Rv*inputs%T**2*const%air_spec_heat)*Rdot + &
        !     (inputs%H*Lv*Ls/(Rv*inputs%T**2*const%air_spec_heat) + &
        !      Rv*inputs%p/(P0_ice*Rd)) * dRdot_dH)
        dHdoti_dH = Hdot_i / Rdot * dRdot_dH

        outputs%Ddot = 2d0 * Rdot
        outputs%Hdot_i = Hdot_i
        outputs%dDdot_dD = dRdot_dR
        outputs%dDdot_dH = 2d0 * dRdot_dH
        outputs%dHdoti_dD = dHdoti_dR / 2d0
        outputs%dHdoti_dH = dHdoti_dH
        !outputs%Hdot_i = 0d0
        !outputs%Ddot = 1d0
        !outputs%Hdot_i = 1d0
        !outputs%dDdot_dD = 1d200
        !outputs%dDdot_dH = -2d200
        !outputs%dHdoti_dD = 3d200
        !outputs%dHdoti_dH = -4d200
        !print*, Hdot_i * P0 / Rv / inputs%T + 4 * const%pi * R**2 * &
        !    inputs%den_ice * Rdot / inputs%V_comp

        !print*, "Ls=",Ls,"Rv=",Rv,"k_a=",k_a,"D_v=", D_v,&
        !    "dR=",outputs%Ddot/2,"Cap/R=", Cap/R,"gtp=",G**(-1d0), &
        !    "si=",(inputs%H * P0 / P0_ice - 1d0),"RhoDep=",spec_den,&
        !    "R=", R

    end if

    

  end subroutine iceGrowth_rates


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_SUNDIALS
  !> Compute the condensation rates (Ddot and Hdot) at the current
  !> value of the state (D and H).
  subroutine condense_vf_f(n_eqn, time, state_p, state_dot_p) bind(c)

    !> Length of state vector.
    integer(kind=c_int), value, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=c_double), value, intent(in) :: time
    !> Pointer to state data.
    type(c_ptr), value, intent(in) :: state_p
    !> Pointer to state_dot data.
    type(c_ptr), value, intent(in) :: state_dot_p

    real(kind=c_double), pointer :: state(:)
    real(kind=c_double), pointer :: state_dot(:)
    real(kind=dp) :: Hdot
    integer :: i_part
    type(condense_rates_inputs_t) :: inputs
    type(condense_rates_outputs_t) :: outputs
    type(iceGrowth_rates_inputs_t) :: inputs_ice
    type(iceGrowth_rates_outputs_t) :: outputs_ice

    condense_count_vf = condense_count_vf + 1

    call c_f_pointer(state_p, state, (/ n_eqn /))
    call c_f_pointer(state_dot_p, state_dot, (/ n_eqn /))

    inputs%T = condense_saved_env_state_initial%temp &
         + time * condense_saved_Tdot
    inputs%p = condense_saved_env_state_initial%pressure &
         + time * condense_saved_pdot
    inputs%Tdot = condense_saved_Tdot
    inputs%pdot = condense_saved_pdot
    inputs%H = state(n_eqn)

    inputs_ice%T = condense_saved_env_state_initial%temp &
         + time * condense_saved_Tdot
    inputs_ice%p = condense_saved_env_state_initial%pressure &
         + time * condense_saved_pdot
    inputs_ice%Tdot = condense_saved_Tdot
    inputs_ice%pdot = condense_saved_pdot
    inputs_ice%H = state(n_eqn)
    !print*, "state(n_eqn) = ", state(1), state(2), state(n_eqn), state_p
    !print*, state_p
    !print*, ""

    Hdot = 0d0
    ! TangWenhan
    do i_part = 1,(n_eqn - 1)
       if (.not. condense_saved_frozen(i_part)) then ! Water droplet
           inputs%D = state(i_part)
           inputs%D_dry = condense_saved_D_dry(i_part)
           inputs%V_comp = (inputs%T &
                * condense_saved_env_state_initial%pressure) &
                / (condense_saved_env_state_initial%temp * inputs%p) &
                / condense_saved_num_conc(i_part)
           inputs%kappa = condense_saved_kappa(i_part)
           call condense_rates(inputs, outputs)
           state_dot(i_part) = outputs%Ddot
           Hdot = Hdot + outputs%Hdot_i
       else ! Ice
           !print*, "frozen3"
           inputs_ice%D = state(i_part)
           inputs_ice%D_dry = condense_saved_D_dry(i_part)
           inputs_ice%V_comp = (inputs%T &
                * condense_saved_env_state_initial%pressure) &
                / (condense_saved_env_state_initial%temp * inputs%p) &
                / condense_saved_num_conc(i_part)
           inputs_ice%den_ice = condense_saved_den_ice(i_part)
           inputs_ice%den_dep = condense_saved_ice_density_dep(i_part)
           inputs_ice%den_sub = condense_saved_ice_density_sub(i_part)
           inputs_ice%ice_shape_phi = condense_saved_ice_shape_phi(i_part)
           inputs_ice%fvi = condense_saved_fv(i_part)
           inputs_ice%fti = condense_saved_ft(i_part)
           inputs_ice%dfv_dR = condense_saved_dfv_dR(i_part)
           inputs_ice%dft_dR = condense_saved_dft_dR(i_part)
           call iceGrowth_rates(inputs_ice, outputs_ice)
           !if (i_part .eq. 1) then
           !    print*, inputs_ice%T, inputs_ice%p, inputs_ice%D, inputs_ice%H, outputs_ice%Ddot
           !    stop
           !    continue
           !end if
           !print*, outputs_ice%Ddot / 2, inputs_ice%D / 2
           state_dot(i_part) = outputs_ice%Ddot
           Hdot = Hdot + outputs_ice%Hdot_i
       end if
    end do
    ! TangWenhan
    if (condense_saved_frozen(n_eqn - 1)) then  
        !print*, "frozen4"
        Hdot = Hdot + outputs_ice%Hdot_env
    else
        Hdot = Hdot + outputs%Hdot_env
    end if

    state_dot(n_eqn) = Hdot
    !print*, "state_dot = ", state_dot(1), state_dot(n_eqn)

  end subroutine condense_vf_f
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_SUNDIALS
  !> Compute the Jacobian given by the derivatives of the condensation
  !> rates (Ddot and Hdot) with respect to the input variables (D and
  !> H).
  subroutine condense_jac(n_eqn, time, state_p, dDdot_dD, dDdot_dH, &
       dHdot_dD, dHdot_dH)

    !> Length of state vector.
    integer(kind=c_int), intent(in) :: n_eqn
    !> Current time (s).
    real(kind=c_double), intent(in) :: time
    !> Pointer to current state vector.
    type(c_ptr), intent(in) :: state_p
    !> Derivative of Ddot with respect to D.
    real(kind=dp), intent(out) :: dDdot_dD(n_eqn - 1)
    !> Derivative of Ddot with respect to H.
    real(kind=dp), intent(out) :: dDdot_dH(n_eqn - 1)
    !> Derivative of Hdot with respect to D.
    real(kind=dp), intent(out) :: dHdot_dD(n_eqn - 1)
    !> Derivative of Hdot with respect to H.
    real(kind=dp), intent(out) :: dHdot_dH

    real(kind=c_double), pointer :: state(:)
    integer :: i_part
    type(condense_rates_inputs_t) :: inputs
    type(condense_rates_outputs_t) :: outputs
    ! TangWenhan
    type(iceGrowth_rates_inputs_t) :: inputs_ice
    type(iceGrowth_rates_outputs_t) :: outputs_ice

    call c_f_pointer(state_p, state, (/ n_eqn /))

    inputs%T = condense_saved_env_state_initial%temp &
         + time * condense_saved_Tdot
    inputs%p = condense_saved_env_state_initial%pressure &
         + time * condense_saved_pdot
    inputs%Tdot = condense_saved_Tdot
    inputs%pdot = condense_saved_pdot
    inputs%H = state(n_eqn)

    inputs_ice%T = condense_saved_env_state_initial%temp &
         + time * condense_saved_Tdot
    inputs_ice%p = condense_saved_env_state_initial%pressure &
         + time * condense_saved_pdot
    inputs_ice%Tdot = condense_saved_Tdot
    inputs_ice%pdot = condense_saved_pdot
    inputs_ice%H = state(n_eqn)

    dHdot_dH = 0d0
    do i_part = 1,(n_eqn - 1)
       if (.not. condense_saved_frozen(i_part)) then
           inputs%D = state(i_part)
           inputs%D_dry = condense_saved_D_dry(i_part)
           inputs%V_comp = (inputs%T &
                * condense_saved_env_state_initial%pressure) &
                / (condense_saved_env_state_initial%temp * inputs%p) &
                / condense_saved_num_conc(i_part)
           inputs%kappa = condense_saved_kappa(i_part)
           call condense_rates(inputs, outputs)
           dDdot_dD(i_part) = outputs%dDdot_dD
           dDdot_dH(i_part) = outputs%dDdot_dH
           dHdot_dD(i_part) = outputs%dHdoti_dD + outputs%dHdotenv_dD
           dHdot_dH = dHdot_dH + outputs%dHdoti_dH
       else
           !print*, "frozen5"
           inputs_ice%D = state(i_part)
           inputs_ice%D_dry = condense_saved_D_dry(i_part)
           inputs_ice%V_comp = (inputs%T &
                * condense_saved_env_state_initial%pressure) &
                / (condense_saved_env_state_initial%temp * inputs%p) &
                / condense_saved_num_conc(i_part)
           inputs_ice%den_ice = condense_saved_den_ice(i_part)
           inputs_ice%den_dep = condense_saved_ice_density_dep(i_part)
           inputs_ice%den_sub = condense_saved_ice_density_sub(i_part)
           inputs_ice%ice_shape_phi = condense_saved_ice_shape_phi(i_part)
           inputs_ice%fvi = condense_saved_fv(i_part)
           inputs_ice%fti = condense_saved_ft(i_part)
           inputs_ice%dfv_dR = condense_saved_dfv_dR(i_part)
           inputs_ice%dft_dR = condense_saved_dft_dR(i_part)
           call iceGrowth_rates(inputs_ice, outputs_ice)
           dDdot_dD(i_part) = outputs_ice%dDdot_dD
           dDdot_dH(i_part) = outputs_ice%dDdot_dH
           dHdot_dD(i_part) = outputs_ice%dHdoti_dD + outputs_ice%dHdotenv_dD
           dHdot_dH = dHdot_dH + outputs_ice%dHdoti_dH
       end if
    end do
    !print*, "dHdot_dH (before) = ", dHdot_dH
    if (condense_saved_frozen(n_eqn - 1)) then
        !print*, "frozen6"
        !print*, "A", outputs_ice%dHdotenv_dH
        dHdot_dH = dHdot_dH + outputs_ice%dHdotenv_dH
    else
        !print*, "B", outputs%dHdotenv_dH
        dHdot_dH = dHdot_dH + outputs%dHdotenv_dH
    end if
    !print*, "dHdot_dH (after) = ", dHdot_dH

  end subroutine condense_jac
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_SUNDIALS
  !> Solve the system \f$ Pz = r \f$ where \f$ P = I - \gamma J \f$
  !> and \f$ J = \partial f / \partial y \f$. The solution is returned
  !> in the \f$ r \f$ vector.
  subroutine condense_jac_solve_f(n_eqn, time, state_p, state_dot_p, &
       rhs_p, gamma) bind(c)

    !> Length of state vector.
    integer(kind=c_int), value, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=c_double), value, intent(in) :: time
    !> Pointer to current state vector.
    type(c_ptr), value, intent(in) :: state_p
    !> Pointer to current state derivative vector.
    type(c_ptr), value, intent(in) :: state_dot_p
    !> Pointer to right-hand-side vector.
    type(c_ptr), value, intent(in) :: rhs_p
    !> Value of \c gamma scalar parameter.
    real(kind=c_double), value, intent(in) :: gamma

    real(kind=c_double), pointer :: state(:), state_dot(:), rhs(:)
    real(kind=c_double) :: soln(n_eqn)
    real(kind=dp) :: dDdot_dD(n_eqn - 1), dDdot_dH(n_eqn - 1)
    real(kind=dp) :: dHdot_dD(n_eqn - 1), dHdot_dH
    real(kind=dp) :: lhs_n, rhs_n
    real(kind=c_double) :: residual(n_eqn)
    real(kind=dp) :: rhs_norm, soln_norm, residual_norm
    integer :: i_part

    
    condense_count_solve = condense_count_solve + 1

    call condense_jac(n_eqn, time, state_p, dDdot_dD, dDdot_dH, &
         dHdot_dD, dHdot_dH)

    call c_f_pointer(state_p, state, (/ n_eqn /))
    call c_f_pointer(state_dot_p, state_dot, (/ n_eqn /))
    call c_f_pointer(rhs_p, rhs, (/ n_eqn /))
    !print*, "-------- condense_jac_solve_f --------"
    !print*, "state = ", state
    !print*, "state_dot = ", state_dot
    !print*, dDdot_dD(1), dDdot_dH(1), dHdot_dD(1), dHdot_dH
    !FIXME: write this all in matrix-vector notation, no i_part looping
    lhs_n = 1d0 - gamma * dHdot_dH
    rhs_n = rhs(n_eqn)
    do i_part = 1,(n_eqn - 1)
       lhs_n = lhs_n - (- gamma * dDdot_dH(i_part)) &
            * (- gamma * dHdot_dD(i_part)) / (1d0 - gamma * dDdot_dD(i_part))
       rhs_n = rhs_n - (- gamma * dHdot_dD(i_part)) * rhs(i_part) &
            / (1d0 - gamma * dDdot_dD(i_part))
    end do
    soln(n_eqn) = rhs_n / lhs_n

    do i_part = 1,(n_eqn - 1)
       soln(i_part) = (rhs(i_part) &
            - (- gamma * dDdot_dH(i_part)) * soln(n_eqn)) &
            / (1d0 - gamma * dDdot_dD(i_part))
    end do

    if (CONDENSE_DO_TEST_JAC_SOLVE) then
       ! (I - g J) soln = rhs

       ! residual = J soln
       residual(n_eqn) = sum(dHdot_dD * soln(1:(n_eqn-1))) &
            + dHdot_dH * soln(n_eqn)
       residual(1:(n_eqn-1)) = dDdot_dD * soln(1:(n_eqn-1)) &
            + dDdot_dH * soln(n_eqn)

       residual = rhs - (soln - gamma * residual)
       rhs_norm = sqrt(sum(rhs**2))
       soln_norm = sqrt(sum(soln**2))
       residual_norm = sqrt(sum(residual**2))
       write(0,*) 'rhs, soln, residual, residual/rhs = ', &
            rhs_norm, soln_norm, residual_norm, residual_norm / rhs_norm
    end if

    !print*, "--------------------------------------"
    rhs = soln
    !print*, "Fortran rhs(1) = ", rhs(1)

  end subroutine condense_jac_solve_f
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the water equilibrium state of a single particle.
  subroutine condense_equilib_particle(env_state, aero_data, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real(kind=dp) :: X, kappa, D_dry, D, g, dg_dD, a_w, daw_dD
    integer :: newton_step

    X = 4d0 * const%water_molec_weight * const%water_surf_eng &
         / (const%univ_gas_const * env_state%temp &
         * const%water_density)
    kappa = aero_particle_solute_kappa(aero_particle, aero_data)
    D_dry = aero_data_vol2diam(aero_data, &
         aero_particle_solute_volume(aero_particle, &
         aero_data))

    D = D_dry
    g = 0d0
    dg_dD = 1d0
    do newton_step = 1,20
       D = D - g / dg_dD
       a_w = (D**3 - D_dry**3) / (D**3 + (kappa - 1d0) * D_dry**3)
       daw_dD = 3d0 * D**2 * kappa * D_dry**3 &
            / (D**3 + (kappa - 1d0) * D_dry**3)**2
       g = env_state%rel_humid - a_w * exp(X / D)
       dg_dD = - daw_dD * exp(X / D) + a_w * exp(X / D) * (X / D**2)
    end do
    call warn_assert_msg(426620001, abs(g) < 1d3 * epsilon(1d0), &
         "convergence problem in equilibration")

    aero_particle%vol(aero_data%i_water) = aero_data_diam2vol(aero_Data, D) &
         - aero_data_diam2vol(aero_data, D_dry)

  end subroutine condense_equilib_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Call condense_equilib_particle() on each particle in the aerosol
  !> to ensure that every particle has its water content in
  !> equilibrium.
  subroutine condense_equilib_particles(env_state, aero_data, aero_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state

    integer :: i_part
    real(kind=dp) :: reweight_num_conc(aero_state_n_part(aero_state))

    ! We're modifying particle diameters, so bin sorting is now invalid
    aero_state%valid_sort = .false.

    call aero_state_num_conc_for_reweight(aero_state, aero_data, &
         reweight_num_conc)
    do i_part = aero_state_n_part(aero_state),1,-1
       call condense_equilib_particle(env_state, aero_data, &
            aero_state%apa%particle(i_part))
    end do
    ! adjust particles to account for weight changes
    call aero_state_reweight(aero_state, aero_data, reweight_num_conc)

  end subroutine condense_equilib_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(kind=dp) function get_Gamma_for_ice_growth(env_state)
    type(env_state_t), intent(in) :: env_state
    real(kind=dp) :: celsius, IGR1, IGR2, weight

    celsius = env_state%temp - const%water_freeze_temp
    if (celsius .gt. -1) then
        celsius = -1
    end if
    if (celsius .lt. -60) then
        celsius = -60
    end if
    weight = (ABS(real(int(celsius))) + 1.0) - ABS(celsius)
    IGR1 = condense_saved_gammaFindTrip(int(celsius)*(-1))
    IGR2 = condense_saved_gammaFindTrip((int(celsius)*(-1))+1)
    !print*, "IGR1", IGR1, "IGR2", IGR2
    get_Gamma_for_ice_growth = weight*IGR1 + (1.0-weight)*IGR2
    
  end function get_Gamma_for_ice_growth

  real(kind=dp) function condense_update_ice_shape_phi(phi0,&
          volume0, volume1, env_state, fva, fvc)
    real(kind=dp), intent(in) :: phi0  
    real(kind=dp), intent(in) :: volume0
    real(kind=dp), intent(in) :: volume1
    real(kind=dp), intent(in) :: fva, fvc
    type(env_state_t), intent(in) :: env_state
    real(kind=dp) :: IGR, II, IGRv
    !integer :: i_part
    IGR = get_Gamma_for_ice_growth(env_state)
    IGRv = IGR * fvc / fva
    II = (IGRv - 1) / (IGRv + 2)
    !print*, volume1, volume0, phi0 * (volume1 / volume0) ** II
    condense_update_ice_shape_phi = phi0 * (volume1 / volume0) ** II

  end function condense_update_ice_shape_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp) function condense_update_ice_density(volume0, volume1, &
          den_ice, den_dep, den_sub)
    real(kind=dp), intent(in) :: volume0, volume1, den_ice, den_dep, den_sub
    real(kind=dp) :: Vmin, beta, den_ice_new

    if (volume1 .ge. volume0) then
        den_ice_new = den_ice * (volume0 / volume1) + den_dep * &
            (1 - volume0 / volume1)
    else
        den_ice_new = den_ice * (volume0 / volume1) + den_sub * &
            (1 - volume0 / volume1)
        den_ice_new = min(den_ice_new, const%reference_ice_density)
        Vmin = 4d0 / 3d0 * const%pi * condense_Rice_avg**3
        beta = const%reference_ice_density / den_ice_new / log(Vmin / volume0)
        den_ice_new = den_ice_new + const%reference_ice_density * &
            (volume1 ** beta  - volume0 ** beta) / (Vmin ** beta)
    end if
    condense_update_ice_density = den_ice_new

  end function condense_update_ice_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp) function condense_iceGrowth_capacitance(R, phi)
    real(kind=dp), intent(in) :: R, phi
    real(kind=dp) :: V, a, c, Cap

    V = 4d0/3d0 * const%pi * R**3d0
    a = (3d0/(4d0*const%pi) * (V / phi)) ** (1d0/3d0)
    c = phi * a
    if (phi .eq. 1) then
        Cap = a
    else if(phi .gt. 1) then
        Cap = sqrt(c**2d0 - a**2d0) / log(  phi * (1d0 + sqrt(1d0 - phi**(-2)))  )
    else
        Cap = sqrt(a**2d0 - c**2d0) / asin(sqrt(1d0 - phi**2d0))
    end if
    condense_iceGrowth_capacitance = Cap

  end function condense_iceGrowth_capacitance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp) function condense_dC_dR(phi)
    real(kind=dp), intent(in) :: phi
    real(kind=dp) :: dC_dR

    if (phi .eq. 1) then
        dC_dR = 1d0
    else if(phi .gt. 1) then
        dC_dR = phi ** (2d0/3d0) * sqrt(1d0 - phi ** (-2d0)) / log( phi *&
                (1d0 + sqrt(1d0 - phi ** (-2d0))) )
    else
        dC_dR = phi ** (-1d0/3d0) * sqrt(1d0 - phi**2d0) / &
            asin(sqrt(1d0 - phi**2d0))
    end if
    condense_dC_dR = dC_dR

  end function condense_dC_dR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_ice_supersat_density(env_state)
    type(env_state_t), intent(in) :: env_state
    real(kind=dp) :: RH, es, ei, T, e, Rv
    T = env_state%temp
    RH = env_state%rel_humid
    es = env_state_saturated_vapor_pressure_water_2(T)
    ei = env_state_saturated_vapor_pressure_ice(T)
    e = es * RH
    Rv = const%univ_gas_const / const%water_molec_weight
    condense_saved_ice_supersat_density = (e - ei) / (Rv * T) * 1000d0
    !print*, "el=",e,"ei=",ei,"Rv=",Rv,"T=",T

  end subroutine condense_ice_supersat_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_ice_density_dep_sub(Vice, den_ice, env_state, dep, sub, &
          do_ice_shape, ice_dep_density_scheme_type, aero_particle)

    type(aero_particle_t), intent(inout) :: aero_particle
    type(env_state_t), intent(in) :: env_state
    real(kind=dp), intent(in) :: Vice, den_ice
    logical, intent(in) :: do_ice_shape
    integer, intent(in) :: ice_dep_density_scheme_type
    real(kind=dp), intent(out) :: dep, sub
    real(kind=dp) :: Vice_avg, T

    
    !T = env_state%temp

    !> Calculate deposition density.
    if (ice_dep_density_scheme_type .eq. &
         CONDENSE_ICE_DEP_DENSITY_SCHEME_CHENLAMB) then
       dep = condense_ice_depden_chenlamb(env_state, do_ice_shape)
    else if (ice_dep_density_scheme_type .eq. &
         CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_CTL) then
       call condense_ice_depden_pokrifka_ctl(aero_particle, env_state, dep)
    else if (ice_dep_density_scheme_type .eq. &
         CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_LOCAL) then
       dep = condense_ice_depden_pokrifka_local(env_state)
    end if

    !!! Temporary
    !dep = const%reference_ice_density * exp(- 3d0 * max( &
    !        0.206045434 - 5d-2, 0d0) / IGR)
    !print*, condense_saved_ice_supersat_density, IGR, dep
    !print*, "hhh", dep, - 3d0 * max( condense_saved_ice_supersat_density - &
    !            5d-2, 0d0) / IGR, IGR

    !> Calculate sublimation density.
    !Vice = 4d0 / 3d0 * const%pi * Rice**3
    Vice_avg = 4d0 / 3d0 * const%pi * condense_Rice_avg**3
    if (Vice .ge. Vice_avg) then
        sub = den_ice
    else
        sub = const%reference_ice_density
    end if

  end subroutine condense_ice_density_dep_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp) function condense_ice_depden_chenlamb(env_state, do_ice_shape)
     type(env_state_t), intent(in) :: env_state

     logical, intent(in) :: do_ice_shape
     real(kind=dp) :: IGR

    if (do_ice_shape) then
       IGR = get_Gamma_for_ice_growth(env_state)
    else
       IGR = 1d0
    end if
    condense_ice_depden_chenlamb = const%reference_ice_density * &
         exp(- 3d0 * max(condense_saved_ice_supersat_density - 5d-2, 0d0) / IGR)
  end function condense_ice_depden_chenlamb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_access_pokrifka_depden_dist()
    integer :: n_si, n_den, i_si, i_den
    open(61, file = 'pokrifka_rhodep.dat')
    read(61, *) n_si
    read(61, *) n_den
    !allocate(condense_saved_pokrifka_depden_dist_si(n_si))
    !allocate(condense_saved_pokrifka_depden_dist_den(n_den))
    !allocate(condense_saved_pokrifka_depden_cdf(n_si, n_den))
    read(61, *) (condense_saved_pokrifka_depden_dist_si(i_si),&
         i_si = 1, n_si)
    read(61, *) (condense_saved_pokrifka_depden_dist_den(i_den),&
         i_den = 1, n_den)
    read(61, *) ((condense_saved_pokrifka_depden_cdf(i_si, i_den),&
         i_si = 1, n_si), i_den = 1, n_den)
    close(61)
    !CONDENSE_POKRIFKA_SI_DIMENSION = n_si
    !CONDENSE_POKRIFKA_DEPDEN_DIMENSION = n_den
    condense_pokrifka_dendist_accessed = .true.
  end subroutine  condense_access_pokrifka_depden_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function condense_pokrifka_find_nearest_si_index(si)
    real(kind=dp) :: si, si_lower, si_upper = -1d0
    integer :: ind_si = -1, i_si

    si_lower = 0d0
    do i_si = 1, CONDENSE_POKRIFKA_SI_DIMENSION
       si_upper = condense_saved_pokrifka_depden_dist_si(i_si)
       if (si_upper .ge. si) then
          if ((i_si .eq. 0) .or. (si - si_lower .gt. si_upper - si)) then
             ind_si = i_si
          else
             ind_si = i_si - 1
          end if
          exit
       end if
       si_lower = si_upper
    end do
    if (si_upper .lt. si) then
       ind_si = CONDENSE_POKRIFKA_SI_DIMENSION
    end if
    condense_pokrifka_find_nearest_si_index = ind_si

    !print*, "Sampling Si = ", si, condense_saved_pokrifka_depden_dist_si(ind_si)

  end function condense_pokrifka_find_nearest_si_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp) function condense_pokrifka_depden_sampling(env_state)
    type(env_state_t), intent(in) :: env_state
    real(kind=dp) :: si, P0, P0_ice, P_vapor, rand_var
    real(kind=dp) :: depden_lower, depden_upper, alpha
    real(kind=dp) :: thres_lower, thres_upper = -1d0, depden
    integer :: ind_si, ind_den = -1, i_den

    P0 = env_state_saturated_vapor_pressure_water_2(env_state%temp)
    P0_ice = env_state_saturated_vapor_pressure_ice(env_state%temp)
    P_vapor = env_state%rel_humid * P0
    si = (P_vapor - P0_ice) / (P0 - P0_ice)
    ind_si = condense_pokrifka_find_nearest_si_index(si)
    rand_var = pmc_random()
    thres_lower = 0d0
    do i_den = 1, CONDENSE_POKRIFKA_DEPDEN_DIMENSION
       thres_upper = condense_saved_pokrifka_depden_cdf(ind_si, i_den) + 1e-12
       if (thres_upper .ge. rand_var) then
          ind_den = i_den
          exit
        end if
        thres_lower = thres_upper
    end do
    if (ind_den .eq. 1) then
       depden = condense_saved_pokrifka_depden_dist_den(ind_den)
    else
       alpha = (rand_var - thres_lower) / (thres_upper - thres_lower)
       depden_upper = condense_saved_pokrifka_depden_dist_den(ind_den)
       depden_lower = condense_saved_pokrifka_depden_dist_den(ind_den - 1)
       depden = depden_lower * (1 - alpha) + depden_upper * alpha
    end if
    condense_pokrifka_depden_sampling = depden

    !print*, "cdf:", thres_lower, rand_var, thres_upper

  end function condense_pokrifka_depden_sampling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !subroutine condense_ice_depden_pokrifka_ctl_initialization(aero_state)
  !  type(aero_state_t), intent(in) :: aero_state
  !  allocate(condense_saved_pokrifka_depden_part(aero_state_n_part(aero_state)))
  !  condense_saved_pokrifka_depden_part(:) = const%nan
  !  condense_pokrifka_depden_part_allocated = .true.
  !end subroutine condense_ice_depden_pokrifka_ctl_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_ice_depden_pokrifka_ctl(aero_particle, env_state, depden)
    type(aero_particle_t), intent(inout) :: aero_particle
    type(env_state_t), intent(in) :: env_state
    real(kind=dp), intent(out) :: depden
    if (aero_particle%depden .ge. 0d0) then
       depden = aero_particle%depden
       !print*, "Already sampled depden = ", depden 
    else
       !print*, "Sample for new;  old = ", &
             !aero_particle%depden
       depden = condense_pokrifka_depden_sampling(env_state)
       aero_particle%depden = depden
    end if
    !condense_ice_depden_pokrifka_ctl = depden
  end subroutine condense_ice_depden_pokrifka_ctl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp) function condense_ice_depden_pokrifka_local(env_state)
    type(env_state_t), intent(in) :: env_state
    condense_ice_depden_pokrifka_local = &
         condense_pokrifka_depden_sampling(env_state)
  end function condense_ice_depden_pokrifka_local

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_ice_depden_pokrifka_fix_initialization(env_state)
    type(env_state_t), intent(in) :: env_state
  end subroutine condense_ice_depden_pokrifka_fix_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp) function condense_ice_depden_pokrifka_fix(env_state)
    type(env_state_t), intent(in) :: env_state
    condense_ice_depden_pokrifka_fix = 0
  end function condense_ice_depden_pokrifka_fix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_ice_ventilation(Rice_array, n_part, env_state)
    implicit none
    type(env_state_t), intent(in) :: env_state
    integer, intent(in) :: n_part
    real(kind=dp), intent(in) :: Rice_array(n_part)
    real(kind=dp) :: T, press, den_ice, r_ice, phi, a_ice, c_ice, A, L, epsl_a
    real(kind=dp) :: xi, am, bm, Nre_i, Nsc, Npr, Vti, Xvent_i, Xtherm_i, m_ice
    real(kind=dp) :: Dv, Kt, kvisc, Rd, Rv, rho_air, bv1, bv2, gv, bt1, bt2, gt
    real(kind=dp) :: fvi, fti, fva_i, fvc_i, dfv_dr, dft_dr
    real(kind=dp) :: dm_dr, dxi_dr, dNre_dr, dXtherm_dr, dXvent_dr
    integer :: i_part
    T = env_state%temp
    press = env_state%pressure
    kt = 5.69 + 0.019 * (T - 273.15) * 418.684 * 1d-5 
    Dv = 0.211d-4 / (press / const%air_std_press) * (T / 273d0)**1.94d0 
    Rd = const%univ_gas_const / const%air_molec_weight
    Rv = const%univ_gas_const / const%water_molec_weight
    rho_air = press / (Rd * T) 
    epsl_a = (1.718 + 0.0049 * (T - 273.15) - 1.2d-5 * (T - 273.15)**2) * 1d-5
    kvisc = epsl_a / rho_air
    Nsc = kvisc / Dv
    Npr = kvisc / Kt
    do i_part = 1, n_part
        if (.not. condense_saved_frozen(i_part)) then
            cycle
        end if
        r_ice = Rice_array(i_part)
        den_ice = condense_saved_den_ice(i_part)
        phi = condense_saved_ice_shape_phi(i_part)
        a_ice = r_ice * phi**(-1d0/3d0)
        c_ice = r_ice * phi**(2d0/3d0)
        if (phi .le. 1d0) then
            A = const%pi * a_ice**2
            L = 2 * a_ice
        else
            A = const%pi * a_ice * c_ice
            L = 2 * c_ice
        end if
        m_ice = 4d0/3d0 * const%pi * r_ice**3 * den_ice
        xi = 2 * m_ice * (den_ice - rho_air) * const%std_grav * rho_air * L**2 &
                / (A * epsl_a**2 * den_ice)
        if (xi .le. 10) then
            am = 0.04394
            bm = 0.970
        else if (xi .le. 585) then
            am = 0.06049
            bm = 0.831
        else if (xi .le. 1.56d5) then
            am = 0.2072
            bm = 0.638
        else if (xi .le. 1e8) then
            am = 1.0865
            bm = 0.499
        else
            !print*, "Ah Ah Ah"
            am = 1.0865
            bm = 0.499
        end if
        Nre_i = am * xi**bm
        Vti = Nre_i * epsl_a / (rho_air * L)
        Xvent_i = Nsc**(1d0/3d0) * Nre_i**(1d0/2d0)
        Xtherm_i = Npr**(1d0/3d0) * Nre_i**(1d0/2d0)
    
        if (Xvent_i .le. 1) then
            bv1 = 1d0
            bv2 = 0.14
            gv = 2d0
        else
            bv1 = 0.86
            bv2 = 0.28
            gv = 1d0
        end if

        if (Xtherm_i .lt. 1.4) then
            bt1 = 1d0
            bt2 = 0.108
            gt = 2d0
        else
            bt1 = 0.78
            bt2 = 0.308
            gt = 1d0
        end if

        fvi = bv1 + bv2 * Xvent_i ** gv
        fti = bt1 + bt2 * Xtherm_i ** gt
        fva_i = bv1 + bv2 * (Xvent_i ** gv) * (a_ice / r_ice) ** (gv / 2d0)
        fvc_i = bv1 + bv2 * (Xvent_i ** gv) * (c_ice / r_ice) ** (gv / 2d0)

        dm_dr = 4 * const%pi * r_ice**2 * den_ice
        if (phi .le. 1d0) then
            dxi_dr = 2 * (den_ice - rho_air) * const%std_grav * rho_air /&
                    (epsl_a**2 * den_ice) * (4d0 / const%std_grav) * dm_dr
        else
            dxi_dr = 2 * (den_ice - rho_air) * const%std_grav * rho_air/& 
                    (epsl_a**2 * den_ice) * (4d0 / const%std_grav * phi) * dm_dr
        end if
        dNre_dr = am * bm * xi**(bm-1) * dxi_dr
        dXvent_dr = 1d0 / 2d0 * Nsc**(1d0/3d0) * Nre_i **(-1d0/2d0) * dNre_dr
        dXtherm_dr = 1d0 / 2d0 * Npr**(1d0/3d0) * Nre_i **(-1d0/2d0) * dNre_dr
        dfv_dr = bv2 * gv * Xvent_i**(gv-1) * dXvent_dr
        dft_dr = bt2 * gt * Xtherm_i**(gt-1) * dXtherm_dr


        condense_saved_fv(i_part) = fvi
        condense_saved_ft(i_part) = fti
        condense_saved_fva(i_part) = fva_i
        condense_saved_fvc(i_part) = fvc_i
        condense_saved_dfv_dr(i_part) = dfv_dr
        condense_saved_dft_dr(i_part) = dft_dr

        !print*,"fv-1=",fvi-1,"ft-1=",fti-1,"fvc-1=",fvc_i-1,"fva-1=",fva_i-1,&
        !    "Xvent=",Xvent_i,"Xtherm=",Xtherm_i, "X=",xi,"epsl=",epsl_a, &
        !    "area=", &
        !    A,"nre=",Nre_i,"npr=",Npr,"am=",am,"bm=",bm,"bv1=",bv1,"bv2=",bv2,"gv=",gv,&
        !    "bt1=",bt1,"bt2=",bt2,"gt=",gt, "nsc=", Nsc

    end do

  end subroutine condense_ice_ventilation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read the specification for a kernel type from a spec file and
  !> generate it.
  subroutine spec_file_read_ice_dep_density_scheme_type(file, ice_dep_density_scheme_type)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Kernel type.
    integer, intent(out) :: ice_dep_density_scheme_type

    character(len=SPEC_LINE_MAX_VAR_LEN) :: dep_scheme

    !> \page input_format_ice_dep_density_scheme Input File Format: ice
    !!   deposition density scheme.
    !!
    !! The ice deposition density  scheme is specified by the parameter:
    !!   - \b ice_dep_density_scheme (string): the type of ice deposition density scheme
    !!   must be one of: \c sedi for the gravitational sedimentation
    !!   kernel; \c additive for the additive kernel; \c constant
    !!   for the constant kernel; \c brown for the Brownian kernel,
    !!   or \c zero for no immersion freezing
    !!
    !! If \c ice_dep_density_scheme is \c additive, the kernel coefficient needs to be
    !! provided using the \c additive_kernel_coeff parameter
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'ice_dep_density_scheme', dep_scheme)
    if (trim(dep_scheme) == 'chenlamb') then
       ice_dep_density_scheme_type = CONDENSE_ICE_DEP_DENSITY_SCHEME_CHENLAMB
    elseif (trim(dep_scheme) == 'pokrifka_ctl') then
       ice_dep_density_scheme_type = &
            CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_CTL
    elseif (trim(dep_scheme) == 'pokrifka_local') then
       ice_dep_density_scheme_type = &
            CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_LOCAL
    elseif (trim(dep_scheme) == 'pokrifka_fix') then
       ice_dep_density_scheme_type = &
            CONDENSE_ICE_DEP_DENSITY_SCHEME_POKRIFKA_FIX
    else
       call spec_file_die_msg(321125901, file, &
        "Unknown ice deposition density scheme: " // trim(dep_scheme))
    end if

  end subroutine spec_file_read_ice_dep_density_scheme_type



end module pmc_condense

