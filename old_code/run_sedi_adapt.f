C Simulation with sedimentation kernel and adaptive timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo

      ! MM = 2e5, n_bin = 160, t_max = 12000: 123 minutes
      ! MM = 1e6, n_bin = 160, t_max = 7200:   43 minutes
      ! MM = 1e7, n_bin = 160, t_max = 1800:   10 minutes
      ! MM = 1e8, n_bin = 220, t_max = 600:     3 minutes
      integer MM, n_bin, n_loop, scal
      real*8 t_max, rho_p, N_0, t_print, t_progress
      real*8 r_samp_max, del_t_max, V_0
      parameter (MM = 10000000)        ! number of particles
      parameter (n_bin = 160)          ! number of bins
      parameter (n_loop = 1)           ! number of loops
      parameter (scal = 3)             ! scale factor for bins
      parameter (t_max = 360d0)        ! total simulation time (seconds)
      parameter (rho_p = 1000d0)       ! particle density (kg/m^3)
      parameter (N_0 = 1d9)            ! particle number concentration (#/m^3)
      parameter (t_print = 60d0)       ! interval between printing (s)
      parameter (t_progress = 1d0)     ! interval between progress (s)
      parameter (r_samp_max = 0.005d0) ! maximum sampling ratio per timestep
      parameter (del_t_max = 1d0)      ! maximum timestep (s)
      parameter (V_0 = 4.1886d-15)     ! mean volume of initial distribution (m^3)

      integer M, i_loop
      real*8 V(MM), V_comp, dlnr
      real*8 bin_v(n_bin), bin_r(n_bin)
      real*8 bin_g(n_bin)
      integer n_ini(n_bin), bin_n(n_bin)

      external kernel_sedi

      open(30,file='out_sedi_adapt.d')
      call print_header(n_loop, n_bin, nint(t_max / t_print) + 1)
      call srand(12)

      do i_loop = 1,n_loop

         call make_grid(n_bin, scal, rho_p, bin_v, bin_r, dlnr)
         call init_exp(MM, V_0, dlnr, n_bin, bin_v, bin_r, n_ini)
         !call init_bidisperse(MM, n_bin, n_ini)
         call compute_volumes(n_bin, MM, n_ini, bin_v, dlnr, V, M)
         V_comp = M / N_0

         call mc_adapt(MM, M, V, V_comp,
     &        n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &        kernel_sedi, t_max, t_print, t_progress,
     &        r_samp_max, del_t_max, i_loop)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC