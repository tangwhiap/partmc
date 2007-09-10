! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Datatype for process_state.

module pmc_process_spec

  use pmc_inout
  use pmc_aero_data

  integer, parameter :: PROCESS_SPEC_TYPE_LEN = 40
  integer, parameter :: PROCESS_SPEC_FILE_LEN = 200

  type process_spec_t
     character(len=PROCESS_SPEC_TYPE_LEN) :: type ! processing type
     character(len=PROCESS_SPEC_FILE_LEN) :: suffix ! output file suffix
     integer :: n_step                  ! number of steps for histograms
     real*8 :: min_real                 ! minimum histogram value (real)
     real*8 :: max_real                 ! maximum histogram value (real)
     integer :: min_int                 ! minimum histogram value (int)
     integer :: max_int                 ! maximum histogram value (int)
     logical :: log_scale               ! use a log-scale?
     character(len=AERO_NAME_LEN), pointer :: a_species(:) ! comp A species
     character(len=AERO_NAME_LEN), pointer :: b_species(:) ! comp B species
  end type process_spec_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_spec_alloc(process_spec)

    ! Allocate and initialize storage.

    type(process_spec_t), intent(out) :: process_spec ! structure to alloc

    process_spec%type = "none"
    process_spec%suffix = "none"
    process_spec%n_step = 0
    process_spec%min_real = 0d0
    process_spec%max_real = 0d0
    process_spec%min_int = 0
    process_spec%max_int = 0
    process_spec%log_scale = .false.
    allocate(process_spec%a_species(0))
    allocate(process_spec%b_species(0))

  end subroutine process_spec_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_spec_free(process_spec)

    ! Free storage.

    type(process_spec_t), intent(inout) :: process_spec ! structure to free

    deallocate(process_spec%a_species)
    deallocate(process_spec%b_species)

  end subroutine process_spec_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_spec_copy(process_spec_from, process_spec_to)

    ! Copy all data.

    type(process_spec_t), intent(in) :: process_spec_from ! source of data
    type(process_spec_t), intent(inout) :: process_spec_to ! dest of data

    process_spec_to%type = process_spec_from%type
    process_spec_to%suffix = process_spec_from%suffix
    process_spec_to%n_step = process_spec_from%n_step
    process_spec_to%min_real = process_spec_from%min_real
    process_spec_to%max_real = process_spec_from%max_real
    process_spec_to%min_int = process_spec_from%min_int
    process_spec_to%max_int = process_spec_from%max_int
    process_spec_to%log_scale = process_spec_from%log_scale
    allocate(process_spec_to%a_species(size(process_spec_from%a_species)))
    allocate(process_spec_to%b_species(size(process_spec_from%b_species)))
    process_spec_to%a_species = process_spec_from%a_species
    process_spec_to%b_species = process_spec_from%b_species

  end subroutine process_spec_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec_list(filename, process_spec_list)

    ! Read a list of process_specs from the named file.

    character(len=*), intent(in) :: filename ! filename to read
    type(process_spec_t), pointer :: process_spec_list(:) ! data to read

    type(inout_file_t) :: file
    type(process_spec_t), pointer :: new_process_spec_list(:)
    type(process_spec_t) :: process_spec
    integer :: i, n_process_spec
    logical :: eof
    
    call inout_open_read(filename, file)
    n_process_spec = 0
    allocate(process_spec_list(n_process_spec))
    call inout_read_process_spec(file, process_spec, eof)
    do while (.not. eof)
       n_process_spec = n_process_spec + 1
       allocate(new_process_spec_list(n_process_spec))
       call process_spec_copy(process_spec, &
            new_process_spec_list(n_process_spec))
       call process_spec_free(process_spec)
       do i = 1,(n_process_spec - 1)
          call process_spec_copy(process_spec_list(i), new_process_spec_list(i))
          call process_spec_free(process_spec_list(i))
       end do
       deallocate(process_spec_list)
       process_spec_list => new_process_spec_list
       nullify(new_process_spec_list)
       call inout_read_process_spec(file, process_spec, eof)
    end do
    
  end subroutine inout_read_process_spec_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec(file, process_spec, eof)

    ! Read from an inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(process_spec_t), intent(out) :: process_spec ! data to read
    logical :: eof                      ! if eof instead of reading data

    type(inout_line_t) :: line

    call inout_read_line(file, line, eof)
    if (.not. eof) then
       call inout_check_line_name(file, line, "process")
       call inout_check_line_length(file, line, 1)
       call process_spec_alloc(process_spec)
       process_spec%type = line%data(1)(1:PROCESS_SPEC_TYPE_LEN)
       call inout_read_string(file, "suffix", process_spec%suffix)
       if (process_spec%type == "kappa") then
          call inout_read_process_spec_kappa(file, process_spec)
       elseif (process_spec%type == "comp") then
          call inout_read_process_spec_comp(file, process_spec)
       elseif (process_spec%type == "n_orig_part") then
          call inout_read_process_spec_n_orig_part(file, process_spec)
       elseif (process_spec%type == "optic_absorb") then
          call inout_read_process_spec_optic(file, process_spec)
       elseif (process_spec%type == "optic_scatter") then
          call inout_read_process_spec_optic(file, process_spec)
       elseif (process_spec%type == "optic_extinct") then
          call inout_read_process_spec_optic(file, process_spec)
       elseif ((process_spec%type /= "env") &
            .and. (process_spec%type /= "gas") &
            .and. (process_spec%type /= "aero")) then
          write(0,*) 'ERROR: unknown process type on line ', &
               file%line_num, ' of file ', trim(file%name)
          call exit(1)
       end if
       call inout_line_free(line)
    end if

  end subroutine inout_read_process_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec_kappa(file, process_spec)

    ! Read kappa spec from an inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(process_spec_t), intent(out) :: process_spec ! data to read

    call inout_read_integer(file, "n_step", process_spec%n_step)
    call inout_read_real(file, "min", process_spec%min_real)
    call inout_read_real(file, "max", process_spec%max_real)

  end subroutine inout_read_process_spec_kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec_comp(file, process_spec)

    ! Read comp spec from an inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(process_spec_t), intent(out) :: process_spec ! data to read

    type(inout_line_t) :: line
    integer :: i

    call inout_read_integer(file, "n_step", process_spec%n_step)
    call inout_read_real(file, "min", process_spec%min_real)
    call inout_read_real(file, "max", process_spec%max_real)

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, "a_species")
    allocate(process_spec%a_species(size(line%data)))
    do i = 1,size(line%data)
       process_spec%a_species(i) = line%data(i)(1:AERO_NAME_LEN)
    end do
    call inout_line_free(line)

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, "b_species")
    allocate(process_spec%b_species(size(line%data)))
    do i = 1,size(line%data)
       process_spec%b_species(i) = line%data(i)(1:AERO_NAME_LEN)
    end do
    call inout_line_free(line)

  end subroutine inout_read_process_spec_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec_n_orig_part(file, process_spec)

    ! Read n_orig_part spec from an inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(process_spec_t), intent(out) :: process_spec ! data to read

    call inout_read_integer(file, "min", process_spec%min_int)
    call inout_read_integer(file, "max", process_spec%max_int)

  end subroutine inout_read_process_spec_n_orig_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec_optic(file, process_spec)

    ! Read optic spec from an inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(process_spec_t), intent(out) :: process_spec ! data to read

    call inout_read_integer(file, "n_step", process_spec%n_step)
    call inout_read_real(file, "min", process_spec%min_real)
    call inout_read_real(file, "max", process_spec%max_real)

  end subroutine inout_read_process_spec_optic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_process_spec