!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of SPECFEM3D_GLOBE version 7.0.0 and ASKI version 1.2.
!
!   SPECFEM3D_GLOBE version 7.0.0 and ASKI version 1.2 are free software: 
!   you can redistribute it and/or modify it under the terms of the GNU 
!   General Public License as published by the Free Software Foundation, 
!   either version 2 of the License, or (at your option) any later version.
!
!   SPECFEM3D_GLOBE version 7.0.0 and ASKI version 1.2 are distributed in 
!   the hope that they will be useful, but WITHOUT ANY WARRANTY; without 
!   even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
!   PURPOSE.  See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with SPECFEM3D_GLOBE version 7.0.0 and ASKI version 1.2.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------

!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
!
! United States and French Government Sponsorship Acknowledged.


  module specfem_for_ASKI_par
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! THE CONTENT OF THIS MODULE ORIGINATES FROM FILES
! specfem3D_par_ASKI.f90, specfem3D_for_ASKI.f90 
! FROM 
!    SPECFEM3D_Cartesian version 2.1 for ASKI 0.3
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

! specific module for ASKI extension package

  implicit none

! parameters for generating ASKI output
  integer :: ASKI_np_local
  integer, dimension(:), allocatable :: ASKI_np_local_all
  integer, dimension(:,:), allocatable :: ASKI_indx_local
  complex(kind=kind(1.d0)), dimension(:,:), allocatable :: ASKI_efactors_tapered !variable is also used in case there is no tapering!
  complex(kind=kind(1.d0)), dimension(:,:,:), allocatable :: ASKI_spectra_local_double
  complex, dimension(:,:,:), allocatable :: ASKI_spectra_local_single

  ! other temporary stuff
  character(len=500), dimension(:), allocatable :: ASKI_val_parfile
  character(len=100), dimension(:), allocatable :: ASKI_key_parfile
  ! chunk rotation stuff
  double precision, dimension(3,3) :: ASKI_Mrot_chunk1,ASKI_Mrot_chunk2,ASKI_Mrot_chunk3
  double precision :: ASKI_xmin_chunk1,ASKI_xmax_chunk1,ASKI_ymin_chunk1,ASKI_ymax_chunk1
  ! Green function source stuff
  logical :: ASKI_gf_is_NEZ,ASKI_gf_is_XYZ
  integer :: ASKI_gf_icomp


! content of Par_file_ASKI is following now:

!------------------------------------------------------------------------------------
! ASKI OUTPUT
!
! producing kernel green tensor / kernel displacement files
!------------------------------------------------------------------------------------

! parameter COMPUTE_ASKI_OUTPUT controls if any ASKI output files (i.e. kernel green tensor / 
! kernel displacement files) are produced
  logical :: COMPUTE_ASKI_OUTPUT

! Decide, whether to ONLY produce the main ASKI output file (at the beginning of this simulation, filename 'ASKI_outfile.main') 
! and then abort the simulation.
! This can be useful if you want to first run initBasics to check the ASKI output volume and the background model etc.
 logical :: ASKI_MAIN_FILE_ONLY

! choose to overwrite existing ASKI output files, or to abort if files already exist
  logical :: OVERWRITE_ASKI_OUTPUT

! absolute output filename including absolute path (will be used to open the output file as is),
! i.e. in case of a green tensor simulation, outfile should contain the extension "_[green_tensor_component]", like "_X","_Y","_Z"
  character(len=500) :: ASKI_outfile

! id of ASKI output (e.g. eventID, station name plus component, etc.)
  integer, parameter :: length_ASKI_output_ID = 13
  character(len=length_ASKI_output_ID) :: ASKI_output_ID

! flag whether deconvolution of a source time function should be applied
  logical :: ASKI_DECONVOLVE_STF

! flag whether the source is a Grenn function single force source or not
  logical :: COMPUTE_ASKI_GREEN_FUNCTION

! Green function source direction (i.e. respective receiver component of back propagation)
  integer, parameter :: length_ASKI_component = 2 ! always make consistent with the supported types of components
  character(len=length_ASKI_component) :: ASKI_GREEN_FUNCTION_COMPONENT

!------------------------------------------------------------------------------------
! FREQUENCY DISCRETIZATION
!
!   the terms df,jf have the following meaning:
!   the spectra are saved for all frequencies f = (jf)*df
!------------------------------------------------------------------------------------

! predefined df, that is used to evaluate spectrum (in case we want to do an inverse FT, we need to choose with care 
! df=1/length_of_time_series and suitably high frequency indices (dependent on frequency content), as
! we could lose periodicity (if in exp^(-i2pi(k)(n)/N) "N" is no integer, these are no roots of 1 anymore))
! save the spectra for frequencies f = (ASKI_jf)*ASKI_df (ASKI_nf many)
  double precision :: ASKI_df
  integer :: ASKI_nf
  integer, dimension(:), allocatable :: ASKI_jf

! choose precision of Discrete Fourier Transform, if there is enough memory available, it is highly recommended
! to use ASKI_DFT_double = .true. in which casedouble precision complex coefficients exp^(-i*2pi*f*t) are used
! and double complex spectra are hold in memory(single precision is written to file, though, but less roundoffs
! during transformation
! otherwise choose ASKI_DFT_double = .false. inwhich case single precision will be used
  logical :: ASKI_DFT_double

! decide whether the (oversampled, noisy, ...) time series should be tapered bya hanning taper(on tail)
! before (i.e. while) applying the discrete fourier transform
! if ASKI_DFT_apply_taper = .true. the value ofASKI_DFT_taper_percentage (between 0.0 and 1.0)definesthe amount of
! total time for which the hanning taper will be applied at thetail ofthe time series
  logical :: ASKI_DFT_apply_taper
  double precision :: ASKI_DFT_taper_percentage

!------------------------------------------------------------------------------------
! INVERSION GRID
!
! ASKI supports several types of inversion grids for FORWARD_METHOD = SPECFEM3D:
!
!   ASKI_type_inversion_grid = 1, (TYPE_INVERSION_GRID = schunkInversionGrid)
!      ASKI internal, but SPECFEM independent spherical inverison grid (simple chunk inversion grid)
!
!   ASKI_type_inversion_grid = 2, (TYPE_INVERSION_GRID = scartInversionGrid) NOT TO BE USED WITH SPECFEM3D GLOBE!
!      ASKI internal, but SPECFEM independent cartesian inversion grid:
!      the values for ASKI output are stored at all inner GLL points of spectral elements which lie
!      inside the ASKI output volume
!      ASKI loactes the coordinates of those points inside the inversion grid cells and computes
!      integration weights for them
!
!   ASKI_type_inversion_grid = 3, (TYPE_INVERSION_GRID = ecartInversionGrid)
!      external inversion grid provided e.g. by CUBIT, which may contain tetrahedra, as well as hexahedra
!      which must be defined in the coordinates used for global wavefield points (i.e. usually global Cartesian coordinates in km).
!      as in case of ASKI_type_inversion_grid = 2, ASKI output is stored at all inner GLL points of elements
!      which are inside the volume defined by ASKI_nchunk, ASKI_(cw)(latlon), ASKI_r(minmax).
!      ASKI locates the wavefield points inside the inversion grid and computes integration weights
!
!   ASKI_type_inversion_grid = 4, (TYPE_INVERSION_GRID = specfem3dInversionGrid)
!      use SPECFEM elements as inversion grid:
!      wavefield points are ALL GLL points of an element for elements which are (at least partly) inside the 
!      volume defined by ASKI_nchunk, ASKI_(cw)(latlon), ASKI_r(minmax)
!      additionally store the jacobians for all wavefield points
!      assume ncell = ntot_wp/(NGLLX*NGLLY*NGLLY) as the number of inversion grid cells, and the order of 
!      wavefield points accordingly (do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX; ip=ip+1 ....)
!
!   ASKI_type_inversion_grid = 5, (TYPE_INVERSION_GRID = chunksInversionGrid)
!      ASKI internal, but SPECFEM independent more elaborate spherical inverison grid supporting several chunks
!
! The following parameters may be used to define a spherical volume within which wavefield points are searched for:
!
! ASKI_nchunk defines the number of chunks of a chunk cubed sphere. 
! SPECFEM3D_GLOBE for ASKI supports 1,2,3 and 6 chunks (like the forward code itself)
! However, not all those values are supported by all types of inversion grids: 
!   the schunkInversionGrid only supports ASKI_nchunk = 1
!   the ecartInversionGrid,specfem3dInversionGrid and chunksInversionGrid support all values 1,2,3,6
!
! ASKI_clat,ASKI_clon define the center of the first chunks (not used for full sphere with 6 chunks)
! ASKI_wlat,ASKI_wlon define the width of the chunk, if ASKI_nchunk = 1. For nchunk = 2,3,6 always wlon=wlat=90.0 is used
! ASKI_rot_gamma defines the azimuthal rotation angle in degrees by which the 1,2 or 3 chunks are rotated (anti-clockwise) about the local vertical axis through the center of the first chunk. Not used for ASKI_nchunk = 6
!------------------------------------------------------------------------------------

! type of inversion grid, in order to produce correct output
  integer :: ASKI_type_inversion_grid

! number of chunks of a chunk cubed sphere (supported values: 1,2,3,6)
  integer :: ASKI_nchunk

! lateral width of ASKI output in degrees, used only for 1 chunk (for 2,3,6 chunks always 90.0 used)
  real :: ASKI_wlat
  real :: ASKI_wlon
! azimuthal rotation angle in degrees about the vertical through center of first chunk (only used for 1,2,3 chunks)
  real :: ASKI_rot_gamma
! center of first ASKI output chunk in geographical coordinates (deg) (only used for 1,2,3 chunks)
  real :: ASKI_clat
  real :: ASKI_clon
! global depht distribution of ASKI output volume, defined by maximum radius and minimum radius in km (used for all numbers of chunks)
  real :: ASKI_rmax
  real :: ASKI_rmin


contains

subroutine prepare_timerun_ASKI()

  use specfem_par,only: CUSTOM_REAL,SIZE_REAL,NPROCTOT,myrank

  implicit none

  integer :: iproc
  integer, parameter :: itag = 100
  integer, dimension(1) :: i_array_one_value

  if (CUSTOM_REAL /= SIZE_REAL) call exit_MPI_without_rank('so far, only single precision supported for ASKI output')

  call read_Par_file_ASKI()

  if (.not.COMPUTE_ASKI_OUTPUT) return

  call define_ASKI_chunk_rotations()

  ! depending on type of inversion grid, collect information on GLL points
  ! where ASKI output should be produced (like number of points, their indices, 
  ! coordinates, model values)

  ASKI_np_local = 0
  select case(ASKI_type_inversion_grid)
  case(1) ! ASKI internal, but SPECFEM independent inversion grids suitable for 1 chunk only (schunkInversionGrid), store at all inner GLL points
     if(ASKI_nchunk /= 1) call exit_MPI_without_rank("for ASKI_type_inversion_grid==1 (schunkInversionGrid), only "//&
          "1 chunk supported (i.e. ASKI_nchunk must be 1)")
     call search_ASKI_wavefield_points_inner_GLL_in_chunks()
  case(3,5) ! ASKI internal, but SPECFEM independent inversion grids, (ecartInversionGrid,chunksInversionGrid) so store at all inner GLL points
     call search_ASKI_wavefield_points_inner_GLL_in_chunks()
  case(4) ! use SPECFEM elements as inversion grid (specfem3dInversionGrid)
     call search_ASKI_wavefield_points_type_invgrid_4()
  case default
     call exit_MPI_without_rank('values for ASKI_type_inversion_grid other than 1,3,4,5 not supported yet')
  end select ! ASKI_type_inversion_grid
  call synchronize_all()

  ! in the routines search_ASKI_wavefield_points*, the following variables are defined:
  !   ASKI_np_local (number of ASKI wavefield points for this proc)
  !   ASKI_indx_local

  ! gather ASKI_np_local from everybody on proc 0
  if(myrank == 0) then
     allocate(ASKI_np_local_all(NPROCTOT))
     ASKI_np_local_all(1) = ASKI_np_local ! this is me, rank 0
     do iproc = 1,NPROCTOT-1
        ! receive ASKI_np_local from rank iproc
        call recv_i(i_array_one_value,1,iproc,itag)
        ASKI_np_local_all(iproc+1) = i_array_one_value(1)
     end do ! iproc

     if(sum(ASKI_np_local_all) .le. 0) then
        call write_ASKI_log('ERROR_ASKI.txt',"no ASKI wavefield points found at all, so no ASKI output can be computed")
        call exit_MPI_without_rank('no ASKI wavefield points found at all, so no ASKI output can be computed')
     end if

  else ! (myrank == 0)
     ! send ASKI_np_local to rank 0
     i_array_one_value(1) = ASKI_np_local
     call send_i(i_array_one_value,1,0,itag)

  end if ! (myrank == 0)

  ! define discrete fourier transform factors exp(...) once here, before time loop, plus checks, plus allcoation
  call prepare_ASKI_output()

  ! if the ASKI output is the wavefield of a Green function, define the respective single force source
  ! by modifying the array sourcearrays
  if(COMPUTE_ASKI_GREEN_FUNCTION) call define_ASKI_gf_source()

  if(myrank == 0) call write_ASKI_log_start()

  ! write wavefield points and kernel reference model (and jacobian, in case type_invgrid = 4)
  ! and frequency info to file
  call write_ASKI_main_file()

  if(ASKI_MAIN_FILE_ONLY) then
     ! wait until the main parfile has been written
     call synchronize_all()
     ! abort this run
     call exit_MPI_without_rank("logical parameter 'ASKI_MAIN_FILE_ONLY' in DATA/Par_file_ASKI requests to only "//&
          "write the main ASKI output file, hence aborting this run")
  end if


end subroutine prepare_timerun_ASKI
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine read_Par_file_ASKI()

  use specfem_par,only: myrank
  use constants,only: IMAIN

  implicit none

  character(len=601) :: line
  character(len=500) :: val
  integer :: npar,ios,IOASKI,eqindx

  ! open Par_file_ASKI and find number of valid lines
  call get_file_unit_ASKI(IOASKI)
  open(unit=IOASKI,file='DATA/'//'Par_file_ASKI',&
       form='formatted',status='old',action='read',iostat=ios)
  if(ios/=0) then
     close(IOASKI)
     COMPUTE_ASKI_OUTPUT = .false.
     if(myrank==0) call write_ASKI_log('LOG_ASKI_start.txt',"could not open file '"//'DATA/'//&
          "Par_file_ASKI', so no ASKI output is produced")
     return
  end if
  ! number of valid lines
  npar = 0
  do while(ios==0)
     read(IOASKI,"(a601)",iostat=ios) line
     if( len_trim(line) > 0 ) then
        line = adjustl(line)
        if(line(1:1) /= '#') then ! ignore comment lines
           eqindx = index(line,'=') ! only allow lines with at least one character in front of '='
           if(eqindx>1) npar = npar + 1
        end if
     end if
  end do
  close(IOASKI)

  if(npar == 0) call exit_MPI_without_rank("no valid lines in file '"//'DATA/'//"Par_file_ASKI'")
  allocate(ASKI_key_parfile(npar),ASKI_val_parfile(npar))

  ! now open again and store key,val pairs of valid lines
  call get_file_unit_ASKI(IOASKI)
  open(unit=IOASKI,file='DATA/'//'Par_file_ASKI',&
       form='formatted',status='old',action='read',iostat=ios)
  npar = 0
  do while(ios==0)
     read(IOASKI,"(a601)",iostat=ios) line
     if( len_trim(line) > 0 ) then
        line = adjustl(line)
        if(line(1:1) /= '#') then ! ignore comment lines
           eqindx = index(line,'=') ! only allow lines with at least one character in front of '='
           if(eqindx>1) then
                npar = npar + 1
                ASKI_key_parfile(npar) = line(1:eqindx-1)
                ASKI_val_parfile(npar) = adjustl(line(eqindx+1:))
           end if
        end if
     end if
  end do
  close(IOASKI)

  ! now set values of variables in module specfem3D_par_ASKI according to content of Par_file_ASKI

  ! COMPUTE_ASKI_OUTPUT
  call get_value_Par_file_ASKI('COMPUTE_ASKI_OUTPUT',val)
  read(val,*,iostat=ios) COMPUTE_ASKI_OUTPUT
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter logical 'COMPUTE_ASKI_OUTPUT' in '"&
       //'DATA/'//"Par_file_ASKI'")

  if(.not.COMPUTE_ASKI_OUTPUT) then
     if(myrank == 0) then
        call write_ASKI_log('LOG_ASKI_start.txt',"in '"//'DATA/'//&
             "Par_file_ASKI': COMPUTE_ASKI_OUTPUT is .false., so no ASKI output is produced")
        write(IMAIN,*) "in '"//'DATA/'//&
             "Par_file_ASKI': COMPUTE_ASKI_OUTPUT is .false., so no ASKI output is produced"
     end if
     deallocate(ASKI_key_parfile,ASKI_val_parfile)
     return
  end if

  ! ASKI_MAIN_FILE_ONLY
  call get_value_Par_file_ASKI('ASKI_MAIN_FILE_ONLY',val)
  read(val,*,iostat=ios) ASKI_MAIN_FILE_ONLY
  if(ios/=0) call exit_MPI_without_rank("invalid value for logical parameter 'ASKI_MAIN_FILE_ONLY' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! OVERWRITE_ASKI_OUTPUT
  call get_value_Par_file_ASKI('OVERWRITE_ASKI_OUTPUT',val)
  read(val,*,iostat=ios) OVERWRITE_ASKI_OUTPUT
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter logical 'OVERWRITE_ASKI_OUTPUT' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! ASKI_outfile
  call get_value_Par_file_ASKI('ASKI_outfile',ASKI_outfile)

  ! ASKI_output_ID
  call get_value_Par_file_ASKI('ASKI_output_ID',val)
  ASKI_output_ID = val(1:length_ASKI_output_ID)

  ! ASKI_DECONVOLVE_STF
  call get_value_Par_file_ASKI('ASKI_DECONVOLVE_STF',val)
  read(val,*,iostat=ios) ASKI_DECONVOLVE_STF
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DECONVOLVE_STF' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! COMPUTE_ASKI_GREEN_FUNCTION
  call get_value_Par_file_ASKI('COMPUTE_ASKI_GREEN_FUNCTION',val)
  read(val,*,iostat=ios) COMPUTE_ASKI_GREEN_FUNCTION
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'COMPUTE_ASKI_GREEN_FUNCTION' in '"&
       //'DATA/'//"Par_file_ASKI'")

  if(COMPUTE_ASKI_GREEN_FUNCTION) then
     ! ASKI_GREEN_FUNCTION_COMPONENT
     call get_value_Par_file_ASKI('ASKI_GREEN_FUNCTION_COMPONENT',val)
     ASKI_GREEN_FUNCTION_COMPONENT = val(1:length_ASKI_component)
     select case(ASKI_GREEN_FUNCTION_COMPONENT)
     case('N','E','UP')
        ASKI_gf_is_NEZ = .true.
        ASKI_gf_is_XYZ = .false.
        select case(ASKI_GREEN_FUNCTION_COMPONENT)
        case('N'); ASKI_gf_icomp = 1
        case('E'); ASKI_gf_icomp = 2
        case('UP'); ASKI_gf_icomp = 3
        end select
     case('CX','CY','CZ')
        ASKI_gf_is_NEZ = .false.
        ASKI_gf_is_XYZ = .true.
        select case(ASKI_GREEN_FUNCTION_COMPONENT)
        case('CX'); ASKI_gf_icomp = 1
        case('CY'); ASKI_gf_icomp = 2
        case('CZ'); ASKI_gf_icomp = 3
        end select
     case default
     call exit_MPI_without_rank("value for 'ASKI_nchunk' in '"&
          //'DATA/'//"Par_file_ASKI' must be one of 1, 2, 3, 6")
        call exit_MPI_without_rank("value '"//trim(ASKI_GREEN_FUNCTION_COMPONENT)//&
             "' of parameter 'ASKI_GREEN_FUNCTION_COMPONENT' in '"//'DATA/'//"Par_file_ASKI'"//&
             " is not supported. Supported values are: N, E, UP, CX, CY, CZ")
     end select
  else
     ASKI_GREEN_FUNCTION_COMPONENT = ''
  end if

  ! ASKI_df
  call get_value_Par_file_ASKI('ASKI_df',val)
  read(val,*,iostat=ios) ASKI_df
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_df' in '"&
       //'DATA/'//"Par_file_ASKI'")
  if(ASKI_df<0.d0) call exit_MPI_without_rank("value for 'ASKI_df' in '"&
       //'DATA/'//"Par_file_ASKI' must be positive")

  ! ASKI_nf
  call get_value_Par_file_ASKI('ASKI_nf',val)
  read(val,*,iostat=ios) ASKI_nf
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_nf' in '"&
       //'DATA/'//"Par_file_ASKI'")
  if(ASKI_nf<1) call exit_MPI_without_rank("value for 'ASKI_nf' in '"&
       //'DATA/'//"Par_file_ASKI' must be positive")

  allocate(ASKI_jf(ASKI_nf))
  ! ASKI_jf
  call get_value_Par_file_ASKI('ASKI_jf',val)
  read(val,*,iostat=ios) ASKI_jf
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_jf' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! ASKI_DFT_double
  call get_value_Par_file_ASKI('ASKI_DFT_double',val)
  read(val,*,iostat=ios) ASKI_DFT_double
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DFT_double' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! ASKI_DFT_apply_taper
  call get_value_Par_file_ASKI('ASKI_DFT_apply_taper',val)
  read(val,*,iostat=ios) ASKI_DFT_apply_taper
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DFT_apply_taper' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! ASKI_DFT_taper_percentage
  call get_value_Par_file_ASKI('ASKI_DFT_taper_percentage',val)
  read(val,*,iostat=ios) ASKI_DFT_taper_percentage
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DFT_taper_percentage' in '"&
       //'DATA/'//"Par_file_ASKI'")
  if(ASKI_DFT_taper_percentage<0.d0 .or. ASKI_DFT_taper_percentage>1.d0) &
       call exit_MPI_without_rank("value for 'ASKI_DFT_taper_percentage' in '"&
       //'DATA/'//"Par_file_ASKI' must be between 0.0 and 1.0")

  ! ASKI_type_inversion_grid
  call get_value_Par_file_ASKI('ASKI_type_inversion_grid',val)
  read(val,*,iostat=ios) ASKI_type_inversion_grid
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_type_inversion_grid' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! ASKI_nchunk
  call get_value_Par_file_ASKI('ASKI_nchunk',val)
  read(val,*,iostat=ios) ASKI_nchunk
  if(ios/=0) call exit_MPI_without_rank("invalid integer value for parameter 'ASKI_nchunk' in '"&
       //'DATA/'//"Par_file_ASKI'")
  select case(ASKI_nchunk)
  case (1,2,3,6)
     ! OK, do nothing
  case default
     call exit_MPI_without_rank("value for 'ASKI_nchunk' in '"&
          //'DATA/'//"Par_file_ASKI' must be one of 1, 2, 3, 6")
  end select

  ! ASKI_wlat
  call get_value_Par_file_ASKI('ASKI_wlat',val)
  read(val,*,iostat=ios) ASKI_wlat
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_wlat' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! ASKI_wlon
  call get_value_Par_file_ASKI('ASKI_wlon',val)
  read(val,*,iostat=ios) ASKI_wlon
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_wlon' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! ASKI_rot_gamma
  call get_value_Par_file_ASKI('ASKI_rot_gamma',val)
  read(val,*,iostat=ios) ASKI_rot_gamma
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_rot_gamma' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! ASKI_clat
  call get_value_Par_file_ASKI('ASKI_clat',val)
  read(val,*,iostat=ios) ASKI_clat
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_clat' in '"&
       //'DATA/'//"Par_file_ASKI'")
  if(ASKI_clat > 90.0 .or. ASKI_clat < -90.0) &
       call exit_MPI_without_rank("value of parameter 'ASKI_clat' in '"//'DATA/'//"Par_file_ASKI'"//&
       " must be between -90.0 and +90.0")

  ! ASKI_clon
  call get_value_Par_file_ASKI('ASKI_clon',val)
  read(val,*,iostat=ios) ASKI_clon
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_clon' in '"&
       //'DATA/'//"Par_file_ASKI'")
  if(ASKI_clat > 360.0 .or. ASKI_clat < -360.0) &
       call exit_MPI_without_rank("value of parameter 'ASKI_clon' in '"//'DATA/'//"Par_file_ASKI'"//&
       " must be between -360.0 and +360.0")

  ! ASKI_rmax
  call get_value_Par_file_ASKI('ASKI_rmax',val)
  read(val,*,iostat=ios) ASKI_rmax
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_rmax' in '"&
       //'DATA/'//"Par_file_ASKI'")

  ! ASKI_rmin
  call get_value_Par_file_ASKI('ASKI_rmin',val)
  read(val,*,iostat=ios) ASKI_rmin
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_rmin' in '"&
       //'DATA/'//"Par_file_ASKI'")

  if(ASKI_rmax <= ASKI_rmin) call exit_MPI_without_rank("value for 'ASKI_rmax' in '"&
       //'DATA/'//"Par_file_ASKI' must be strictly larger than value for 'ASKI_rmin'")


  if(allocated(ASKI_key_parfile)) deallocate(ASKI_key_parfile)
  if(allocated(ASKI_val_parfile)) deallocate(ASKI_val_parfile)
end subroutine read_Par_file_ASKI
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine get_value_Par_file_ASKI(key,val)
  character(len=*), intent(in) :: key
  character(len=500), intent(out) :: val
  integer :: ipar
  logical :: found
  found = .false.
  do ipar = 1,size(ASKI_key_parfile)
     if(key == ASKI_key_parfile(ipar)) then
        val = ASKI_val_parfile(ipar)
        found = .true.
        exit
     end if
  end do ! ipar
  if(.not.found) call exit_MPI_without_rank("definition of parameter '"//trim(key)//"' not found in '"&
          //'DATA/'//"Par_file_ASKI'")
end subroutine get_value_Par_file_ASKI
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine search_ASKI_wavefield_points_inner_GLL_in_chunks()

  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,R_EARTH_KM
  use specfem_par_crustmantle,only: ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  implicit none

  ! local variables
  integer, dimension(:,:), allocatable :: ASKI_indx_local_tmp
  real(kind=CUSTOM_REAL) :: xstore,ystore,zstore
  double precision :: rgll,xgll,ygll,zgll,rmin,rmax
  integer :: ispec,i,j,k,iglob

  allocate(ASKI_indx_local_tmp((NGLLX-2)*(NGLLY-2)*(NGLLZ-2)*NSPEC_CRUST_MANTLE,4))
  ASKI_indx_local_tmp(:,:) = 0

  ASKI_np_local = 0

  rmin = ASKI_rmin / R_EARTH_KM
  rmax = ASKI_rmax / R_EARTH_KM

  ! loop only on points inside the element. That results in a sufficiently uniform scatter of gridpoints in case of NGLL = 5
  do ispec=1,NSPEC_CRUST_MANTLE
     do k=2,NGLLZ-1
        do j=2,NGLLY-1
           do i=2,NGLLX-1

              iglob = ibool_crust_mantle(i,j,k,ispec)
              xstore = xstore_crust_mantle(iglob)
              ystore = ystore_crust_mantle(iglob)
              zstore = zstore_crust_mantle(iglob)

              ! THESE VALUES ARE SPHERICAL COORDINATES:
              !   xstore = radius in range [0,1]
              !   ystore = theta (radians measured from north pole) in range [0,pi]
              !   zstore = phi (radians longitude measured from 0th meridian to the east) in range [0,2pi]

              rgll = dble(xstore)

              ! exclude this point if it is not in the requested depth range
              if(rgll < rmin .or. rgll > rmax) cycle

              ! if point is in requested depth range, check if it is laterally in the chunk, thererfore
              ! compute unit vector (xgll,ygll,zgll)
              xgll = sin(ystore) ! first factor of x = sin(theta) cos(phi), since ystore=theta
              ygll = xgll*sin(zstore) ! xgll contains the first part of y = sin(theta) sin(phi)  , this is done to save a computation of sin(theta)
              xgll = xgll*cos(zstore) ! last factor of x = sin(theta) cos(phi)
              zgll = cos(ystore)

              if(GLL_point_contained_laterally_in_any_ASKI_output_chunk(xgll,ygll,zgll)) then

                 ! increment index of points found in kernel chunk
                 ASKI_np_local = ASKI_np_local + 1 

                 ! store index of element
                 ASKI_indx_local_tmp(ASKI_np_local,1) = ispec
                 ! store index of x - gridpoint in that element
                 ASKI_indx_local_tmp(ASKI_np_local,2) = i
                 ! store index of y - gridpoint in that element
                 ASKI_indx_local_tmp(ASKI_np_local,3) = j
                 ! store index of z - gridpoint in that element
                 ASKI_indx_local_tmp(ASKI_np_local,4) = k  

              end if ! current point lies in ASKI output volume

           end do ! i
        end do ! j
     enddo ! k
  enddo ! ispec

  if(ASKI_np_local .gt. 0) then
     allocate(ASKI_indx_local(ASKI_np_local,4))
     ASKI_indx_local(:,:) = ASKI_indx_local_tmp(1:ASKI_np_local,:)
  end if

  if(allocated(ASKI_indx_local_tmp)) deallocate(ASKI_indx_local_tmp)

end subroutine search_ASKI_wavefield_points_inner_GLL_in_chunks
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine search_ASKI_wavefield_points_type_invgrid_4()

  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,R_EARTH_KM
  use specfem_par_crustmantle,only: ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  implicit none

  ! local variables
  real(kind=CUSTOM_REAL) :: xstore,ystore,zstore
  double precision :: rgll,xgll,ygll,zgll,rmin,rmax
  integer :: nspec_in_ASKI_volume,ispec,i,j,k,iglob,ip,jspec
  integer, dimension(:), allocatable :: ispec_in_ASKI_volume
  logical :: kji_loop_exit

  allocate(ispec_in_ASKI_volume(NSPEC_CRUST_MANTLE))

  ASKI_np_local = 0

  rmin = ASKI_rmin / R_EARTH_KM
  rmax = ASKI_rmax / R_EARTH_KM

  ! loop on all points inside an element. If any point is contained in the ASKI volume, remember the element index and
  ! later add THE WHOLE element (and all points contained in it) to the list of inversion grid cells / wavefield points
  nspec_in_ASKI_volume = 0
  do ispec=1,NSPEC_CRUST_MANTLE
     kji_loop_exit = .false.
     do k=1,NGLLZ
        do j=1,NGLLY
           do i=1,NGLLX

              iglob = ibool_crust_mantle(i,j,k,ispec)
              xstore = xstore_crust_mantle(iglob)
              ystore = ystore_crust_mantle(iglob)
              zstore = zstore_crust_mantle(iglob)

              ! THESE VALUES ARE SPHERICAL COORDINATES:
              !   xstore = radius in range [0,1]
              !   ystore = theta (radians measured from north pole) in range [0,pi]
              !   zstore = phi (radians longitude measured from 0th meridian to the east) in range [0,2pi]

              rgll = dble(xstore)

              ! exclude this point if it is not in the requested depth range
              if(rgll < rmin .or. rgll > rmax) cycle

              ! if point is in requested depth range, check if it is laterally in the chunk, thererfore
              ! compute unit vector (xgll,ygll,zgll)
              xgll = sin(ystore) ! first factor of x = sin(theta) cos(phi), since ystore=theta
              ygll = xgll*sin(zstore) ! xgll contains the first part of y = sin(theta) sin(phi)  , this is done to save a computation of sin(theta)
              xgll = xgll*cos(zstore) ! last factor of x = sin(theta) cos(phi)
              zgll = cos(ystore)

              if(GLL_point_contained_laterally_in_any_ASKI_output_chunk(xgll,ygll,zgll)) then

                 nspec_in_ASKI_volume = nspec_in_ASKI_volume + 1
                 ispec_in_ASKI_volume(nspec_in_ASKI_volume) = ispec
                 kji_loop_exit = .true.

              end if ! current point lies in ASKI output volume

              if(kji_loop_exit) exit
           end do ! i
           if(kji_loop_exit) exit
        end do ! j
        if(kji_loop_exit) exit
     enddo ! k
  enddo ! ispec

  if(nspec_in_ASKI_volume > 0) then
     ! store wavefield point information at all GLL points of elements found above (i.e. for which at least one point is inside the ASKI output volume
     ASKI_np_local = NGLLX*NGLLY*NGLLZ*nspec_in_ASKI_volume
     allocate(ASKI_indx_local(ASKI_np_local,4))

     ip = 0
     do jspec=1,nspec_in_ASKI_volume; ispec = ispec_in_ASKI_volume(jspec)
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 ! increment point index
                 ip = ip + 1 

                 ! store index of element
                 ASKI_indx_local(ip,1) = ispec
                 ! store index of x - gridpoint in that element
                 ASKI_indx_local(ip,2) = i
                 ! store index of y - gridpoint in that element
                 ASKI_indx_local(ip,3) = j
                 ! store index of z - gridpoint in that element
                 ASKI_indx_local(ip,4) = k  
              end do ! i
           end do ! j
        enddo ! k
     enddo ! ispec

  else ! nspec_in_ASKI_volume > 0
     ASKI_np_local = 0
  end if ! nspec_in_ASKI_volume > 0

  if(allocated(ispec_in_ASKI_volume)) deallocate(ispec_in_ASKI_volume)
end subroutine search_ASKI_wavefield_points_type_invgrid_4
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine define_ASKI_chunk_rotations()
  double precision, parameter :: deg2rad = 0.017453292519943295
  double precision :: ctheta,costheta,cosphi,cosgam,sintheta,sinphi,singam
  double precision, dimension(3,3) :: Mrot_tmp

  ! in case of 6 chunks, there is nothing to rotate
  if(ASKI_nchunk == 6) return

  ! prepare values for the matrix computation below
  cosgam = dcos(dble(ASKI_rot_gamma)*deg2rad)
  singam = dsin(dble(ASKI_rot_gamma)*deg2rad)
  ctheta = 90.d0-dble(ASKI_clat)
  costheta = dcos(ctheta*deg2rad)
  sintheta = dsin(ctheta*deg2rad)
  cosphi = dcos(dble(ASKI_clon)*deg2rad)
  sinphi = dsin(dble(ASKI_clon)*deg2rad)

  ! WHEN PROGRAM COMES TO THIS POINT IN THE CODE, ASKI_nchunk IS ONE OF 1,2,3
  ! Hence, always define ASKI_Mrot_chunk1 in this case:
  ! the (inverse) rotation matrix ASKI_Mrot_chunk1 rotates a point in Cartesian coordinates in the following way:
  !   - the center of a chunk is rotated to the "North Pole" (0,0,1) and the rotation by gamma about the local 
  !     vertical axis of the chunk is reversed
  !   - that means: in this reference frame, x points to local south (i.e. "lat" corresponds to x) and y points to 
  !     local east (i.e. "lon" corresponds to y)

  ASKI_Mrot_chunk1(1,1) = cosgam*costheta*cosphi-singam*sinphi
  ASKI_Mrot_chunk1(1,2) = cosgam*costheta*sinphi+singam*cosphi
  ASKI_Mrot_chunk1(1,3) = -cosgam*sintheta
  ASKI_Mrot_chunk1(2,1) = -singam*costheta*cosphi-cosgam*sinphi
  ASKI_Mrot_chunk1(2,2) = -singam*costheta*sinphi+cosgam*cosphi
  ASKI_Mrot_chunk1(2,3) = singam*sintheta
  ASKI_Mrot_chunk1(3,1) = sintheta*cosphi
  ASKI_Mrot_chunk1(3,2) = sintheta*sinphi
  ASKI_Mrot_chunk1(3,3) = costheta

  if(ASKI_nchunk == 1) then
     ! In case of one chunk, the lateral width of the chunk is accounted for, so define here the parameters
     ! ASKI_xmin_chunk1,ASKI_xmax_chunk1,ASKI_ymin_chunk1,ASKI_ymax_chunk1
     ! NOTE:  in the "north pole" reference frame, x points to local south (i.e. "lat" corresponds to x) and 
     !        y points to local east (i.e. "lon" corresponds to y)
     ASKI_xmin_chunk1 = dtan( - 0.5d0*dble(ASKI_wlat)*deg2rad)
     ASKI_xmax_chunk1 = dtan(   0.5d0*dble(ASKI_wlat)*deg2rad)
     ASKI_ymin_chunk1 = dtan( - 0.5d0*dble(ASKI_wlon)*deg2rad)
     ASKI_ymax_chunk1 = dtan(   0.5d0*dble(ASKI_wlon)*deg2rad)
  end if

  ! In case of 2 or 3 chunks, additionally define the rotation matrix ASKI_Mrot_chunk2
  if(ASKI_nchunk >= 2) then
     ! define ASKI_Mrot_chunk2 as ... (according to the SPECFEM3D 2-chunk simulations and ASKI 2-chunk conventions)
     ! chunk 2 touches chunk 1 on its "left" side (looking from above onto unrotated chunk 1)
     Mrot_tmp(:,:) = 0.d0
     Mrot_tmp(1,1) = 1.d0
     Mrot_tmp(2,3) = 1.d0
     Mrot_tmp(3,2) = -1.d0
     ASKI_Mrot_chunk2(:,:) = matmul(Mrot_tmp,ASKI_Mrot_chunk1)
  end if

  ! In case of 3 chunks, additionally define the rotation matrix ASKI_Mrot_chunk3
  if(ASKI_nchunk == 3) then
     ! define ASKI_Mrot_chunk3 as ... (according to the SPECFEM3D 2-chunk simulations and ASKI 2-chunk conventions)
     ! chunk 3 touches chunk 1 on its "top" side (looking from above onto unrotated chunk 1)
     Mrot_tmp(:,:) = 0.d0
     Mrot_tmp(1,3) = 1.d0
     Mrot_tmp(2,2) = 1.d0
     Mrot_tmp(3,1) = -1.d0
     ASKI_Mrot_chunk3(:,:) = matmul(Mrot_tmp,ASKI_Mrot_chunk1)
  end if
end subroutine define_ASKI_chunk_rotations
!
! ----------------------------------------------------------------------------------------------------------
!
function GLL_point_contained_laterally_in_any_ASKI_output_chunk(xgll,ygll,zgll) result(l)
  double precision :: xgll,ygll,zgll
  logical :: l
  ! local
!
  ! initiate return value to .true.
  ! then simply return in the clause below, when successful
  l  = .true.
!
  ! do these extensive decision clauses in order to reduce the number of calls to function GLL_point_contained_in_this_ASKI_chunk (i.e. do not concetenate those calls by .or. in one if() )
!
  select case(ASKI_nchunk)
  case(1)
     ! check if the point is contained in the only one chunk
     if(GLL_point_contained_laterally_in_this_ASKI_chunk(xgll,ygll,zgll,ASKI_Mrot_chunk1,&
          ASKI_xmin_chunk1,ASKI_xmax_chunk1,ASKI_ymin_chunk1,ASKI_ymax_chunk1)) return
  case(2)
     ! check if the point is contained in the first chunk, assuming width of 90 deg
     if(GLL_point_contained_laterally_in_this_ASKI_chunk(xgll,ygll,zgll,ASKI_Mrot_chunk1,&
          -1.00001d0,1.00001d0,-1.00001d0,1.00001d0)) return
     ! if not, check second chunk, assuming width of 90 deg
     if(GLL_point_contained_laterally_in_this_ASKI_chunk(xgll,ygll,zgll,ASKI_Mrot_chunk2,&
          -1.00001d0,1.00001d0,-1.00001d0,1.00001d0)) return
  case(3)
     ! check if the point is contained in the first chunk, assuming width of 90 deg
     if(GLL_point_contained_laterally_in_this_ASKI_chunk(xgll,ygll,zgll,ASKI_Mrot_chunk1,&
          -1.00001d0,1.00001d0,-1.00001d0,1.00001d0)) return
     ! if not, check second chunk, assuming width of 90 deg
     if(GLL_point_contained_laterally_in_this_ASKI_chunk(xgll,ygll,zgll,ASKI_Mrot_chunk2,&
          -1.00001d0,1.00001d0,-1.00001d0,1.00001d0)) return
     ! if not, check third chunk, assuming width of 90 deg
     if(GLL_point_contained_laterally_in_this_ASKI_chunk(xgll,ygll,zgll,ASKI_Mrot_chunk3,&
          -1.00001d0,1.00001d0,-1.00001d0,1.00001d0)) return
     case(6)
        ! in case of 6 chunks, EVERY point is laterally contained in the ASKI output chunks,
        ! so return l = .true.
        return
  end select
!
  ! if code comes here, no above checks have been successfull, so return .false.
  l = .false.
end function GLL_point_contained_laterally_in_any_ASKI_output_chunk
!
! ----------------------------------------------------------------------------------------------------------
!
function GLL_point_contained_laterally_in_this_ASKI_chunk(xgll,ygll,zgll,Mrot,xmin,xmax,ymin,ymax) result(l)
  double precision :: xgll,ygll,zgll
  double precision :: xmin,xmax,ymin,ymax
  double precision, dimension(3,3) :: Mrot
  logical :: l
  ! local
  double precision :: xtmp,ytmp,ztmp
!
  ! FIRST rotate the incoming GLL point xgll,ygll,zgll (assumed in global Cartesian Coordinates w.r.t EARTH of 
  ! radius 1) to the reference frame, as if the chunk was centered on the north pole
  ! THEN, the rotated point (xtmp,ytmp,ztmp) lies laterally in the chunk if
  !   1) ztmp > 0 
  !   AND if the for the projections onto the tangential plane which touches the chunk at the north pole (xtmp/ztmp , ytmp/ztmp) we have
  !   2) xmin <= xtmp/ztmp <= xmax  AND  ymin <= ytmp/ztmp <= ymax
  ! the implementation below tries to avoid any unnecessary computations and tries to check things
  ! as soon as possible.

  ! initialize return value to false; then simply return below whenever a condition for "contained" is not met
  l = .false.

  ! first compute z-cooordinate of rotated point and check whether it is >0. If not, return
  ztmp = Mrot(3,1)*xgll + Mrot(3,2)*ygll + Mrot(3,3)*zgll
  if(ztmp <= 0.d0) return

  ! compute the rotated x-coordinate and project it onto the tangential plane which touches the chunk at the north pole
  xtmp = (Mrot(1,1)*xgll + Mrot(1,2)*ygll + Mrot(1,3)*zgll)/ztmp
  ! check whether it is in required range between xmin and xmax. If not, return
  if(xtmp < xmin .or. xtmp > xmax) return

  ! now do the same for the y-coordinate
  ytmp = (Mrot(2,1)*xgll + Mrot(2,2)*ygll + Mrot(2,3)*zgll)/ztmp
  if(ytmp < ymin .or. ytmp > ymax) return
!
  ! when the program comes here, the point is laterally contained in this chunk. So return the value l=.true.
  l = .true.
end function GLL_point_contained_laterally_in_this_ASKI_chunk
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine prepare_ASKI_output()

  use specfem_par,only: DT,NSTEP,PI,TWO_PI,myrank

  implicit none

  double precision :: wtaper
  integer :: jt,jf,IOASKI,ios,ntaper
  logical :: file_exists
  character(len=509) :: filename


  if(myrank == 0) then

     do jf=0,ASKI_nf
        if(jf==0) then
           filename = trim(ASKI_outfile)//".main"
        else
           write(filename,"(a,i6.6)") trim(ASKI_outfile)//".jf",ASKI_jf(jf)
        end if

        inquire(file=filename,exist=file_exists)

        ! check here, according to OVERWRITE_ASKI_OUTPUT: is there any conflict with existing files?
        ! if so, abort the current run
        if(file_exists .and. .not.OVERWRITE_ASKI_OUTPUT) &
             call exit_MPI(myrank,"output file '"//trim(filename)//"' "//&
             "already exists and must not be overwritten according to OVERWRITE_ASKI_OUTPUT")

        ! check if files can be opened to write
        ! open file to write with status='unknown' in order not to delete the files when opening
        call get_file_unit_ASKI(IOASKI)
        open(unit=IOASKI, file=filename, status='unknown', access="stream", &
             form="unformatted", action='write',iostat=ios)
        if(ios/=0) call exit_MPI(myrank,"will not be able to write ASKI output file '"//trim(filename)//"'")
        if(file_exists) then
           close(IOASKI,status='keep')
        else
           close(IOASKI,status='delete')
        end if

     end do ! jf        

  end if ! myrank == 0

  ! check if frequency discretization as defined in Par_file_ASKI is valid:
  if(ASKI_df < 0.d0) call exit_MPI_without_rank('ASKI_df in  must not be smaller than zero')
  if(ASKI_nf < 1) call exit_MPI_without_rank('ASKI_nf in Par_file_ASKI must be a strictly positive number')
  if(size(ASKI_jf) /= ASKI_nf) call exit_MPI_without_rank('size(ASKI_jf) must equal ASKI_nf in Par_file_ASKI')

  ! the following parameters and variables are only used for procs which compute any ASKI output at local wavefield points
  if (ASKI_np_local .gt. 0) then

     ! allocate and compute efactors
     allocate(ASKI_efactors_tapered(ASKI_nf,NSTEP))
     do jt = 1,NSTEP
        do jf = 1,ASKI_nf
           ASKI_efactors_tapered(jf,jt) = exp(-(0.d0,1.d0)*TWO_PI*(ASKI_jf(jf)*ASKI_df)*(dble(jt-1)*DT))*DT
        end do ! jf
     end do ! jt

     if(ASKI_DFT_apply_taper) then
        if(ASKI_DFT_taper_percentage<0.d0 .or. ASKI_DFT_taper_percentage>1.d0) &
             call exit_MPI_without_rank('ASKI_DFT_taper_percentage must be between 0.0 and 1.0 in Par_file_ASKI')
        ! set hanning taper parameters, taper the last ASKI_DFT_taper_percentage of the time series
        wtaper = DT*dble(NSTEP-1)*ASKI_DFT_taper_percentage
        ntaper = wtaper/DT

        if(ntaper>0) then
           do jt = NSTEP-ntaper+1,NSTEP
              ASKI_efactors_tapered(:,jt) = &
                   ASKI_efactors_tapered(:,jt)*(0.5d0*(1.d0-cos(PI*DT*dble(NSTEP-jt)/wtaper)))
           end do ! jt
        end if ! ntaper>0
     end if ! ASKI_DFT_apply_taper

     ! allocate for spectra, first rank: 3 underived components, plus 6 strains = 9
     if(ASKI_DFT_double) then
        allocate(ASKI_spectra_local_double(9,ASKI_nf,ASKI_np_local))
        ASKI_spectra_local_double = (0.d0,0.d0)
     else
        allocate(ASKI_spectra_local_single(9,ASKI_nf,ASKI_np_local))
        ASKI_spectra_local_single = (0.0,0.0)
     end if
  end if ! ASKI_np_local .gt. 0

end subroutine prepare_ASKI_output
!
!-------------------------------------------------------------------------------------------
!
subroutine define_ASKI_gf_source()

  use specfem_par,only: myrank,NSOURCES,islice_selected_source,sourcearrays,CUSTOM_REAL,SIZE_REAL,&
       NGLLX,NGLLY,NGLLZ,xigll,yigll,zigll,xi_source,eta_source,gamma_source,nu_source
  use constants, only: R_EARTH,GRAV,RHOAV,PI

  implicit none
  
  character(len=400) :: errstr
  integer :: i,j,k,isource
  double precision :: factor_force_for_one_newton
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  ! FOR WRITING THIS ROUTINE: compare files setup_sources_receivers.f90, from lines 704 on
  ! also compare: SPECFEM CARTESIAN 2.1 src/specfem3d/setup_sources_receivers.f90 from lines 597 on

  ! if NSOURCES is not 1: raise error
  if(NSOURCES /= 1) then
     write(errstr,*) "NSOURCES = ",NSOURCES,"; ASKI Green function computation is only supported for 1 source"
     call exit_MPI_without_rank(trim(errstr))
  end if

  isource = 1

  ! islice_selected_source has size NSOURCES (so only one entry, if code comes here), 
  ! and contains the value myrank for that proc which contains the source
  if(myrank /= islice_selected_source(isource)) return

  ! compute Lagrange polynomials at the source location
  call lagrange_any(xi_source(isource),NGLLX,xigll,hxis,hpxis)
  call lagrange_any(eta_source(isource),NGLLY,yigll,hetas,hpetas)
  call lagrange_any(gamma_source(isource),NGLLZ,zigll,hgammas,hpgammas)

  ! reset the sourcearray for the first (and only) source before defining (some of) its values anew
  sourcearrays(:,:,:,:,isource) = 0._CUSTOM_REAL

  ! calculate source array entries for interpolated single force location
  ! NOTE THAT THE APPLIED FORCE MUST BE A UNIT FORCE. WE NEED TO NON-DIMENSIONALIIZE 1 Newton HERE
  ! BY APPLYING AN ADDITIONAL NON-DIMENSIONALIZATION FACTOR BELOW!
  ! for the regular setup of sourcearrays, the non-dimensionalization of the 
  ! moment tensor is done at the very end of src/specfem3d/get_cmt.f90.
  ! So if we want to have a force of 1 Newton here, divide by 
  ! (RHOAV * R_EARTH**3) * R_EARTH * (PI * GRAV * RHOAV)
  ! the first factor is for Kg (which is density * volume), the second factor is for m and the last factor is for 1/s^2 
  ! (compare SPECFEM3D_GLOBE manual, time dimensionalization)
  factor_force_for_one_newton = 1.d0 / (RHOAV**2 * R_EARTH**4 * PI * GRAV)

  ! presumably there is a bug (?) in form of a nearly planar wave arriving from the deep model
  ! (observed in a 1chunk simulation). The amplitude and waveform of that wave is INDEPENDENT of 
  ! the exciting seismic source. When applying a single force with 1 N, the resulting aplitudes
  ! of the seismic waves are VERY small, whereby the strange artefact "wave" becomes highly 
  ! significant and contaminates the actual waveforms tremendously.
  ! (Florian Schumacher is still unsure whether this is a "bug" which is known to the developer of SPECFEM3D)
  ! In order to account for this phenomenom, we boost the amplitude from 1 Newton to 1e12 Newton.
  ! We do NOT divide by 1e12 before writing ASKI output to file, but instead apply a unit factor of 1e-12 
  ! in the ASKI main package for kernel green tensor output. 
  ! Since the amplitude of the strange artefact signal 
  ! is independent of the applied seismic source, this removes the artefact from the waveforms.
  factor_force_for_one_newton = factor_force_for_one_newton * 1.d12 

  if(ASKI_gf_is_NEZ) then
     do k=1,NGLLZ
        do j=1,NGLLY
           do i=1,NGLLX
              if(CUSTOM_REAL == SIZE_REAL) then
                 ! for ASKI_gf_icomp (first dimension of array nu_source) : E/N/Z = 1/2/3
                 sourcearrays(:,i,j,k,isource) = real( factor_force_for_one_newton * &
                      hxis(i)*hetas(j)*hgammas(k)*nu_source(ASKI_gf_icomp,:,isource) )
              else
                 sourcearrays(:,i,j,k,isource) = factor_force_for_one_newton * &
                      hxis(i)*hetas(j)*hgammas(k)*nu_source(ASKI_gf_icomp,:,isource)
              end if
           end do ! i
        end do ! j
     end do ! k
  end if ! ASKI_gf_is_NEZ

  if(ASKI_gf_is_XYZ) then
     do k=1,NGLLZ
        do j=1,NGLLY
           do i=1,NGLLX
              if(CUSTOM_REAL == SIZE_REAL) then
                 ! in case of Cartesian X,Y,Z components, only one entry in the first dimension of sourcearrays is /= 0
                 sourcearrays(ASKI_gf_icomp,i,j,k,isource) = real( factor_force_for_one_newton * &
                      hxis(i)*hetas(j)*hgammas(k) )
              else
                 sourcearrays(ASKI_gf_icomp,i,j,k,isource) = factor_force_for_one_newton * &
                      hxis(i)*hetas(j)*hgammas(k)
              end if
           end do ! i
        end do ! j
     end do ! k
  end if ! ASKI_gf_is_XYZ

end subroutine define_ASKI_gf_source
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_main_file()

  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,myrank,NPROCTOT
  use specfem_par_crustmantle,only: ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle,&
       rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle,& ! assume isotropy, i.e. use only arrays for isotroppic elements (vertical properties kappav, muv)
       xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
       etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
       gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle
  use constants, only: R_EARTH_KM,GRAV,RHOAV,PI

  implicit none

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: xyz_local,xyz
  real(kind=CUSTOM_REAL) :: r,theta,phi
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: model_local,model
  double precision :: dimens_veloc
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: jacob_local,jacob
  integer :: ip,iglob,ispec,i,j,k
  integer :: iproc,specfem_version
  integer :: IOASKI
  logical :: file_exists
  character (len=7) :: open_status
  integer, dimension(:,:), allocatable :: neighbours
  integer :: nnb,ncell
  integer, parameter :: itag = 100

  ! DIMENSIONALIZE THE RETURN VALUES:
  ! constants GRAV,RHOAV are defined in constants.h, hence contained in module constants
  dimens_veloc = R_EARTH_KM*dsqrt(PI*GRAV*RHOAV) ! time scaling (s^{-1}) is done with dsqrt(PI*GRAV*RHOAV)

  ! first, locally collect all coordinates of wavefield points and compute kernel reference model
  if(ASKI_np_local > 0) then
     allocate(xyz_local(ASKI_np_local,3),model_local(ASKI_np_local,3))

     do ip = 1,ASKI_np_local
        ! set ispec,i,j,k
        ispec = ASKI_indx_local(ip,1)
        i = ASKI_indx_local(ip,2)
        j = ASKI_indx_local(ip,3)
        k = ASKI_indx_local(ip,4)
        iglob = ibool_crust_mantle(i,j,k,ispec)

        ! r, theta, phi were stored in xstore,ystore,zstore (r non-dimensionalized  in m/R_EARTH_m)
        r = R_EARTH_KM * xstore_crust_mantle(iglob)
        theta = ystore_crust_mantle(iglob)
        phi = zstore_crust_mantle(iglob)

        ! x, y, z
        xyz_local(ip,1) = r * sin(theta) * cos(phi)
        xyz_local(ip,2) = r * sin(theta) * sin(phi)
        xyz_local(ip,3) = r * cos(theta)

        ! DIMENSIONALIZE THE RETURN VALUES:

        ! density rho was non-dimensionalized in  (Kg/m^3)/RHOAV
        model_local(ip,1) = (RHOAV/1000.0d0)* rhostore_crust_mantle(i,j,k,ispec)  ! now [g/cm^3]

        ! velocities were non-dimensionalized in  (m/s)/(R_EARTH_m*sqrt(PI*GRAV*RHOAV))
        model_local(ip,2) = dimens_veloc * sqrt((kappavstore_crust_mantle(i,j,k,ispec) + &
             4.*muvstore_crust_mantle(i,j,k,ispec)/3.)/rhostore_crust_mantle(i,j,k,ispec))  ! now [km/s]
        ! vs
        model_local(ip,3) = dimens_veloc * sqrt(muvstore_crust_mantle(i,j,k,ispec)/rhostore_crust_mantle(i,j,k,ispec))  ! now [km/s]

     end do ! ip

     ! in case of specfem3dInversionGrid, additionally compute jacobian at all GLL points
     select case(ASKI_type_inversion_grid)
     case(4)
        allocate(jacob_local(ASKI_np_local))
        do ip = 1,ASKI_np_local
           ! set ispec,i,j,k of this point
           ispec = ASKI_indx_local(ip,1)
           i = ASKI_indx_local(ip,2)
           j = ASKI_indx_local(ip,3)
           k = ASKI_indx_local(ip,4)

           ! get values from xix... arrays
           xixl = xix_crust_mantle(i,j,k,ispec)
           xiyl = xiy_crust_mantle(i,j,k,ispec)
           xizl = xiz_crust_mantle(i,j,k,ispec)
           etaxl = etax_crust_mantle(i,j,k,ispec)
           etayl = etay_crust_mantle(i,j,k,ispec)
           etazl = etaz_crust_mantle(i,j,k,ispec)
           gammaxl = gammax_crust_mantle(i,j,k,ispec)
           gammayl = gammay_crust_mantle(i,j,k,ispec)
           gammazl = gammaz_crust_mantle(i,j,k,ispec)

           ! jacobian
           ! NEED TO "DIMENSIONALIZE" JACOBIAN BY R_EARTH_KM^3, since the values x,y,z 
           ! as functions of (xi,eta,gamma) in the code here refer to values in an 
           ! earth with radius 1.0 . The actual coordinate functions x,y,z need to be scaled
           ! by R_EARTH_KM in order to match the ASKI wavefield point definition.
           ! But this means that in the products below of the form  (xixl*(etayl*gammazl-etazl*gammayl)
           ! there is a factor 1./R_EARTH_KM in every term xix,etay,gammz, .. etc. (it's 1/R_EARTH_KM 
           ! in this formula, because xix is the derivatives of the the inverse transformation xi(x,y,z) 
           ! from x,y,z to xi (same for eta,gamma) ). This is the reason for the factor R_EARTH_KM**3 in the next formula.
           jacob_local(ip) = R_EARTH_KM**3 / & ! 1.0_CUSTOM_REAL / &
                (xixl*(etayl*gammazl-etazl*gammayl) &
                - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                + xizl*(etaxl*gammayl-etayl*gammaxl))
        end do ! ip
     end select

  end if ! ASKI_np_local > 0

  call synchronize_all()

  ! gather wavefield points on proc 0
  if(myrank == 0) then
     allocate(xyz(sum(ASKI_np_local_all),3))
     ip = 0

     ! this is me, rank 0
     if(ASKI_np_local > 0) then
        xyz(1:ASKI_np_local,:) = xyz_local(:,:)
        ip = ip + ASKI_np_local
     end if
     do iproc = 1,NPROCTOT-1
        if(ASKI_np_local_all(iproc+1) > 0) then
           call recv_cr(xyz(ip+1:ip+ASKI_np_local_all(iproc+1),:), ASKI_np_local_all(iproc+1)*3, iproc, itag)
           ip = ip + ASKI_np_local_all(iproc+1)
        end if
     end do ! iproc

  else ! myrank == 0

     if(ASKI_np_local > 0) call send_cr(xyz_local, ASKI_np_local*3, 0, itag)

  end if ! myrank == 0

  call synchronize_all()

  select case(ASKI_type_inversion_grid)
  case(4)
     ! gather jacobian on proc 0 and compute neighbours
     if(myrank == 0) then
        allocate(jacob(sum(ASKI_np_local_all)))
        ip = 0

        ! this is me, rank 0
        if(ASKI_np_local > 0) then
           jacob(1:ASKI_np_local) = jacob_local(:)
           ip = ip + ASKI_np_local
        end if
        do iproc = 1,NPROCTOT-1
           if(ASKI_np_local_all(iproc+1) > 0) then
              call recv_cr(jacob(ip+1:ip+ASKI_np_local_all(iproc+1)), ASKI_np_local_all(iproc+1), iproc, itag)
              ip = ip + ASKI_np_local_all(iproc+1)
           end if
        end do ! iproc

        ! find neighbouring elements
        ncell = sum(ASKI_np_local_all)/(NGLLX*NGLLY*NGLLZ)
        allocate(neighbours(7,ncell)); neighbours = 0
        call find_ASKI_neighbours_type_invgrid_4(neighbours,ncell,xyz,sum(ASKI_np_local_all))

     else ! myrank == 0

        if(ASKI_np_local > 0) &
             call send_cr(jacob_local, ASKI_np_local, 0, itag)

     end if ! myrank == 0
  end select ! ASKI_type_inversion_grid

  ! gather model on proc 0
  if(myrank == 0) then
     allocate(model(sum(ASKI_np_local_all),3))
     ip = 0

     ! this is me, rank 0
     if(ASKI_np_local > 0) then
        model(1:ASKI_np_local,:) = model_local(:,:)
        ip = ip + ASKI_np_local
     end if
     do iproc = 1,NPROCTOT-1
        if(ASKI_np_local_all(iproc+1) > 0) then
           call recv_cr(model(ip+1:ip+ASKI_np_local_all(iproc+1),:), ASKI_np_local_all(iproc+1)*3, iproc, itag)
           ip = ip + ASKI_np_local_all(iproc+1)
        end if
     end do ! iproc

  else ! myrank == 0

     if(ASKI_np_local > 0) &
          call send_cr(model_local, ASKI_np_local*3, 0, itag)

  end if ! myrank == 0

  call synchronize_all()

  if(myrank == 0) then
     ! do not worry about whether to overwrite or not, if program comes here, it already has been checked 
     ! if the value of OVERWRITE_ASKI_OUTPUT is in conflict with existing files
     inquire(file=trim(ASKI_outfile)//".main",exist=file_exists)
     if(file_exists) then
        open_status = 'replace'
     else
        open_status = 'new'
     end if

     ! specfem_version
     ! 1 = SPECFEM3D_GLOBE, 2 = SPECFEM3D_Cartesian
     specfem_version = 1

     ! master proc writes out file containing general info and wavefield points and kernel reference model
     ! also write NPROCTOT to file, as a safety feature in order to assure that wavefield points (written only in ASKI_main_file
     ! and kernel values  (written only in kernel files) are not confused
     call get_file_unit_ASKI(IOASKI)
     open(unit=IOASKI, file=trim(ASKI_outfile)//".main", status=open_status, access="stream", form="unformatted")
     write(IOASKI) specfem_version,length_ASKI_output_ID,ASKI_output_ID,NPROCTOT,ASKI_type_inversion_grid,&
          sum(ASKI_np_local_all),ASKI_df,ASKI_nf
     write(IOASKI) ASKI_jf
     write(IOASKI) xyz(:,1)
     write(IOASKI) xyz(:,2)
     write(IOASKI) xyz(:,3)
     write(IOASKI) model(:,1)
     write(IOASKI) model(:,2)
     write(IOASKI) model(:,3)
     select case(ASKI_type_inversion_grid)
     case(4)
        write(IOASKI) NGLLX,NGLLY,NGLLZ
        write(IOASKI) jacob
        write(IOASKI) ncell,sum(neighbours(1,:))+ncell !sum(neighbours(1,:))+ncell is the total number of integers following now, before model values
        do ispec = 1,ncell
           nnb = neighbours(1,ispec)
           write(IOASKI) nnb
           if(nnb > 0) write(IOASKI) neighbours(2:1+nnb,ispec)
        end do
     end select
     close(IOASKI)
  end if

  ! deallocate local stuff
  if(allocated(model)) deallocate(model)
  if(allocated(model_local)) deallocate(model_local)
  if(allocated(xyz)) deallocate(xyz)
  if(allocated(xyz_local)) deallocate(xyz_local)
  if(allocated(jacob)) deallocate(jacob)
  if(allocated(jacob_local)) deallocate(jacob_local)
  if(allocated(neighbours)) deallocate(neighbours)

end subroutine write_ASKI_main_file
!
!-------------------------------------------------------------------------------------------
!
subroutine find_ASKI_neighbours_type_invgrid_4(neighbours,ncell,xyz,np)

  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer :: ncell,np
  integer, dimension(7,ncell) :: neighbours
  real(kind=CUSTOM_REAL), dimension(np,3) :: xyz
  ! local
  integer :: icell,jcell,iface,nnb,cell_shift,np_cell,&
       nb1_shift,nb2_shift,nb3_shift,nb4_shift,nb5_shift,nb6_shift
  real(kind=CUSTOM_REAL), dimension(6*ncell) :: xnb,ynb,znb
  real(kind=CUSTOM_REAL), dimension(3) :: v
  real :: eps

  if(mod(NGLLX,2)/=1 .or. mod(NGLLY,2)/=1 .or. mod(NGLLZ,2)/=1) &
       call exit_MPI_without_rank("Neighbour search for specfem3dInversionGrid for "//&
       "ASKI only supported for odd NGLLX,NGLLY,NGLLZ")

  neighbours = 0

  ! loop on all cells and store their 6 face-center points in arrays xnb,ynb,znb

  nb1_shift = NGLLX*NGLLY*(NGLLZ-1)/2 + (NGLLX+1)/2 ! point-index of ymin-face relative to cell
  nb2_shift = NGLLX*NGLLY*(NGLLZ-1)/2 + NGLLX*(NGLLY+1)/2 ! xmax-face
  nb3_shift = NGLLX*NGLLY*(NGLLZ-1)/2 + NGLLX*(NGLLY-1) + (NGLLX+1)/2 ! ymax-face
  nb4_shift = NGLLX*NGLLY*(NGLLZ-1)/2 + NGLLX*(NGLLY-1)/2 + 1 ! xmin-face
  nb5_shift = NGLLX*(NGLLY-1)/2 + (NGLLX+1)/2 ! zmin-face
  nb6_shift = NGLLX*NGLLY*(NGLLZ-1) + NGLLX*(NGLLY-1)/2 + (NGLLX+1)/2 ! zmax-face
  np_cell = NGLLX*NGLLY*NGLLZ

  cell_shift = 0
  do icell = 1,ncell
     ! coordinates of first face-center
     xnb((icell-1)*6+1) = xyz(cell_shift+nb1_shift,1)
     ynb((icell-1)*6+1) = xyz(cell_shift+nb1_shift,2)
     znb((icell-1)*6+1) = xyz(cell_shift+nb1_shift,3)

     ! coordinates of second face-center
     xnb((icell-1)*6+2) = xyz(cell_shift+nb2_shift,1)
     ynb((icell-1)*6+2) = xyz(cell_shift+nb2_shift,2)
     znb((icell-1)*6+2) = xyz(cell_shift+nb2_shift,3)

     ! coordinates of third face-center
     xnb((icell-1)*6+3) = xyz(cell_shift+nb3_shift,1)
     ynb((icell-1)*6+3) = xyz(cell_shift+nb3_shift,2)
     znb((icell-1)*6+3) = xyz(cell_shift+nb3_shift,3)

     ! coordinates of fourth face-center
     xnb((icell-1)*6+4) = xyz(cell_shift+nb4_shift,1)
     ynb((icell-1)*6+4) = xyz(cell_shift+nb4_shift,2)
     znb((icell-1)*6+4) = xyz(cell_shift+nb4_shift,3)

     ! coordinates of fifth face-center
     xnb((icell-1)*6+5) = xyz(cell_shift+nb5_shift,1)
     ynb((icell-1)*6+5) = xyz(cell_shift+nb5_shift,2)
     znb((icell-1)*6+5) = xyz(cell_shift+nb5_shift,3)

     ! coordinates of sixth face-center
     xnb((icell-1)*6+6) = xyz(cell_shift+nb6_shift,1)
     ynb((icell-1)*6+6) = xyz(cell_shift+nb6_shift,2)
     znb((icell-1)*6+6) = xyz(cell_shift+nb6_shift,3)

     ! go to next NGLLX*NGLLY*NGLLZ points in array xyz, which form the next cell
     cell_shift = cell_shift + np_cell
  end do ! icell


  ! again loop on all cells and check for all faces if any other face-cell has the same coordinates (i.e. the faces are a match)

  do icell = 1,ncell
     ! define threashold for this cell to test distances to other face-centers as the minimum over all dimensions of
     ! face width per number of GLL points (multiplied scaled by 1e-3)
     v(1) = xnb((icell-1)*6+1)-xnb((icell-1)*6+3) ! vector v contains xmin-face-center minus xmax-face-center
     v(2) = ynb((icell-1)*6+1)-ynb((icell-1)*6+3)
     v(3) = znb((icell-1)*6+1)-znb((icell-1)*6+3)
     eps = sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )/real(NGLLX)
     v(1) = xnb((icell-1)*6+2)-xnb((icell-1)*6+4) ! vector v contains ymin-face-center minus ymax-face-center
     v(2) = ynb((icell-1)*6+2)-ynb((icell-1)*6+4)
     v(3) = znb((icell-1)*6+2)-znb((icell-1)*6+4)
     eps = min(eps, sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )/real(NGLLY) )
     v(1) = xnb((icell-1)*6+5)-xnb((icell-1)*6+6) ! vector v contains zmin-face-center minus zmax-face-center
     v(2) = ynb((icell-1)*6+5)-ynb((icell-1)*6+6)
     v(3) = znb((icell-1)*6+5)-znb((icell-1)*6+6)
     eps = min(eps, sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )/real(NGLLZ) )
     eps = eps * 0.001

     ! loop on all faces of cell icell
     do iface = 1,6
        ! v contains current face-center
        v(1) = xnb((icell-1)*6+iface)
        v(2) = ynb((icell-1)*6+iface)
        v(3) = znb((icell-1)*6+iface)

        ! loop on all other cells now and compare all face centers of those faces to current face center v
        do jcell = 1,ncell
           if(jcell==icell) cycle ! don't compare to own face-centers

           nnb = count( sqrt( (xnb((jcell-1)*6+1:jcell*6)-v(1))**2 + &
                (ynb((jcell-1)*6+1:jcell*6)-v(2))**2 + &
                (znb((jcell-1)*6+1:jcell*6)-v(3))**2) < eps)

           if(nnb>1) then
              write(*,*) "ASKI ERROR: for ",iface,"'th face of cell ",icell," there are ",nnb,&
                   " neighbouring faces in cell ",jcell
              call exit_MPI_without_rank("ASKI ERROR finding neighbours for specfem3dInversionGrid: "//&
                   "more than 1 face neighbours for one face")
           end if
           if(nnb==1) then
              ! increase the count of found neighbours
              neighbours(1,icell) = neighbours(1,icell) + 1
              ! store cell index of this found neighbour
              neighbours(neighbours(1,icell)+1,icell) = jcell
              exit ! loop over jcell, goto next face
           end if ! nnb==1

        end do ! jcell
     end do ! iface
  end do ! icell

end subroutine find_ASKI_neighbours_type_invgrid_4
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_output()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ, R_EARTH_KM,GRAV,RHOAV,PI
  use specfem_par,only: it,NSTEP,hprime_xx,hprime_yy,hprime_zz,scale_veloc,scale_t_inv
  use specfem_par_crustmantle,only: xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
       etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
       gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,&
       veloc_crust_mantle,ibool_crust_mantle

  implicit none

  integer :: ip,ispec,i,j,k,l,iglob,jf
  real(kind=CUSTOM_REAL) :: xix_ip,xiy_ip,xiz_ip, &
       etax_ip,etay_ip,etaz_ip, &
       gammax_ip,gammay_ip,gammaz_ip
  real(kind=CUSTOM_REAL) :: sumxxi,sumyxi,sumzxi, &
       sumxeta,sumyeta,sumzeta, &
       sumxgamma,sumygamma,sumzgamma
  real(kind=CUSTOM_REAL) :: ux,uy,uz,&
       uxdx,uxdy,uxdz, &
       uydx,uydy,uydz, &
       uzdx,uzdy,uzdz, &
       e11,e22,e33,e23,e13,e12
  complex(kind=kind(1.d0)) :: efactor_tapered
  double precision :: dimensionalize_wavefield,dimensionalize_strain

  if(.not.COMPUTE_ASKI_OUTPUT) return

  ! fourier transform to spectra only, if there are local wavefield points
  if(ASKI_np_local > 0) then

     do ip=1,ASKI_np_local

        ! set ispec,i,j,k
        ispec = ASKI_indx_local(ip,1)
        i = ASKI_indx_local(ip,2)
        j = ASKI_indx_local(ip,3)
        k = ASKI_indx_local(ip,4)

        ! calculate derivatives

        ! set derivatives of xi,eta,gamma
        xix_ip = xix_crust_mantle(i,j,k,ispec)
        xiy_ip = xiy_crust_mantle(i,j,k,ispec)
        xiz_ip = xiz_crust_mantle(i,j,k,ispec)
        etax_ip = etax_crust_mantle(i,j,k,ispec)
        etay_ip = etay_crust_mantle(i,j,k,ispec)
        etaz_ip = etaz_crust_mantle(i,j,k,ispec)
        gammax_ip = gammax_crust_mantle(i,j,k,ispec)
        gammay_ip = gammay_crust_mantle(i,j,k,ispec)
        gammaz_ip = gammaz_crust_mantle(i,j,k,ispec)

        ! calculate sum with h'_i(xi_l)=hprime_xx(i,l)
        sumxxi = 0.0
        sumyxi = 0.0
        sumzxi = 0.0

        do l=1,NGLLX

           iglob = ibool_crust_mantle(l,j,k,ispec)
           sumxxi = sumxxi + veloc_crust_mantle(1,iglob)*hprime_xx(i,l)
           sumyxi = sumyxi + veloc_crust_mantle(2,iglob)*hprime_xx(i,l)
           sumzxi = sumzxi + veloc_crust_mantle(3,iglob)*hprime_xx(i,l)

        end do ! l

        ! calculate sum with h'_j(eta_l)=hprime_yy(j,l)
        sumxeta = 0.0
        sumyeta = 0.0
        sumzeta = 0.0

        do l=1,NGLLY

           iglob = ibool_crust_mantle(i,l,k,ispec)
           sumxeta = sumxeta + veloc_crust_mantle(1,iglob)*hprime_yy(j,l)
           sumyeta = sumyeta + veloc_crust_mantle(2,iglob)*hprime_yy(j,l)
           sumzeta = sumzeta + veloc_crust_mantle(3,iglob)*hprime_yy(j,l)

        end do ! l

        ! calculate sum with h'_k(gamma_l)=hprime_zz(k,l)
        sumxgamma = 0.0
        sumygamma = 0.0
        sumzgamma = 0.0

        do l=1,NGLLZ

           iglob = ibool_crust_mantle(i,j,l,ispec)
           sumxgamma = sumxgamma + veloc_crust_mantle(1,iglob)*hprime_zz(k,l)
           sumygamma = sumygamma + veloc_crust_mantle(2,iglob)*hprime_zz(k,l)
           sumzgamma = sumzgamma + veloc_crust_mantle(3,iglob)*hprime_zz(k,l)

        end do ! l

        ! now calculate the derivative of veloc w.r.t. x, y and z with help of the sums calculated above
        ! call it u, since in the context of ASKI, this velocity field w.r.t Heaviside excitation is 
        ! interpreted as the DISPLACEMENT field w.r.t Delta-impulse excitation!

        ! derivative by x
        uxdx = xix_ip*sumxxi + etax_ip*sumxeta + gammax_ip*sumxgamma
        uydx = xix_ip*sumyxi + etax_ip*sumyeta + gammax_ip*sumygamma
        uzdx = xix_ip*sumzxi + etax_ip*sumzeta + gammax_ip*sumzgamma

        ! derivative by y
        uxdy = xiy_ip*sumxxi + etay_ip*sumxeta + gammay_ip*sumxgamma
        uydy = xiy_ip*sumyxi + etay_ip*sumyeta + gammay_ip*sumygamma
        uzdy = xiy_ip*sumzxi + etay_ip*sumzeta + gammay_ip*sumzgamma

        ! derivative by z
        uxdz = xiz_ip*sumxxi + etaz_ip*sumxeta + gammaz_ip*sumxgamma
        uydz = xiz_ip*sumyxi + etaz_ip*sumyeta + gammaz_ip*sumygamma
        uzdz = xiz_ip*sumzxi + etaz_ip*sumzeta + gammaz_ip*sumzgamma

        ! store underived velocity wavefield
        ! call it u, since in the context of ASKI, this velocity field w.r.t Heaviside excitation is 
        ! interpreted as the DISPLACEMENT field w.r.t Delta-impulse excitation!
        iglob = ibool_crust_mantle(i,j,k,ispec)
        ux = veloc_crust_mantle(1,iglob)
        uy = veloc_crust_mantle(2,iglob)
        uz = veloc_crust_mantle(3,iglob)

        ! strains
        e11 = uxdx
        e22 = uydy
        e33 = uzdz
        e23 = 0.5*(uydz + uzdy)
        e13 = 0.5*(uxdz + uzdx)
        e12 = 0.5*(uxdy + uydx)

        ! conduct DFT
        if(ASKI_DFT_double) then

           do jf = 1,ASKI_nf
              ! reduce array access by reading from array efactors_tapered only once
              efactor_tapered = ASKI_efactors_tapered(jf,it)

              ! underived wavefiled spectra
              ASKI_spectra_local_double(1,jf,ip) = ASKI_spectra_local_double(1,jf,ip) + ux*efactor_tapered
              ASKI_spectra_local_double(2,jf,ip) = ASKI_spectra_local_double(2,jf,ip) + uy*efactor_tapered
              ASKI_spectra_local_double(3,jf,ip) = ASKI_spectra_local_double(3,jf,ip) + uz*efactor_tapered

              ! strain spectra in order e11,e22,e33,e23,e13,e12
              ASKI_spectra_local_double(4,jf,ip) = ASKI_spectra_local_double(4,jf,ip) + e11*efactor_tapered
              ASKI_spectra_local_double(5,jf,ip) = ASKI_spectra_local_double(5,jf,ip) + e22*efactor_tapered
              ASKI_spectra_local_double(6,jf,ip) = ASKI_spectra_local_double(6,jf,ip) + e33*efactor_tapered
              ASKI_spectra_local_double(7,jf,ip) = ASKI_spectra_local_double(7,jf,ip) + e23*efactor_tapered
              ASKI_spectra_local_double(8,jf,ip) = ASKI_spectra_local_double(8,jf,ip) + e13*efactor_tapered
              ASKI_spectra_local_double(9,jf,ip) = ASKI_spectra_local_double(9,jf,ip) + e12*efactor_tapered
           end do ! jf

        else ! ASKI_DFT_double

           do jf = 1,ASKI_nf
              ! reduce array access by reading from array efactors_tapered only once
              efactor_tapered = ASKI_efactors_tapered(jf,it)

              ! underived wavefiled spectra
              ASKI_spectra_local_single(1,jf,ip) = ASKI_spectra_local_single(1,jf,ip) + ux*efactor_tapered
              ASKI_spectra_local_single(2,jf,ip) = ASKI_spectra_local_single(2,jf,ip) + uy*efactor_tapered
              ASKI_spectra_local_single(3,jf,ip) = ASKI_spectra_local_single(3,jf,ip) + uz*efactor_tapered

              ! strain spectra in order e11,e22,e33,e23,e13,e12
              ASKI_spectra_local_single(4,jf,ip) = ASKI_spectra_local_single(4,jf,ip) + e11*efactor_tapered
              ASKI_spectra_local_single(5,jf,ip) = ASKI_spectra_local_single(5,jf,ip) + e22*efactor_tapered
              ASKI_spectra_local_single(6,jf,ip) = ASKI_spectra_local_single(6,jf,ip) + e33*efactor_tapered
              ASKI_spectra_local_single(7,jf,ip) = ASKI_spectra_local_single(7,jf,ip) + e23*efactor_tapered
              ASKI_spectra_local_single(8,jf,ip) = ASKI_spectra_local_single(8,jf,ip) + e13*efactor_tapered
              ASKI_spectra_local_single(9,jf,ip) = ASKI_spectra_local_single(9,jf,ip) + e12*efactor_tapered
           end do ! jf

        end if ! ASKI_DFT_double

     end do ! ip

  end if ! ASKI_np_local > 0

  if (it == NSTEP) then

     if(ASKI_np_local > 0) then
        ! dimensionalize the output wavefield spectra!
        ! since particle VELOCITIES were used as wavefield output, dimensionalize here to velocity in m/s.
        ! dimensionalize the strains to 1/s.
        dimensionalize_wavefield = scale_veloc ! scale_veloc = R_EARTH * dsqrt(PI*GRAV*RHOAV)  ->  see subroutine prepare_timerun_constants()
        dimensionalize_strain = scale_t_inv    ! scale_t_inv = dsqrt(PI*GRAV*RHOAV)

        ! for Green function simulation, one could do the back-scaling by 1.d-12 here (compensating for the boosted
        ! 1.d12 N single fource source). 
        ! HOWEVER, this results in very small amplitudes, which might get to the edge of single precision representation 
        ! (around 1.18e-39). HENCE, keep the factor of 1.d12, resulting in pico meter wavefields for 1N source.
        ! This is compensated later on when computing kernels, setting the unit factors of specfem3d_kernel_green_tensor
        ! objects to 1.e-12 for both, underived wavefield and strains.

        !
        ! if(COMPUTE_ASKI_GREEN_FUNCTION) then
        !    dimensionalize_wavefield = dimensionalize_wavefield * 1.d-12
        !    dimensionalize_strain = dimensionalize_strain * 1.d-12
        ! end if

        if(ASKI_DFT_double) then
           ASKI_spectra_local_double(1:3,:,:) = ASKI_spectra_local_double(1:3,:,:) * dimensionalize_wavefield
           ASKI_spectra_local_double(4:9,:,:) = ASKI_spectra_local_double(4:9,:,:) * dimensionalize_strain
        else
           ! is it better to do multiplication in double precision, i.e. not to do real(dimensionalize...) ?! 
           ! It is anyway recommended to use ASKI_DFT_double = .true. !! (so this clause is not recommended)
           ASKI_spectra_local_single(1:3,:,:) = ASKI_spectra_local_single(1:3,:,:) * real(dimensionalize_wavefield)
           ASKI_spectra_local_single(4:9,:,:) = ASKI_spectra_local_single(4:9,:,:) * real(dimensionalize_strain)
        end if

        if(ASKI_DECONVOLVE_STF) call deconvolve_diff_stf_from_ASKI_output()
     end if ! ASKI_np_local > 0

     call write_ASKI_output_files()

     ! deallocate everything, simulation is over
     if(allocated(ASKI_np_local_all)) deallocate(ASKI_np_local_all)
     if(allocated(ASKI_indx_local)) deallocate(ASKI_indx_local)
     if(allocated(ASKI_efactors_tapered)) deallocate(ASKI_efactors_tapered)
     if(allocated(ASKI_spectra_local_double)) deallocate(ASKI_spectra_local_double)
     if(allocated(ASKI_spectra_local_single)) deallocate(ASKI_spectra_local_single)
     if(allocated(ASKI_jf)) deallocate(ASKI_jf)
  end if ! it == NSTEP

end subroutine write_ASKI_output
!
!-------------------------------------------------------------------------------------------
!
subroutine deconvolve_diff_stf_from_ASKI_output()
  use specfem_par,only: myrank,NSTEP,DT,t0,tshift_cmt,hdur_gaussian
  use constants,only: IMAIN
  implicit none
  double precision, external :: comp_source_time_function
  complex(kind=kind(1.d0)), dimension(:), allocatable :: stf_spectrum_double
  double precision, dimension(:), allocatable :: stf_time_diff
  integer :: jt,jf,IOASKI

  ! even if rank 0 does not hold any ASKI output, let it in, since it writes the logs
  if(ASKI_np_local <= 0 .and. myrank/=0) return

  ! ################################################################################################
  ! in the future: make this routine independent of the knowledge about the chosen stf 
  ! by reading in output OUTPUT_FILES/plot_source_time_function.txt (assuming it is produced)
  ! ################################################################################################

  allocate(stf_time_diff(NSTEP))

  ! FIRST COMPUTE COMPLETE SOURCE TIME FUNCTION IN TIME DOMAIN
  do jt = 1,NSTEP
     !comp_source_time_function(dble(it-1)*DT-t0-tshift_cmt(isource),hdur_gaussian(isource))
     stf_time_diff(jt) = comp_source_time_function(dble(jt-1)*DT-t0-tshift_cmt(1),hdur_gaussian(1))
  end do
  ! THEN DIFFERENTIATE SOURCE TIME FUNCTION BY SIMPLE FINITE DIFFERENCES
  do jt = 1,NSTEP-1
     stf_time_diff(jt) = (stf_time_diff(jt+1)-stf_time_diff(jt))/DT
  end do
  ! make the derivative continuous at the end
  ! however, when using this properly, it should be zero anyway AND ADDITIONALLY it should be tapered to zero
  stf_time_diff(NSTEP) = stf_time_diff(NSTEP-1)

  ! NOW TRANSFORM STF TO SPECTRUM AT THE CHOSEN DISCRETE FREQUENCIES

  allocate(stf_spectrum_double(ASKI_nf))
!!$  stf_spectrum_double = (0.d0,0.d0)
!!$  do jt = 1,NSTEP
!!$     do jf = 1,ASKI_nf
!!$        stf_spectrum_double(jf) = stf_spectrum_double(jf) + stf_time_diff(jt)*ASKI_efactors_tapered(jf,jt)
!!$     end do ! jf
!!$  end do ! jt
  stf_spectrum_double = matmul(ASKI_efactors_tapered,stf_time_diff)

  ! RANK 0 WRITES OUT LOGS CONTAINING THE DIFFERENTIATED STF AND THE SPECTRUM WHICH IS DECONVOLVED
  if(myrank==0) then
     write(IMAIN,*) "WILL DECONVOLVE DERIVATIVE OF QUASI-HEAVISIDE SOURCE-TIME-FUNCTION FROM VELOCITY FIELD ",&
          "(LOGS CONTAINING THIS DIFFERENTIATED TIME-SERIES AND ITS SPECTRUM ARE WRITTEN NOW TO OUTPUT_FILES/)"
     call write_ASKI_log("LOG_ASKI_DECONVOLVE_stf.txt","As indicated in Par_file_ASKI by key 'ASKI_DECONVOLVE_STF', "//&
          "the derivative of the source-time-function is deconvolved from the velocity field. "//&
          "However, for now a quasi-Heaviside stf is ASSUMED for this purpose, any other actually chosen stf is "//&
          "NOT accounted for! File 'LOG_ASKI_DECONVOLVE_stf_diff.dat' contains the derivative of the quasi-Heaviside "//&
          "and file 'LOG_ASKI_DECONVOLVE_stf_diff_spectrum.dat' contains its spectrum at the ASKI frequencies, "//&
          "which is deconvolved from the kernel spectra (by simple division of the complex numbers)")

     call get_file_unit_ASKI(IOASKI)
     open(unit=IOASKI,file='OUTPUT_FILES/'//'LOG_ASKI_DECONVOLVE_stf_diff.dat',&
          form='formatted',status='unknown',action='write')
     do jt = 1,NSTEP
        write(IOASKI,*) real(dble(jt-11)*DT),real(stf_time_diff(jt))
     end do ! jt
     close(IOASKI)
     call get_file_unit_ASKI(IOASKI)
     open(unit=IOASKI,file='OUTPUT_FILES/'//'LOG_ASKI_DECONVOLVE_stf_diff_spectrum.dat',&
          form='formatted',status='unknown',action='write')
     do jf = 1,ASKI_nf
        write(IOASKI,*) (jf-1)*ASKI_df, stf_spectrum_double(jf), abs(stf_spectrum_double(jf))
     end do ! jf
     close(IOASKI)

     ! IF myrank==0 IS ONLY HERE TO WRITE THE DECONVOLVE LOG OUTPUT, BUT DOES NOT HAVE ANY ASKI OUTPUT, 
     ! THEN DEALLOCATE THE DECONVOLVE STUFF AND LEAVE THIS ROUTINE
     if(ASKI_np_local <= 0) then
        deallocate(stf_time_diff,stf_spectrum_double)
        return
     end if
  end if

  if(ASKI_DFT_double) then
     do jf = 1,ASKI_nf
        ASKI_spectra_local_double(:,jf,:) = ASKI_spectra_local_double(:,jf,:) / stf_spectrum_double(jf)
     end do ! jf
  else ! ASKI_DFT_double
     do jf = 1,ASKI_nf
        ASKI_spectra_local_single(:,jf,:) = ASKI_spectra_local_single(:,jf,:) / stf_spectrum_double(jf)
     end do ! jf
  end if ! ASKI_DFT_double

  deallocate(stf_time_diff,stf_spectrum_double)
end subroutine deconvolve_diff_stf_from_ASKI_output
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_output_files()

  use specfem_par,only: myrank,NPROCTOT

  implicit none

  integer :: jf,ip,iproc,IOASKI,specfem_version
  complex, dimension(:,:), allocatable :: spectrum_one_frequency
  character(len=509) :: filename
  logical :: file_exists
  character (len=7) :: open_status
  integer, parameter :: itag = 100

  if(myrank == 0) then

     ! THE OUTPUT SPECTRA ARE WRITTEN IN SINGLE PRECISION ONLY, REGARDLES IF ASKI_DFT_double OR NOT
     allocate(spectrum_one_frequency(sum(ASKI_np_local_all),9))

     do jf = 1,ASKI_nf

        ip = 0

        ! this is me, rank 0
        if(ASKI_np_local > 0) then
           if(ASKI_DFT_double) then
              spectrum_one_frequency(1:ASKI_np_local,:) = transpose(cmplx(ASKI_spectra_local_double(:,jf,:)))
           else
              spectrum_one_frequency(1:ASKI_np_local,:) = transpose(ASKI_spectra_local_single(:,jf,:))
           end if
           ip = ip + ASKI_np_local
        end if
        do iproc = 1,NPROCTOT-1
           if(ASKI_np_local_all(iproc+1) > 0) then
              call recv_c(spectrum_one_frequency(ip+1:ip+ASKI_np_local_all(iproc+1),:), ASKI_np_local_all(iproc+1)*9, iproc,itag)
              ip = ip + ASKI_np_local_all(iproc+1)
           end if
        end do ! iproc

        ! write spectrum_one_frequency to file file basename.jf###### , where '######' contains the frequency index
        write(filename,"(a,i6.6)") trim(ASKI_outfile)//".jf",ASKI_jf(jf)
        ! do not worry about whether to overwrite or not, if program comes here, it already has been checked
        ! if the value of OVERWRITE_ASKI_OUTPUT is in conflict with existing files
        inquire(file=filename,exist=file_exists)
        if(file_exists) then
           open_status = 'replace'
        else
           open_status = 'new'
        end if

        ! specfem_version
        ! 1 = SPECFEM3D_GLOBE, 2 = SPECFEM3D_Cartesian
        specfem_version = 1

        ! also write NPROCTOT to file, as a safety feature in order to assure that wavefield points (written only in ASKI_main_file
        ! and kernel values  (written only in kernel files) are not confused
        call get_file_unit_ASKI(IOASKI)
        open(unit=IOASKI, file=filename, status=open_status, access="stream", form="unformatted", action = 'write')
        write(IOASKI) specfem_version,length_ASKI_output_ID,ASKI_output_ID,NPROCTOT,sum(ASKI_np_local_all),&
             ASKI_df,ASKI_jf(jf)
        write(IOASKI) spectrum_one_frequency
        close(IOASKI)

     end do ! jf

     deallocate(spectrum_one_frequency)

     ! write LOG_ASKI_finish.txt
     call write_ASKI_log('LOG_ASKI_finish.txt',"successfully created ASKI output, as specified in 'LOG_ASKI_start.txt'")

  else ! myrank == 0
     
     if(ASKI_np_local > 0) then
        do jf = 1,ASKI_nf
           if(ASKI_DFT_double) then
              ! convert to single precision before sending
              call send_c(transpose(cmplx(ASKI_spectra_local_double(:,jf,:))), ASKI_np_local*9, 0, itag)
           else
              call send_c(transpose(ASKI_spectra_local_single(:,jf,:)), ASKI_np_local*9, 0, itag)
           end if
        end do
     end if

  end if ! myrank == 0

end subroutine write_ASKI_output_files
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_log(filename,log_message)
  implicit none
  character(len=*) :: filename,log_message
  integer :: fu
  call get_file_unit_ASKI(fu)
  open(unit=fu,file='OUTPUT_FILES/'//trim(filename),&
       form='formatted',status='unknown',action='write')
  write(fu,*) trim(log_message)
  close(fu)
end subroutine write_ASKI_log
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_log_start()
  use specfem_par,only: NPROCTOT,NGLLX,NGLLY,NGLLZ
  integer :: fu,i
  character(len=509) :: filename
  character(len=100) :: overwrite_message
  logical :: file_exists
  real :: size_main_file,size_jf_file
  call get_file_unit_ASKI(fu)
  open(unit=fu,file='OUTPUT_FILES/'//'LOG_ASKI_start.txt',&
       form='formatted',status='unknown',action='write')

  if(ASKI_MAIN_FILE_ONLY) then
     write(fu,*) "ONLY THE MAIN ASKI OUTPUT FILE WILL BE PRODUCED, as indicated by the logical parameter "//&
          "'ASKI_MAIN_FILE_ONLY' in DATA/Par_file_ASKI"
     write(fu,*) "HENCE, NO FREQUENCY KERNEL OUTPUT FILES WILL BE WRITTEN, EVEN IF INDICATED BELOW IN THIS LOGFILE!"
     write(fu,*) "For reasons of debugging and checking, this output was kept nevertheless."
     write(fu,*) ""
  end if

  write(fu,*) "Hello, this is SPECFEM3D_GLOBE for ASKI"
  write(fu,*) ""
  write(fu,*) "computing ASKI output now on ",NPROCTOT," procs with following parameters:"
  write(fu,*) ""
  if(COMPUTE_ASKI_GREEN_FUNCTION) then
     write(fu,*) "computing Green function component '",trim(ASKI_GREEN_FUNCTION_COMPONENT),"'"
     write(fu,*) "   ATTENTION: this means that a single force source is placed at the source position in this direction,"
     write(fu,*) "              ACTING WITH 1.e12 Newton FOR STABILITY REASONS (1 N is just too small). "
     write(fu,*) "              This additional factor of 1.e12 is compensated later on in the kernel computation!"
     write(fu,*) "              Note that the STANDARD SEISMOGRAM OUTPUT is displacement in meter w.r.t. 1.e12 N force!"
     write(fu,*) "              (or equivalently displacement in pico meters w.r.t. 1 N force)."
  else
     write(fu,*) "this is NO Green function, i.e. computing wavefield from seismic source as indicated by file ",&
          "DATA/CMTSOLUTION"
  end if
  write(fu,*) ""
  write(fu,*) "ASKI_type_inversion_grid = ",ASKI_type_inversion_grid
  write(fu,*) ""
  select case(ASKI_nchunk)
  case(1)
     write(fu,*) "ASKI output volume consists of ",ASKI_nchunk," chunk"
     write(fu,*) "   chunk centered in"
     write(fu,*) "     ASKI_clat = ",ASKI_clat
     write(fu,*) "     ASKI_clon = ",ASKI_clon
     write(fu,*) "   having width"
     write(fu,*) "     ASKI_wlat = ",ASKI_wlat
     write(fu,*) "     ASKI_wlon = ",ASKI_wlon
     write(fu,*) "   and azimuth rotation angle"
     write(fu,*) "     ASKI_rot_gamma = ",ASKI_rot_gamma
  case(2,3)
     write(fu,*) "ASKI output volume consists of ",ASKI_nchunk," chunks"
     write(fu,*) "   first chunk centered in"
     write(fu,*) "     ASKI_clat = ",ASKI_clat
     write(fu,*) "     ASKI_clon = ",ASKI_clon
     write(fu,*) "   having azimuth rotation angle"
     write(fu,*) "     ASKI_rot_gamma = ",ASKI_rot_gamma
  end select
  write(fu,*) "   depth range"
  write(fu,*) "     ASKI_rmin = ",ASKI_rmin
  write(fu,*) "     ASKI_rmax = ",ASKI_rmax
  select case(ASKI_type_inversion_grid)
  case(4)
     write(fu,*) "   NGLLX = ",NGLLX
     write(fu,*) "   NGLLY = ",NGLLY
     write(fu,*) "   NGLLZ = ",NGLLZ
     write(fu,*) "   total number of inversion grid cells = ",sum(ASKI_np_local_all)/(NGLLX*NGLLY*NGLLZ)
  end select
  write(fu,*) ""
  write(fu,*) "local number of wavefield points (at which ASKI output is computed):"
  do i = 1,NPROCTOT
     write(fu,*) "   proc ",i-1," : ",ASKI_np_local_all(i)
  end do ! iproc
  write(fu,*) "in total : ",sum(ASKI_np_local_all)
  write(fu,*) ""
  write(fu,*) "output spectra will be computed for ",ASKI_nf," frequencies f = df*jf defined by "
  write(fu,*) "   df = ",ASKI_df
  write(fu,*) "   jf = ",ASKI_jf
  write(fu,*) ""
  size_main_file = (length_ASKI_output_ID + 4.*(7+ASKI_nf+6*sum(ASKI_np_local_all)))/1048576.
  size_jf_file = (length_ASKI_output_ID + 4.*(6+2*9*sum(ASKI_np_local_all)))/1048576.
  write(fu,*) "in total ",ASKI_nf*size_jf_file+size_main_file,&
       " MiB output will be written to ",ASKI_nf+1," files with base_filename = "
  write(fu,*) "'"//trim(ASKI_outfile)//"' :"

  overwrite_message = ''
  if(OVERWRITE_ASKI_OUTPUT) then
     inquire(file=trim(ASKI_outfile)//".main",exist=file_exists)
     if(file_exists) then
        overwrite_message = 'exists and will be overwritten'
     else
        overwrite_message = 'does not exist and will be newly created'
     end if
  else
     ! if the file existed, this would have been detected above already, so if program comes here, the file does not exist
     overwrite_message = 'does not exist and will be newly created'
  end if
  write(fu,*) "   base_filename.main  (",size_main_file," MiB)  "//trim(overwrite_message)

  do i = 1,ASKI_nf
     write(filename,"(a,i6.6)") trim(ASKI_outfile)//".jf",ASKI_jf(i)
     overwrite_message = ''
     if(OVERWRITE_ASKI_OUTPUT) then
        inquire(file=trim(filename),exist=file_exists)
        if(file_exists) then
           overwrite_message = 'exists and will be overwritten'
        else
           overwrite_message = 'does not exist and will be newly created'
        end if
     else
        ! if the file existed, this would have been detected above already, so if program comes here, the file does not exist
        overwrite_message = 'does not exist and will be newly created'
     end if
     write(filename,"('   base_filename.jf',i6.6,'  ')") ASKI_jf(i)
     write(fu,*) trim(filename),"  (",size_jf_file," MiB)  "//trim(overwrite_message)
  end do ! i
  close(fu)
end subroutine write_ASKI_log_start
!
!-------------------------------------------------------------------------------------------
!
  subroutine get_file_unit_ASKI(unit_out)
   implicit none
   integer :: unit_out,fu
   logical :: is_open
   integer, parameter :: min_unit = 20
   integer, parameter :: max_unit = 99

   unit_out = -1
   do fu = min_unit, max_unit
      inquire(unit = fu, opened = is_open)
      if (.not. is_open) then
         unit_out = fu
         return
      end if
   end do
   call exit_MPI_without_rank('no file unit between 20 and 99 available')
 end subroutine get_file_unit_ASKI


  end module specfem_for_ASKI_par
