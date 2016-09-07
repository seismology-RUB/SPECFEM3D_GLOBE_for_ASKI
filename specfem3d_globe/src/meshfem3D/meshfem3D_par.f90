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

!
!-------------------------------------------------------------------------------------------------
!

  module meshfem3D_models_par

!---
!
! ADD YOUR MODEL HERE
!
!---

  use constants

  use shared_input_parameters,only: &
    ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS, &
    ATTENUATION,USE_FULL_TISO_MANTLE

  use shared_compute_parameters, only: &
    REFERENCE_1D_MODEL,THREE_D_MODEL, &
    HONOR_1D_SPHERICAL_MOHO,CRUSTAL,ONE_CRUST,CASE_3D,TRANSVERSE_ISOTROPY, &
    ISOTROPIC_3D_MANTLE,ANISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
    ATTENUATION_3D, &
    ANISOTROPIC_INNER_CORE, &
    CEM_REQUEST, CEM_ACCEPT

  implicit none

! model_attenuation_variables
  type model_attenuation_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinuities Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    integer                                   :: Qn                 ! Number of points
    integer dummy_pad ! padding 4 bytes to align the structure
  end type model_attenuation_variables
  type (model_attenuation_variables) AM_V
! model_attenuation_variables

! model_attenuation_storage_var
  type model_attenuation_storage_var
    sequence
    double precision, dimension(:,:), pointer :: tau_e_storage
    double precision, dimension(:), pointer :: Qmu_storage
    integer Q_resolution
    integer Q_max
  end type model_attenuation_storage_var
  type (model_attenuation_storage_var) AM_S
! model_attenuation_storage_var

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables
  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables

! GLL model_variables
  type model_gll_variables
    sequence
    ! tomographic iteration model on GLL points
    double precision :: scale_velocity,scale_density
    ! isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vs_new,vp_new,rho_new
    ! transverse isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vsv_new,vpv_new, &
      vsh_new,vph_new,eta_new
    logical :: MODEL_GLL
    logical,dimension(3) :: dummy_pad ! padding 3 bytes to align the structure
  end type model_gll_variables
  type (model_gll_variables) MGLL_V

! bathymetry and topography: use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

! for ellipticity
  double precision,dimension(NR) :: rspl,espl,espl2
  integer :: nspl

  end module meshfem3D_models_par


!
!-------------------------------------------------------------------------------------------------
!


  module meshfem3D_par

! main parameter module for specfem simulations

  use constants

  use shared_parameters

  implicit none

  ! correct number of spectral elements in each block depending on chunk type
  integer :: npointot

  ! proc number for MPI process
  integer :: myrank

  ! check area and volume of the final mesh
  double precision :: volume_total

  ! check Earth mass and Earth center of mass computed in the final mesh
  double precision :: Earth_mass_total
  double precision :: Earth_center_of_mass_x_total,Earth_center_of_mass_y_total,Earth_center_of_mass_z_total

  ! arrays containing the positions of the observation points in non-dimensionalized value for gravity integrals
  ! the 1D equivalenced versions are for the FORCE_VECTORIZATION version of the loops
  double precision, dimension(NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX) :: x_observation,y_observation,z_observation
  double precision, dimension(NTOTAL_OBSERVATION) :: x_observation1D,y_observation1D,z_observation1D
  equivalence(x_observation,x_observation1D)
  equivalence(y_observation,y_observation1D)
  equivalence(z_observation,z_observation1D)

  double precision, dimension(NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX) :: lon_observation,lat_observation

  ! arrays containing the computed fields for gravity integrals at the observation points
  ! the 1D equivalenced versions are for the FORCE_VECTORIZATION version of the loops
  double precision, dimension(NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX) :: g_x,g_y,g_z,G_xx,G_yy,G_zz,G_xy,G_xz,G_yz, &
                                                                            temporary_array_for_sum
  double precision, dimension(NTOTAL_OBSERVATION) :: g_x1D,g_y1D,g_z1D,G_xx1D,G_yy1D,G_zz1D,G_xy1D,G_xz1D,G_yz1D
  equivalence(g_x,g_x1D)
  equivalence(g_y,g_y1D)
  equivalence(g_z,g_z1D)
  equivalence(G_xx,G_xx1D)
  equivalence(G_yy,G_yy1D)
  equivalence(G_zz,G_zz1D)
  equivalence(G_xy,G_xy1D)
  equivalence(G_xz,G_xz1D)
  equivalence(G_yz,G_yz1D)

  ! for loop on all the slices
  integer :: iregion_code
  integer :: iproc_xi,iproc_eta,ichunk

  ! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix
  double precision :: ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

  ! for some statistics for the mesh
  integer :: numelem_crust_mantle,numelem_outer_core,numelem_inner_core
  integer :: numelem_total

  ! timer MPI
  double precision :: time_start

  ! addressing for all the slices
  integer, dimension(:), allocatable :: ichunk_slice,iproc_xi_slice,iproc_eta_slice
  integer, dimension(:,:,:), allocatable :: addressing

  ! memory size of all the static arrays
  double precision :: static_memory_size

  integer :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
             NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
             NSPEC_INNER_CORE_ATTENUATION, &
             NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
             NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
             NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
             NSPEC_CRUST_MANTLE_ADJOINT, &
             NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
             NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
             NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
             NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
             NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION

  ! this for the different corners of the slice (which are different if the superbrick is cut)
  ! 1 : xi_min, eta_min
  ! 2 : xi_max, eta_min
  ! 3 : xi_max, eta_max
  ! 4 : xi_min, eta_max
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER

  ! 1 -> min, 2 -> max
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE

  ! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore
  ! parameters needed to store the radii of the grid points
  ! in the spherically symmetric Earth
  integer, dimension(:), allocatable :: idoubling
  integer, dimension(:,:,:,:), allocatable :: ibool

  ! this for non blocking MPI
  logical, dimension(:), allocatable :: is_on_a_slice_edge

  end module meshfem3D_par

!
!-------------------------------------------------------------------------------------------------
!

  module create_regions_mesh_par

  use constants,only: NGLLX,NGLLY,NGLLZ,NGNOD,NGNOD2D,NDIM,NDIM2D

  implicit none

  ! topology of the elements
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

  ! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

  ! 3D shape functions and their derivatives
  double precision, dimension(NGNOD,NGLLX,NGLLY,NGLLZ) :: shape3D
  double precision, dimension(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ) :: dershape3D

  ! 2D shape functions and their derivatives
  double precision, dimension(NGNOD2D,NGLLY,NGLLZ) :: shape2D_x
  double precision, dimension(NGNOD2D,NGLLX,NGLLZ) :: shape2D_y
  double precision, dimension(NGNOD2D,NGLLX,NGLLY) :: shape2D_bottom,shape2D_top
  double precision, dimension(NDIM2D,NGNOD2D,NGLLY,NGLLZ) :: dershape2D_x
  double precision, dimension(NDIM2D,NGNOD2D,NGLLX,NGLLZ) :: dershape2D_y
  double precision, dimension(NDIM2D,NGNOD2D,NGLLX,NGLLY) :: dershape2D_bottom,dershape2D_top

  end module create_regions_mesh_par

!
!-------------------------------------------------------------------------------------------------
!

  module create_regions_mesh_par2

  use constants,only: CUSTOM_REAL,N_SLS,MAX_STRING_LEN

  implicit none

  integer :: nspec_stacey,nspec_actually,nspec_att

  integer :: ifirst_region,ilast_region
  integer, dimension(:), allocatable :: perm_layer

  ! for model density and anisotropy
  integer :: nspec_ani
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,dvpstore, &
    kappavstore,kappahstore,muvstore,muhstore,eta_anisostore

  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

  ! boundary locator
  logical, dimension(:,:), allocatable :: iboun

  ! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  ! mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx,rmassy,rmassz
  integer :: nglob_xy

  ! mass matrices for backward simulation when ROTATION is .true.
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_rmassx,b_rmassy

  ! mass matrix and bathymetry for ocean load
  integer :: nglob_oceans
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

  ! number of elements on the boundaries
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  ! boundary parameters locator
  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax, &
    ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

  ! 2-D Jacobians and normals
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    jacobian2D_xmin,jacobian2D_xmax, &
    jacobian2D_ymin,jacobian2D_ymax,jacobian2D_bottom,jacobian2D_top

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top

  ! MPI cut-planes parameters along xi and along eta
  logical, dimension(:,:), allocatable :: iMPIcut_xi,iMPIcut_eta

  ! Stacey, indices for Clayton-Engquist absorbing conditions
  integer, dimension(:,:), allocatable :: nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

  ! attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),   allocatable :: Qmu_store
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: tau_e_store
  double precision, dimension(N_SLS) :: tau_s
  double precision :: T_c_source

  ! element layers
  integer :: NUMBER_OF_MESH_LAYERS,layer_shift, &
    first_layer_aniso,last_layer_aniso
  logical :: USE_ONE_LAYER_SB

  ! layer stretching
  double precision, dimension(:,:), allocatable :: stretch_tab

  ! Boundary Mesh
  integer :: NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho
  integer, dimension(:), allocatable :: ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
    ibelm_670_top,ibelm_670_bot
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_moho,normal_400,normal_670
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_moho,jacobian2D_400,jacobian2D_670

  integer :: ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot, &
    ispec2D_670_top,ispec2D_670_bot
  double precision :: r_moho,r_400,r_670

  ! flags for transverse isotropic elements
  logical, dimension(:), allocatable :: ispec_is_tiso

  ! name of the database file
  character(len=MAX_STRING_LEN) :: prname, prname_adios

  end module create_regions_mesh_par2

!
!-------------------------------------------------------------------------------------------------
!

  module create_MPI_interfaces_par

  use constants,only: &
    CUSTOM_REAL,NDIM,IMAIN, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE, &
    NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR

  implicit none

  ! indirect addressing for each message for faces and corners of the chunks
  ! a given slice can belong to at most one corner and at most two faces

  ! number of faces between chunks
  integer :: NUMMSGS_FACES

  ! number of corners between chunks
  integer :: NCORNERSCHUNKS

  ! number of message types
  integer :: NUM_MSG_TYPES

  !-----------------------------------------------------------------
  ! assembly
  !-----------------------------------------------------------------

  ! ---- arrays to assemble between chunks
  ! communication pattern for faces between chunks
  integer, dimension(:),allocatable :: iprocfrom_faces,iprocto_faces,imsg_type
  ! communication pattern for corners between chunks
  integer, dimension(:),allocatable :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

  ! indirect addressing for each corner of the chunks
  integer, dimension(:,:),allocatable :: iboolcorner

  ! chunk faces
  integer, dimension(:,:),allocatable :: iboolfaces
  integer, dimension(NUMFACES_SHARED) :: npoin2D_faces
  integer :: NGLOB2DMAX_XY

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(:),allocatable :: iboolleft_xi,iboolright_xi
  integer, dimension(:),allocatable :: iboolleft_eta,iboolright_eta
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_all,npoin2D_eta_all
  integer :: npoin2D_xi,npoin2D_eta

  ! 1-D addressing
  integer :: NGLOB1D_RADIAL_MAX
  integer,dimension(:),allocatable :: ibool1D_leftxi_lefteta,ibool1D_rightxi_lefteta, &
                                    ibool1D_leftxi_righteta,ibool1D_rightxi_righteta

  double precision,dimension(:,:),allocatable :: xyz1D_leftxi_lefteta,xyz1D_rightxi_lefteta, &
                                    xyz1D_leftxi_righteta,xyz1D_rightxi_righteta

  ! this for non blocking MPI

  ! buffers for send and receive between faces of the slices and the chunks
  ! we use the same buffers to assemble scalars and vectors because vectors are
  ! always three times bigger and therefore scalars can use the first part
  ! of the vector buffer in memory even if it has an additional index here
  integer :: npoin2D_max_all_CM_IC

  ! buffers for send and receive between corners of the chunks
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: &
    buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar

  ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: &
     buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector

  end module create_MPI_interfaces_par

!
!-------------------------------------------------------------------------------------------------
!

  module MPI_crust_mantle_par

  use constants,only: CUSTOM_REAL,NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR

  implicit none

  ! collected MPI interfaces
  !--------------------------------------
  ! MPI crust/mantle mesh
  !--------------------------------------
  integer :: num_interfaces_crust_mantle
  integer :: max_nibool_interfaces_cm
  integer, dimension(:), allocatable :: my_neighbours_crust_mantle,nibool_interfaces_crust_mantle
  integer, dimension(:,:), allocatable :: ibool_interfaces_crust_mantle

  !--------------------------------------
  ! crust mantle
  !--------------------------------------
  integer :: NSPEC_CRUST_MANTLE
  integer :: NGLOB_CRUST_MANTLE

  integer :: NGLOB1D_RADIAL_CM
  integer :: NGLOB2DMAX_XMIN_XMAX_CM
  integer :: NGLOB2DMAX_YMIN_YMAX_CM
  integer :: NSPEC2DMAX_XMIN_XMAX_CM
  integer :: NSPEC2DMAX_YMIN_YMAX_CM
  integer :: NSPEC2D_BOTTOM_CM
  integer :: NSPEC2D_TOP_CM

  real(kind=CUSTOM_REAL), dimension(:),allocatable :: &
    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  ! assembly
  integer, dimension(NUMFACES_SHARED) :: npoin2D_faces_crust_mantle
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle

  ! indirect addressing for each corner of the chunks
  integer, dimension(:,:),allocatable :: iboolcorner_crust_mantle

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(:),allocatable :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(:),allocatable :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle

  integer, dimension(:,:),allocatable :: iboolfaces_crust_mantle

  ! inner / outer elements crust/mantle region
  integer :: num_phase_ispec_crust_mantle
  integer :: nspec_inner_crust_mantle,nspec_outer_crust_mantle
  integer, dimension(:,:), allocatable :: phase_ispec_inner_crust_mantle

  ! mesh coloring
  integer :: num_colors_outer_crust_mantle,num_colors_inner_crust_mantle
  integer,dimension(:),allocatable :: num_elem_colors_crust_mantle

  end module MPI_crust_mantle_par

!
!-------------------------------------------------------------------------------------------------
!

  module MPI_inner_core_par

  use constants,only: CUSTOM_REAL,NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR

  implicit none

  !--------------------------------------
  ! MPI inner core mesh
  !--------------------------------------
  integer :: num_interfaces_inner_core
  integer :: max_nibool_interfaces_ic
  integer, dimension(:), allocatable :: my_neighbours_inner_core,nibool_interfaces_inner_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_inner_core

  !--------------------------------------
  ! inner core
  !--------------------------------------
  integer :: NSPEC_INNER_CORE
  integer :: NGLOB_INNER_CORE

  integer :: NGLOB1D_RADIAL_IC
  integer :: NGLOB2DMAX_XMIN_XMAX_IC
  integer :: NGLOB2DMAX_YMIN_YMAX_IC
  integer :: NSPEC2DMAX_XMIN_XMAX_IC
  integer :: NSPEC2DMAX_YMIN_YMAX_IC
  integer :: NSPEC2D_BOTTOM_IC
  integer :: NSPEC2D_TOP_IC

  real(kind=CUSTOM_REAL), dimension(:),allocatable :: &
    xstore_inner_core,ystore_inner_core,zstore_inner_core

  ! for matching with central cube in inner core
  integer, dimension(:), allocatable :: sender_from_slices_to_cube
  integer, dimension(:,:), allocatable :: ibool_central_cube
  double precision, dimension(:,:), allocatable :: buffer_slices,buffer_slices2
  double precision, dimension(:,:,:), allocatable :: buffer_all_cube_from_slices
  integer :: nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube, &
    npoin2D_cube_from_slices,receiver_cube_from_slices

  ! bottom inner core / top central cube
  integer, dimension(:),allocatable :: ibelm_bottom_inner_core

  integer, dimension(NUMFACES_SHARED) :: npoin2D_faces_inner_core
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_inner_core,npoin2D_eta_inner_core

  ! indirect addressing for each corner of the chunks
  integer, dimension(:,:),allocatable :: iboolcorner_inner_core

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(:),allocatable :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(:),allocatable :: iboolleft_eta_inner_core,iboolright_eta_inner_core

  integer, dimension(:,:),allocatable :: iboolfaces_inner_core

  ! inner / outer elements inner core region
  integer :: num_phase_ispec_inner_core
  integer :: nspec_inner_inner_core,nspec_outer_inner_core
  integer, dimension(:,:), allocatable :: phase_ispec_inner_inner_core

  ! mesh coloring
  integer :: num_colors_outer_inner_core,num_colors_inner_inner_core
  integer,dimension(:),allocatable :: num_elem_colors_inner_core

  end module MPI_inner_core_par

!
!-------------------------------------------------------------------------------------------------
!

  module MPI_outer_core_par

  use constants,only: CUSTOM_REAL,NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR

  implicit none

  !--------------------------------------
  ! MPI outer core mesh
  !--------------------------------------
  integer :: num_interfaces_outer_core
  integer :: max_nibool_interfaces_oc
  integer, dimension(:), allocatable :: my_neighbours_outer_core,nibool_interfaces_outer_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_outer_core

  !--------------------------------------
  ! outer core
  !--------------------------------------
  integer :: NSPEC_OUTER_CORE
  integer :: NGLOB_OUTER_CORE

  integer :: NGLOB1D_RADIAL_OC
  integer :: NGLOB2DMAX_XMIN_XMAX_OC
  integer :: NGLOB2DMAX_YMIN_YMAX_OC
  integer :: NSPEC2DMAX_XMIN_XMAX_OC
  integer :: NSPEC2DMAX_YMIN_YMAX_OC
  integer :: NSPEC2D_BOTTOM_OC
  integer :: NSPEC2D_TOP_OC

  real(kind=CUSTOM_REAL), dimension(:),allocatable :: &
    xstore_outer_core,ystore_outer_core,zstore_outer_core

  ! assembly
  integer, dimension(NUMFACES_SHARED) :: npoin2D_faces_outer_core
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_outer_core,npoin2D_eta_outer_core

  ! indirect addressing for each corner of the chunks
  integer, dimension(:,:),allocatable :: iboolcorner_outer_core

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(:),allocatable :: iboolleft_xi_outer_core,iboolright_xi_outer_core
  integer, dimension(:),allocatable :: iboolleft_eta_outer_core,iboolright_eta_outer_core

  integer, dimension(:,:),allocatable :: iboolfaces_outer_core

  ! inner / outer elements outer core region
  integer :: num_phase_ispec_outer_core
  integer :: nspec_inner_outer_core,nspec_outer_outer_core
  integer, dimension(:,:), allocatable :: phase_ispec_inner_outer_core

  ! mesh coloring
  integer :: num_colors_outer_outer_core,num_colors_inner_outer_core
  integer,dimension(:),allocatable :: num_elem_colors_outer_core

  end module MPI_outer_core_par
!----------------------------------------------------------------------------
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
!
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.2 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!
!
  module ASKI_external_model
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! THIS MODULE ORIGINATES FROM  FILE model_external_values.f90 
! FROM 
!    SPECFEM3D_Cartesian version 3.0 for ASKI 1.0
! WITH ITS MAIN STRUCTURE BEING BASED ON THE ORIGINAL FILE model_external_values.f90
! OF THE SPECFEM3D_Cartesian SOFTWARE PACKAGE BY NOVEMBER 2010.
! In addition, the checker-board model "model_ASKI_checker" was added, here. 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

    use constants

    implicit none

    ! defined by flag "USE_ASKI_BACKGROUND_MODEL" in Par_file_ASKI, logical "use_ASKI_background_model" indicates
    ! whether a global 1D background model should be set.
    ! (this will replace the SPECFEM3D_GLOBE model everywhere, before possibly imposing the inverted model)
    ! If IMPOSE_ASKI_BACKGROUND_MODEL = .true. , file_ASKI_background_model then contains the value of
    ! FILE_ASKI_BACKGROUND_MODEL, i.e. the model file, relative to DATA/
    logical :: use_ASKI_background_model = .false.
    character(len=MAX_STRING_LEN) :: file_ASKI_background_model = ''

    ! defined by flag "IMPOSE_ASKI_INVERTED_MODEL" in Par_file_ASKI, logical "impose_ASKI_inverted_model" indicates
    ! whether the ASKI external model should be imposed.
    ! If IMPOSE_ASKI_INVERTED_MODEL = .true. , file_ASKI_inverted_model then contains the value of
    ! FILE_ASKI_INVERTED_MODEL, i.e. the model file, relative to DATA/
    logical :: impose_ASKI_inverted_model = .false.
    character(len=MAX_STRING_LEN) :: file_ASKI_inverted_model = ''

    ! defined by flag "IMPOSE_ASKI_CHECKER_MODEL" in Par_file_ASKI, logical "impose_ASKI_checker_model" indicates
    ! whether SPECFEM3D_GLOBE should impose relative checkerboard anomalies onto the final model.
    ! If IMPOSE_ASKI_CHECKER_MODEL = .true. , file_ASKI_checker_model then contains the value of
    ! FILE_ASKI_CHECKER_MODEL, i.e. the checker model file, relative to DATA/
    logical :: impose_ASKI_checker_model = .false.
    character(len=MAX_STRING_LEN) :: file_ASKI_checker_model = ''


    integer :: model_ASKI_myrank


    type model_ASKI_checker
       ! checker board definition of relative model anomalies

       ! bchecker_* have size (nchecker_*,2) and contains for each checker interval the lower and upper boundary

       integer :: nchecker_radius ! number of depth layers with checkers
       double precision, dimension(:,:), pointer :: bchecker_radius ! normalized (max radius = 1)

       integer :: nchecker_lat ! number of lateral checker cells in latitude direction of chunk
       double precision, dimension(:,:), pointer ::  bchecker_lat

       integer :: nchecker_lon ! number of lateral checker cells in latitude direction of chunk
       double precision, dimension(:,:), pointer ::  bchecker_lon

       ! model anomaly values
       integer :: nparam = 5 ! the parameterrs are: rho,vp,vs,Qmu,Qkappa (isotropic)
       real, dimension(5) :: anomaly ! 1 POSITIVE RELATIVE anomaly value (in %) for each parameter (will be used +-)
    end type model_ASKI_checker
    type (model_ASKI_checker) :: mAchk
    double precision :: ASKI_Mrot_chunk1(3,3)


    type model_ASKI_1Dbackground
       ! 1D spline-interpolated background model to replace the overall SPECFEM3D background model,
       ! before possibly imposing the ASKI inverted model 
       ! This mechanism was implemented in order to properly continue an inversion where the first
       ! iteration was done by a 1D method. Since the simulation domain is usually larger than the inversion
       ! domain, the inverted 3D model should be extendet to the rest of the inversion domain by the very 1D
       ! reference model that was used before by the 1D method.
       integer :: nlayers
       integer, dimension(:), pointer :: nnodes
       real, dimension(:,:), pointer :: depth,rho,vp,vs,Qmu,Qkappa  ! depth and parameter arrays
       real, dimension(:,:), pointer :: sprho,spvp,spvs,spQmu,spQkappa  ! spline parameters p = s''(xj)
    end type model_ASKI_1Dbackground
    type (model_ASKI_1Dbackground) :: mA1Db

    type model_ASKI_cells
       ! inversion grid cells
       integer :: ncell
       real, dimension(:,:), pointer :: cc
       real, dimension(:), pointer :: r
       integer :: max_nnb
       integer, dimension(:,:), pointer :: nb ! of size (max_nnb+1,ncell)
    end type model_ASKI_cells
    type (model_ASKI_cells) :: mAc

    type model_ASKI_isoLame
       ! model values for parametrization isoLame
       real :: maxr_rho,maxr_lambda,maxr_mu
       integer :: nval_rho,nval_lambda,nval_mu
       integer, dimension(:), pointer :: idx_rho,idx_lambda,idx_mu
       real, dimension(:), pointer :: rho,lambda,mu
    end type model_ASKI_isoLame
    type (model_ASKI_isoLame) :: mAisoL

    type model_ASKI_isoVelocity
       ! model values for parametrization isoVelocity
       real :: maxr_rho,maxr_vp,maxr_vs
       integer :: nval_rho,nval_vp,nval_vs
       integer, dimension(:), pointer :: idx_rho,idx_vp,idx_vs
       real, dimension(:), pointer :: rho,vp,vs
    end type model_ASKI_isoVelocity
    type (model_ASKI_isoVelocity) :: mAisoV

    ! model_ASKI_pmtrz has values according to integer parameters 
    ! ipmtrz_isoLame,ipmtrz_isoVelocity,... below (selects which variable is used: mAisoV,mAisoL,...)
    integer :: model_ASKI_pmtrz

    ! type interpolation between of external model that is used
    ! type = 1 : case "shepard_standard"
    !            interpolates model values by modified 3D Shepard interpolation with standard factor for influence radius
    ! type = 2 : case "shepard_factor_radius"
    !            interpolates model values by modified 3D Shepard interpolation with given factor for influence radius
    integer :: model_ASKI_interpolation_type

    ! additional parameters dependent on model_ASKI_interpolation_type
    real :: model_ASKI_factor_shepard_radius


    ! other definitions
    integer, parameter :: ipmtrz_isoLame = 1
    integer, parameter :: ipmtrz_isoVelocity = 2

  contains
!
!-------------------------------------------------------------------------------------------------
!
  subroutine broadcast_ASKI_external_model(myrank) !! ASKI FS FS: THIS ROUTINE IS CALLED FROM OUTSIDE THE MODULE IN setup_model()

! standard routine to setup model

  integer :: myrank,IOASKI

  ! local parameters
  character(len=800) :: error_message
  integer :: maxnnodes
  integer, dimension(1) :: i_array_one_value
  real, dimension(1) :: r_array_one_value
  logical, dimension(1) :: l_array_one_value

  ! dummy to ignore compiler warnings
  model_ASKI_myrank = myrank

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! the variables read are declared and stored in structure MEXT_V
  !if(myrank == 0) call read_external_model()

  ! broadcast the information read on the master to the nodes
  !call bcast_all_dp(MEXT_V%dvs, size(MEXT_V%dvs))

  if(myrank == 0) then

     ! read the values of variables use_ASKI_background_model, impose_ASKI_inverted_model
     call read_Par_file_ASKI()

     ! write logfile about the values of use_ASKI_background_model, impose_ASKI_inverted_model to OUTPUT_FILES
     call get_file_unit_model_ASKI(IOASKI)
     open(unit=IOASKI,file='OUTPUT_FILES/'//'LOG_ASKI_external_model.txt',&
          form='formatted',status='unknown',action='write')
     if(.not.use_ASKI_background_model) then
        write(IOASKI,*) "according to flag 'USE_ASKI_BACKGROUND_MODEL' in file '"//'DATA/'//&
          "Par_file_ASKI' (or the fact that the file is missing), this simulation will NOT use an external ASKI background model"
     else
        write(IOASKI,*) "according to flag 'USE_ASKI_BACKGROUND_MODEL' in file '"//'DATA/'//&
          "Par_file_ASKI', this simulation will use an external ASKI background model defined by file '"//&
          'DATA/'//trim(file_ASKI_background_model)//"'"
     end if
     if(.not.impose_ASKI_inverted_model) then
        write(IOASKI,*) "according to flag 'IMPOSE_ASKI_INVERTED_MODEL' in file '"//'DATA/'//&
          "Par_file_ASKI' (or the fact that the file is missing), this simulation will NOT impose an external ASKI inverted model"
     else
        write(IOASKI,*) "according to flag 'IMPOSE_ASKI_INVERTED_MODEL' in file '"//'DATA/'//&
          "Par_file_ASKI', this simulation will impose an external ASKI inverted model defined by file '"//&
          'DATA/'//trim(file_ASKI_inverted_model)//"'"
     end if
     if(.not.impose_ASKI_checker_model) then
        write(IOASKI,*) "according to flag 'IMPOSE_ASKI_CHECKER_MODEL' in file '"//'DATA/'//&
          "Par_file_ASKI' (or the fact that the file is missing), this simulation will NOT impose an external ASKI checker model"
     else
        write(IOASKI,*) "according to flag 'IMPOSE_ASKI_CHECKER_MODEL' in file '"//'DATA/'//&
          "Par_file_ASKI', this simulation will impose an external ASKI checker model defined by file '"//&
          'DATA/'//trim(file_ASKI_checker_model)//"'"
     end if
     close(IOASKI)

     if(use_ASKI_background_model) then
        call read_ASKI_external_background_model()
        maxnnodes = maxval(mA1Db%nnodes)
     end if
     if(impose_ASKI_inverted_model) call read_ASKI_external_inverted_model()
     if(impose_ASKI_checker_model)  call read_ASKI_external_checker_model()

  end if ! myrank == 0

  call synchronize_all()

  if(myrank == 0) l_array_one_value(1) = use_ASKI_background_model
  call bcast_all_l(l_array_one_value,1)
  if(myrank /= 0) use_ASKI_background_model = l_array_one_value(1)

  if(myrank == 0) l_array_one_value(1) = impose_ASKI_inverted_model
  call bcast_all_l(l_array_one_value,1)
  if(myrank /= 0) impose_ASKI_inverted_model = l_array_one_value(1)

  if(myrank == 0) l_array_one_value(1) = impose_ASKI_checker_model
  call bcast_all_l(l_array_one_value,1)
  if(myrank /= 0) impose_ASKI_checker_model = l_array_one_value(1)

  if(use_ASKI_background_model) then

     call bcast_all_ch(file_ASKI_background_model,len_trim(file_ASKI_background_model))

     if(myrank == 0) i_array_one_value(1) = mA1Db%nlayers
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mA1Db%nlayers = i_array_one_value(1)

     if(myrank == 0) i_array_one_value(1) = maxnnodes
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) maxnnodes = i_array_one_value(1)

     ! allocate for model values if I'm not rank 0
     if(myrank .ne. 0) allocate(mA1Db%nnodes(mA1Db%nlayers),mA1Db%depth(mA1Db%nlayers,maxnnodes), &
          mA1Db%rho(mA1Db%nlayers,maxnnodes),mA1Db%vp(mA1Db%nlayers,maxnnodes),mA1Db%vs(mA1Db%nlayers,maxnnodes), &
          mA1Db%Qmu(mA1Db%nlayers,maxnnodes),mA1Db%Qkappa(mA1Db%nlayers,maxnnodes),&
          mA1Db%sprho(mA1Db%nlayers,maxnnodes),mA1Db%spvp(mA1Db%nlayers,maxnnodes),mA1Db%spvs(mA1Db%nlayers,maxnnodes),&
          mA1Db%spQmu(mA1Db%nlayers,maxnnodes),mA1Db%spQkappa(mA1Db%nlayers,maxnnodes))
 
     call bcast_all_i(mA1Db%nnodes,size(mA1Db%nnodes))
     call bcast_all_r(mA1Db%depth,size(mA1Db%depth))
     call bcast_all_r(mA1Db%rho,size(mA1Db%rho))
     call bcast_all_r(mA1Db%vp,size(mA1Db%vp))
     call bcast_all_r(mA1Db%vs,size(mA1Db%vs))
     call bcast_all_r(mA1Db%Qmu,size(mA1Db%Qmu))
     call bcast_all_r(mA1Db%Qkappa,size(mA1Db%Qkappa))
     call bcast_all_r(mA1Db%sprho,size(mA1Db%sprho))
     call bcast_all_r(mA1Db%spvp,size(mA1Db%spvp))
     call bcast_all_r(mA1Db%spvs,size(mA1Db%spvs))
     call bcast_all_r(mA1Db%spQmu,size(mA1Db%spQmu))
     call bcast_all_r(mA1Db%spQkappa,size(mA1Db%spQkappa))
  end if ! use_ASKI_background_model

  if(impose_ASKI_inverted_model) then

     call bcast_all_ch(file_ASKI_inverted_model,len_trim(file_ASKI_inverted_model))

     if(myrank == 0) i_array_one_value(1) = model_ASKI_interpolation_type
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) model_ASKI_interpolation_type = i_array_one_value(1)

     select case(model_ASKI_interpolation_type)
     case( 2 )
        if(myrank == 0) r_array_one_value(1) = model_ASKI_factor_shepard_radius
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) model_ASKI_factor_shepard_radius = r_array_one_value(1)
     end select

     if(myrank == 0) i_array_one_value(1) = mAc%ncell
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mAc%ncell = i_array_one_value(1)

     if(myrank == 0) i_array_one_value(1) = mAc%max_nnb
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mAc%max_nnb = i_array_one_value(1)

     if(myrank .ne. 0) then
        allocate(mAc%cc(3,mAc%ncell),mAc%r(mAc%ncell),&
             mAc%nb(mAc%max_nnb+1,mAc%ncell))
     end if
     call bcast_all_r(mAc%cc,size(mAc%cc))
     call bcast_all_r(mAc%r,size(mAc%r))
     call bcast_all_i(mAc%nb,size(mAc%nb))

     if(myrank == 0) i_array_one_value(1) = model_ASKI_pmtrz
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) model_ASKI_pmtrz = i_array_one_value(1)

     select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLame )

        if(myrank == 0) r_array_one_value(1) = mAisoL%maxr_rho
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoL%maxr_rho = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoL%maxr_lambda
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoL%maxr_lambda = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoL%maxr_mu
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoL%maxr_mu = r_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoL%nval_rho
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoL%nval_rho = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoL%nval_lambda
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoL%nval_lambda = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoL%nval_mu
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoL%nval_mu = i_array_one_value(1)

        if(myrank .ne. 0) then
           if(mAisoL%nval_rho>0) allocate(mAisoL%idx_rho(mAisoL%nval_rho),&
                mAisoL%rho(mAisoL%nval_rho))
           if(mAisoL%nval_lambda>0) allocate(mAisoL%idx_lambda(mAisoL%nval_lambda),&
                mAisoL%lambda(mAisoL%nval_lambda))
           if(mAisoL%nval_mu>0) allocate(mAisoL%idx_mu(mAisoL%nval_mu),&
                mAisoL%mu(mAisoL%nval_mu))
        end if

        if(mAisoL%nval_rho>0) then
           call bcast_all_i(mAisoL%idx_rho,size(mAisoL%idx_rho))
           call bcast_all_r(mAisoL%rho,size(mAisoL%rho))
        end if
        if(mAisoL%nval_lambda>0) then
           call bcast_all_i(mAisoL%idx_lambda,size(mAisoL%idx_lambda))
           call bcast_all_r(mAisoL%lambda,size(mAisoL%lambda))
        end if
        if(mAisoL%nval_mu>0) then
           call bcast_all_i(mAisoL%idx_mu,size(mAisoL%idx_mu))
           call bcast_all_r(mAisoL%mu,size(mAisoL%mu))
        end if

     case ( ipmtrz_isoVelocity )

        if(myrank == 0) r_array_one_value(1) = mAisoV%maxr_rho
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoV%maxr_rho = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoV%maxr_vp
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoV%maxr_vp = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoV%maxr_vs
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoV%maxr_vs = r_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoV%nval_rho
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoV%nval_rho = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoV%nval_vp
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoV%nval_vp = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoV%nval_vs
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoV%nval_vs = i_array_one_value(1)

        if(myrank .ne. 0) then
           if(mAisoV%nval_rho>0) allocate(mAisoV%idx_rho(mAisoV%nval_rho),&
                mAisoV%rho(mAisoV%nval_rho))
           if(mAisoV%nval_vp>0) allocate(mAisoV%idx_vp(mAisoV%nval_vp),&
                mAisoV%vp(mAisoV%nval_vp))
           if(mAisoV%nval_vs>0) allocate(mAisoV%idx_vs(mAisoV%nval_vs),&
                mAisoV%vs(mAisoV%nval_vs))
        end if

        if(mAisoV%nval_rho>0) then
           call bcast_all_i(mAisoV%idx_rho,size(mAisoV%idx_rho))
           call bcast_all_r(mAisoV%rho,size(mAisoV%rho))
        end if
        if(mAisoV%nval_vp>0) then
           call bcast_all_i(mAisoV%idx_vp,size(mAisoV%idx_vp))
           call bcast_all_r(mAisoV%vp,size(mAisoV%vp))
        end if
        if(mAisoV%nval_vs>0) then
           call bcast_all_i(mAisoV%idx_vs,size(mAisoV%idx_vs))
           call bcast_all_r(mAisoV%vs,size(mAisoV%vs))
        end if
        
     case default

        ! write error message to file and stop
        write(error_message,*) 'in broadcast_ASKI_external_model: model_ASKI_pmtrz = '&
             ,model_ASKI_pmtrz,';  this parametrization index is not known: routines '//&
             'in model_external_values.f90 are inconsistent!'
        call stop_error_model_ASKI(error_message)

     end select

  end if ! impose_ASKI_inverted_model


  if(impose_ASKI_checker_model) then

     call bcast_all_ch(file_ASKI_checker_model,len_trim(file_ASKI_checker_model))

     if(myrank == 0) i_array_one_value(1) = mAchk%nchecker_radius
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mAchk%nchecker_radius = i_array_one_value(1)

     if(myrank == 0) i_array_one_value(1) = mAchk%nchecker_lat
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mAchk%nchecker_lat = i_array_one_value(1)

     if(myrank == 0) i_array_one_value(1) = mAchk%nchecker_lon
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mAchk%nchecker_lon = i_array_one_value(1)

     ! allocate for model values if I'm not rank 0
     if(myrank .ne. 0) allocate(mAchk%bchecker_radius(mAchk%nchecker_radius,2), &
          mAchk%bchecker_lat(mAchk%nchecker_lat,2),mAchk%bchecker_lon(mAchk%nchecker_lon,2))

     call bcast_all_dp(mAchk%bchecker_radius,mAchk%nchecker_radius*2)
     call bcast_all_dp(mAchk%bchecker_lat,mAchk%nchecker_lat*2)
     call bcast_all_dp(mAchk%bchecker_lon,mAchk%nchecker_lon*2)

     call bcast_all_r(mAchk%anomaly,mAchk%nparam)

     call bcast_all_dp(ASKI_Mrot_chunk1,9)
  end if ! impose_ASKI_checker_model


  end subroutine broadcast_ASKI_external_model

!
!-------------------------------------------------------------------------------------------------
!
  subroutine read_ASKI_external_checker_model()
    use shared_input_parameters,only: ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
                      CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,&
                      NCHUNKS
    use constants,only: R_EARTH_KM
    implicit none
    character(len=800) :: error_message
    character(len=MAX_STRING_LEN) :: modelfile,line
    integer :: IOASKI,ier,j,nlayer
    double precision :: wlat,gap_lat,wlon,gap_lon
    double precision, dimension(:), allocatable :: depth_upper,depth_thick
    double precision :: h,h_gap,h_cell,hmax
    logical :: file_exists
    double precision, parameter :: deg2rad = 0.017453292519943295
    double precision :: ctheta,costheta,cosphi,cosgam,sintheta,sinphi,singam
!
    if(NCHUNKS /= 1) then
       call stop_error_model_ASKI("NCHUNKS is not 1. This kind of model is only supported for 1 chunk simulations")
    end if
!
    modelfile = file_ASKI_checker_model
!
    inquire(file='DATA/'//trim(modelfile),exist=file_exists)
    if(.not.file_exists) then
       write(error_message,*) "in read_ASKI_external_checker_model: file '"//'DATA/'//&
            trim(modelfile)//"' does not exist"
       call stop_error_model_ASKI(error_message)
    end if
    ier = -1
    call get_file_unit_model_ASKI(IOASKI)
    open(unit=IOASKI,file='DATA/'//trim(modelfile), &
         status='unknown',action='read',iostat=ier)
    if(ier .ne. 0) then
       close(IOASKI) ! ?? necessary? sensible?
       ! write error message to file and stop
       write(error_message,*) "in read_ASKI_external_checker_model: could not open file '"//'DATA/'//&
            trim(modelfile)//"' to read"
       call stop_error_model_ASKI(error_message)
    end if

    write(IMAIN,*) ""
    write(IMAIN,*) "opened model file '",'DATA/'//trim(modelfile),"' for checker model successfully"

    ! ignore first line!
    read(IOASKI,*,iostat=ier) line
    if(ier/=0) then
       write(error_message,*) "in read_ASKI_external_checker_model: could not read first line ",&
            "of file'",'DATA/'//trim(modelfile),"'"
       call stop_error_model_ASKI(error_message)
    end if
    write(IMAIN,*) 'ignoring first line (comment line)'

    ! wlat,gap_lat
    read(IOASKI,*,iostat=ier) wlat,gap_lat
    if(ier/=0) then
       write(error_message,*) "in read_ASKI_external_checker_model: could not read lat width and gap values from line 2 ",&
            "of file'",'DATA/'//trim(modelfile),"'"
       call stop_error_model_ASKI(error_message)
    end if
    write(IMAIN,*) 'lat width and gap = ',wlat,gap_lat
    if(wlat < 0. .or. gap_lat < 0.) then
       write(IMAIN,*) 'ERROR: lat width or lat gap must not be negative'
       write(error_message,*) "in read_ASKI_external_checker_model: lat width or lat gap must not be negative"
       call stop_error_model_ASKI(error_message)
    end if

    ! wlon,gap_lon
    read(IOASKI,*,iostat=ier) wlon,gap_lon
    if(ier/=0) then
       write(error_message,*) "in read_ASKI_external_checker_model: could not read lon width and gap values from line 3 ",&
            "of file'",'DATA/'//trim(modelfile),"'"
       call stop_error_model_ASKI(error_message)
    end if
    write(IMAIN,*) 'lon width and gap = ',wlon,gap_lon
    if(wlon < 0. .or. gap_lon < 0.) then
       write(IMAIN,*) 'ERROR: lon width or lon gap must not be negative'
       write(error_message,*) "in read_ASKI_external_checker_model: lon width or lon gap must not be negative"
       call stop_error_model_ASKI(error_message)
    end if

    ! nlayer
    read(IOASKI,*,iostat=ier) nlayer
    if(ier/=0) then
       write(error_message,*) "in read_ASKI_external_checker_model: could not read number of checker layers from line 4 ",&
            "of file'",'DATA/'//trim(modelfile),"'"
       call stop_error_model_ASKI(error_message)
    end if
    write(IMAIN,*) 'number of checker layers = ',nlayer
    if(nlayer < 1) then
       write(IMAIN,*) 'ERROR: number of checker layers must be at least 1'
       write(error_message,*) "in read_ASKI_external_checker_model: number of checker layers be at least 1"
       call stop_error_model_ASKI(error_message)
    end if

    allocate(depth_upper(nlayer),depth_thick(nlayer))

    ! depth_upper
    read(IOASKI,*,iostat=ier) depth_upper
    if(ier/=0) then
       write(error_message,*) "in read_ASKI_external_checker_model: could not read ",nlayer," upper depth values ",&
            "from line 5 of file'",'DATA/'//trim(modelfile),"'"
       call stop_error_model_ASKI(error_message)
    end if
    write(IMAIN,*) 'upper dephts of layers = ',depth_upper
    if(any(depth_upper < 0.)) then
       write(IMAIN,*) 'ERROR: upper depths must not be negative'
       write(error_message,*) "in read_ASKI_external_checker_model: upper depths of checker layers must not be negative"
       call stop_error_model_ASKI(error_message)
    end if

    ! depth_thick
    read(IOASKI,*,iostat=ier) depth_thick
    if(ier/=0) then
       write(error_message,*) "in read_ASKI_external_checker_model: could not read ",nlayer," thickness values ",&
            "from line 6 of file'",'DATA/'//trim(modelfile),"'"
       call stop_error_model_ASKI(error_message)
    end if
    write(IMAIN,*) 'thicknesses of checker layers = ',depth_thick
    if(any(depth_thick < 0.)) then
       write(IMAIN,*) 'ERROR: thicknesses must not be negative'
       write(error_message,*) "in read_ASKI_external_checker_model: thicknesses of checker layers must not be negative"
       call stop_error_model_ASKI(error_message)
    end if

    ! mAchk%anomaly
    read(IOASKI,*,iostat=ier) mAchk%anomaly
    if(ier/=0) then
       write(error_message,*) "in read_ASKI_external_checker_model: could not read ",mAchk%nparam,&
            " real values from line 7 of file'",'DATA/'//trim(modelfile),"'"
       call stop_error_model_ASKI(error_message)
    end if
    write(IMAIN,*) 'relative model anomalies = ',mAchk%anomaly
    if(any(mAchk%anomaly < 0.)) then
       write(IMAIN,*) 'ERROR: relative model anomalies must not be negative'
       write(error_message,*) "in read_ASKI_external_checker_model: relative model anomalies must not be negative"
       call stop_error_model_ASKI(error_message)
    end if

    ! COMPUTE THE DEPTH AND LATERAL DISTRIBUTION OF CHECKERS

    ! TEST CONSISTENCY OF THE READ DEPTH VALUES AND DEFINE LAYER BOUNDARIES
    mAchk%nchecker_radius = nlayer
    allocate(mAchk%bchecker_radius(mAchk%nchecker_radius,2))
    write(IMAIN,*) "the radius boundaries of the checker layers in depth are:"
    ! treat the first layer separately
    mAchk%bchecker_radius(1,2) = 1.d0 - depth_upper(1)/R_EARTH_KM  ! max radius of this checker layer
    mAchk%bchecker_radius(1,1) = mAchk%bchecker_radius(1,2) - depth_thick(1)/R_EARTH_KM  ! min radius of this checker layer
    write(IMAIN,*) "    ",mAchk%bchecker_radius(1,1)*R_EARTH_KM,"   <->   ",mAchk%bchecker_radius(1,2)*R_EARTH_KM
    do j = 2,mAchk%nchecker_radius
       h = 1.d0 - depth_upper(j)/R_EARTH_KM
       if(h > mAchk%bchecker_radius(j-1,1)) then
          write(error_message,*) "in read_ASKI_external_checker_model: upper boundary of ",j,&
               "'th checker layer is not below the lower boundary of the ",j-1,"'th layer"
          call stop_error_model_ASKI(error_message)
       end if
       mAchk%bchecker_radius(j,2) = h  ! max radius of this checker layer
       mAchk%bchecker_radius(j,1) = mAchk%bchecker_radius(j,2) - depth_thick(j)/R_EARTH_KM  ! min radius of this checker layer
       write(IMAIN,*) "    ",mAchk%bchecker_radius(j,1)*R_EARTH_KM,"   <->   ",mAchk%bchecker_radius(j,2)*R_EARTH_KM
    end do ! j


    ! LATERAL LATITUDE DISTRIBUTION (equiangular distribution on chunk of width ANGULAR_WIDTH_ETA_IN_DEGREES)
    if(ANGULAR_WIDTH_ETA_IN_DEGREES <= 0.d0 .or. ANGULAR_WIDTH_ETA_IN_DEGREES > 90.d0) then
       write(error_message,*) "in read_ASKI_external_checker_model: ANGULAR_WIDTH_ETA_IN_DEGREES = ",&
            ANGULAR_WIDTH_ETA_IN_DEGREES," of the SPECFEM chunk is out of permitted range (must be > 0.0 and <= 90.0)"
       call stop_error_model_ASKI(error_message)
    end if
    h_gap = gap_lat/R_EARTH_KM ! angle in radians corresponding to distance wlat on surface of the Earth
    h_cell = wlat/R_EARTH_KM ! angle in radians corresponding to distance wlat on surface of the Earth
    ! FIND OUT HOW MANY CHECKERS
    mAchk%nchecker_lat = 0
    ! create a background cell in the center of the chunk, i.e. start with checker search at angle delta/2
    h = 0.5d0*h_gap
    hmax = 0.5d0* ANGULAR_WIDTH_ETA_IN_DEGREES*0.017453292519943295 ! half of latitudinal chunk width in rad
    do while(.true.)
       ! if there is no space for a complete checker cell, define a partial one (if there can 
       ! be a partial one) and exit
       if(h+h_cell > hmax) then
          if(hmax-h > 0.) then
             mAchk%nchecker_lat = mAchk%nchecker_lat + 2 ! the same cell can be found on the negative side of the chunk, hence +2
          end if
          exit
       end if

       ! otherwise, there is space for a complete checker layer
       mAchk%nchecker_lat = mAchk%nchecker_lat + 2 ! the same cell can be found on the negative side of the chunk, hence +2

       ! shift auxiliary value h to the beginning of the next potential checker cell, accounting for
       ! one background cell of width delta (i.e. + 2*delta)
       h = h + h_cell+h_gap
    end do
    write(IMAIN,*) 'there were ',mAchk%nchecker_lat," checker layers in latitudinal direction found"
    if(mAchk%nchecker_lat == 0) then
       write(error_message,*) "in read_ASKI_external_checker_model: there are no checker cells in latitudinal direcion"
       call stop_error_model_ASKI(error_message)
    end if
    allocate(mAchk%bchecker_lat(mAchk%nchecker_lat,2))
    
    ! NOW DEFINE LATITUDINAL BOUNDARIES OF CHECKER CELL
    h = 0.5d0*h_gap ! top of first checker layer
    do j = 1, mAchk%nchecker_lat/2-1
       mAchk%bchecker_lat(mAchk%nchecker_lat/2+j,1) = tan(h)
       mAchk%bchecker_lat(mAchk%nchecker_lat/2+j,2) = tan(h+h_cell)

       mAchk%bchecker_lat(mAchk%nchecker_lat/2+1-j,1) = -tan(h+h_cell)
       mAchk%bchecker_lat(mAchk%nchecker_lat/2+1-j,2) = -tan(h)

       h = h + h_cell+h_gap
    end do ! j
    ! treat the first and the last cell separately, in order to account for possible partial cells
    mAchk%bchecker_lat(mAchk%nchecker_lat,1) = tan(h)
    mAchk%bchecker_lat(mAchk%nchecker_lat,2) = min( tan(hmax) , tan(h+h_cell) )
    mAchk%bchecker_lat(1,1) = max( -tan(hmax) , -tan(h+h_cell) )
    mAchk%bchecker_lat(1,2) = -tan(h)

    write(IMAIN,*) "the  boundaries of the checker layers in latitudinal direction are (on normalized tangential plane):"
    do j = 1, mAchk%nchecker_lat
       write(IMAIN,*) "    ",mAchk%bchecker_lat(j,1),"   <->   ",mAchk%bchecker_lat(j,2)
    end do
    
    ! LATERAL LONGITUDE DISTRIBUTION (equiangular distribution on chunk of width ANGULAR_WIDTH_XI_IN_DEGREES)
    if(ANGULAR_WIDTH_XI_IN_DEGREES <= 0.d0 .or. ANGULAR_WIDTH_XI_IN_DEGREES > 90.d0) then
       write(error_message,*) "in read_ASKI_external_checker_model: ANGULAR_WIDTH_XI_IN_DEGREES = ",&
            ANGULAR_WIDTH_XI_IN_DEGREES," of the SPECFEM chunk is out of permitted range (must be > 0.0 and <= 90.0)"
       call stop_error_model_ASKI(error_message)
    end if
    h_gap = gap_lon/R_EARTH_KM ! angle in radians corresponding to distance wlat on surface of the Earth
    h_cell = wlon/R_EARTH_KM ! angle in radians corresponding to distance wlat on surface of the Earth
    ! FIND OUT HOW MANY CHECKERS
    mAchk%nchecker_lon = 0
    ! create a background cell in the center of the chunk, i.e. start with checker search at angle delta/2
    h = 0.5d0*h_gap
    hmax = 0.5d0* ANGULAR_WIDTH_XI_IN_DEGREES*0.017453292519943295 ! half of longitudinal chunk width in rad
    do while(.true.)
       ! if there is no space for a complete checker cell, define a partial one (if there can 
       ! be a partial one) and exit
       if(h+h_cell > hmax) then
          if(hmax-h > 0.) then
             mAchk%nchecker_lon = mAchk%nchecker_lon + 2 ! the same cell can be found on the negative side of the chunk, hence +2
          end if
          exit
       end if

       ! otherwise, there is space for a complete checker layer
       mAchk%nchecker_lon = mAchk%nchecker_lon + 2 ! the same cell can be found on the negative side of the chunk, hence +2

       ! shift auxiliary value h to the beginning of the next potential checker cell, accounting for
       ! one background cell of width delta (i.e. + 2*delta)
       h = h + h_cell+h_gap
    end do
    write(IMAIN,*) 'there were ',mAchk%nchecker_lon," checker layers in longitudinal direction"
    if(mAchk%nchecker_lon == 0) then
       write(error_message,*) "in read_ASKI_external_checker_model: there are no checker cells in longitudinal direction"
       call stop_error_model_ASKI(error_message)
    end if
    allocate(mAchk%bchecker_lon(mAchk%nchecker_lon,2))
    
    ! NOW DEFINE LONGITUDINAL BOUNDARIES OF CHECKER CELL
    h = 0.5d0*h_gap ! top of first checker layer
    do j = 1, mAchk%nchecker_lon/2-1
       mAchk%bchecker_lon(mAchk%nchecker_lon/2+j,1) = tan(h)
       mAchk%bchecker_lon(mAchk%nchecker_lon/2+j,2) = tan(h+h_cell)

       mAchk%bchecker_lon(mAchk%nchecker_lon/2+1-j,1) = -tan(h+h_cell)
       mAchk%bchecker_lon(mAchk%nchecker_lon/2+1-j,2) = -tan(h)

       h = h + h_cell+h_gap
    end do ! j
    ! treat the first and the last cell separately, in order to account for possible partial cells
    mAchk%bchecker_lon(mAchk%nchecker_lon,1) = tan(h)
    mAchk%bchecker_lon(mAchk%nchecker_lon,2) = min( tan(hmax) , tan(h+h_cell) )
    mAchk%bchecker_lon(1,1) = max( -tan(hmax) , -tan(h+h_cell) )
    mAchk%bchecker_lon(1,2) = -tan(h)

    write(IMAIN,*) "the  boundaries of the checker layers in longitudinal direction are (on normalized tangential plane):"
    do j = 1, mAchk%nchecker_lon
       write(IMAIN,*) "    ",mAchk%bchecker_lon(j,1),"   <->   ",mAchk%bchecker_lon(j,2)
    end do
    write(IMAIN,*) ""
    

    ! define chunk rotation matrix for testing an arbirtrary point to be located laterally in the checkers
    ! prepare values for the matrix computation below
    cosgam = dcos(GAMMA_ROTATION_AZIMUTH*deg2rad)
    singam = dsin(GAMMA_ROTATION_AZIMUTH*deg2rad)
    ctheta = 90.d0-CENTER_LATITUDE_IN_DEGREES
    costheta = dcos(ctheta*deg2rad)
    sintheta = dsin(ctheta*deg2rad)
    cosphi = dcos(CENTER_LONGITUDE_IN_DEGREES*deg2rad)
    sinphi = dsin(CENTER_LONGITUDE_IN_DEGREES*deg2rad)

    ! the (inverse) rotation matrix ASKI_Mrot_chunk1 rotates a point in Cartesian coordinates in the following way:
    !   - the center of a chunk is rotated to the "North Pole" (0,0,1) and the rotation by gamma about the local 
    !   - vertical axis of the chunk is reversed (i.e. chunk at north pole is not rotated about the z-axis)
    !   - in this reference frame, x points to local south (i.e. "lat" corresponds to x) and y points to local 
    !     east (i.e. "lon" corresponds to y)
    ASKI_Mrot_chunk1(1,1) = cosgam*costheta*cosphi-singam*sinphi
    ASKI_Mrot_chunk1(1,2) = cosgam*costheta*sinphi+singam*cosphi
    ASKI_Mrot_chunk1(1,3) = -cosgam*sintheta
    ASKI_Mrot_chunk1(2,1) = -singam*costheta*cosphi-cosgam*sinphi
    ASKI_Mrot_chunk1(2,2) = -singam*costheta*sinphi+cosgam*cosphi
    ASKI_Mrot_chunk1(2,3) = singam*sintheta
    ASKI_Mrot_chunk1(3,1) = sintheta*cosphi
    ASKI_Mrot_chunk1(3,2) = sintheta*sinphi
    ASKI_Mrot_chunk1(3,3) = costheta

  end subroutine read_ASKI_external_checker_model
!
!-------------------------------------------------------------------------------------------------
!
  subroutine read_ASKI_external_background_model()
    implicit none
  integer :: IOASKI,ier
  character(len=800) :: error_message
  integer :: maxnnodes,ilayer,inode,i
  real, dimension(:), allocatable :: d,u,wrho,wvp,wvs,wQmu,wQkappa
  character(len=MAX_STRING_LEN) :: modelfile,line
  logical :: file_exists
!
  modelfile = file_ASKI_background_model
!
  inquire(file='DATA/'//trim(modelfile),exist=file_exists)
  if(.not.file_exists) then
     write(error_message,*) "in read_ASKI_external_background_model: file '"//'DATA/'//&
          trim(modelfile)//"' does not exist"
     call stop_error_model_ASKI(error_message)
  end if
!
  ier = -1
  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file='DATA/'//trim(modelfile), &
       status='unknown',action='read',iostat=ier)
  if(ier .ne. 0) then
     close(IOASKI) ! ?? necessary? sensible?
     ! write error message to file and stop
     write(error_message,*) "in read_ASKI_external_background_model: could not open file '"//'DATA/'//&
          trim(modelfile)//"' to read"
     call stop_error_model_ASKI(error_message)
  end if


  write(IMAIN,*) "opened model file '",'DATA/'//trim(modelfile),"' for layered spline gradients successfully"

  ! ignore first line!
  read(IOASKI,*) line
  write(IMAIN,*) 'ignoring first line (comment line)'

  read(IOASKI,*) mA1Db%nlayers
  write(IMAIN,*) 'nlayers = ',mA1Db%nlayers

  allocate(mA1Db%nnodes(mA1Db%nlayers))
  read(IOASKI,*) mA1Db%nnodes
  write(IMAIN,*) 'number of nodes for each layer = ',mA1Db%nnodes

  if(minval(mA1Db%nnodes) .le. 1) &
       call exit_MPI_without_rank('in read_external_model_layered_spline_gradients: number of nodes in a layer must be at least 2')

  write(IMAIN,*) 'the model values are:'
  write(IMAIN,*) '*************************************************'

  maxnnodes = maxval(mA1Db%nnodes)
  allocate(mA1Db%depth(mA1Db%nlayers,maxnnodes),mA1Db%rho(mA1Db%nlayers,maxnnodes), &
           mA1Db%vp(mA1Db%nlayers,maxnnodes),mA1Db%vs(mA1Db%nlayers,maxnnodes),mA1Db%Qmu(mA1Db%nlayers,maxnnodes),&
           mA1Db%Qkappa(mA1Db%nlayers,maxnnodes),&
           mA1Db%sprho(mA1Db%nlayers,maxnnodes),mA1Db%spvp(mA1Db%nlayers,maxnnodes),mA1Db%spvs(mA1Db%nlayers,maxnnodes),&
           mA1Db%spQmu(mA1Db%nlayers,maxnnodes),mA1Db%spQkappa(mA1Db%nlayers,maxnnodes))
  mA1Db%depth(:,:) = 0.
  mA1Db%rho(:,:) = 0.
  mA1Db%vp(:,:) = 0.
  mA1Db%vs(:,:) = 0.
  mA1Db%Qmu(:,:) = 0.
  mA1Db%Qkappa(:,:) = 0.
  mA1Db%sprho(:,:) = 0.
  mA1Db%spvp(:,:) = 0.
  mA1Db%spvs(:,:) = 0.
  mA1Db%spQmu(:,:) = 0.
  mA1Db%spQkappa(:,:) = 0.
  
  do ilayer = 1,mA1Db%nlayers
     do inode = 1,mA1Db%nnodes(ilayer)
        read(IOASKI,*) mA1Db%depth(ilayer,inode),mA1Db%rho(ilayer,inode),mA1Db%vp(ilayer,inode),mA1Db%vs(ilayer,inode),&
             mA1Db%Qmu(ilayer,inode),mA1Db%Qkappa(ilayer,inode)
           write(IMAIN,*) mA1Db%depth(ilayer,inode),mA1Db%rho(ilayer,inode),mA1Db%vp(ilayer,inode),mA1Db%vs(ilayer,inode),&
                mA1Db%Qmu(ilayer,inode),mA1Db%Qkappa(ilayer,inode)
     end do ! inode
  enddo ! ilayer

  close(IOASKI)

!!$  write(IMAIN,*) "Hello, this is read_external_model. initiating mA1Db now."
!!$
!!$  !find out how many layers
!!$  mA1Db%nlayers = 3
!!$
!!$  write(IMAIN,*) 'mA1Db%nlayers == ',mA1Db%nlayers
!!$
!!$  !find out how many nodes in each layer
!!$  allocate(mA1Db%nnodes(mA1Db%nlayers))
!!$  mA1Db%nnodes(1) = 6
!!$  mA1Db%nnodes(2) = 24
!!$  mA1Db%nnodes(3) = 1 ! 1 means halfspace
!!$
!!$  write(IMAIN,*) 'mA1Db%nnodes == ',mA1Db%nnodes
!!$
!!$  !find out maximum number of nodes in a layer (in this case 24)
!!$  maxnnodes = maxval(mA1Db%nnodes)
!!$
!!$  write(IMAIN,*) 'maxnnodes == ',maxnnodes
!!$
!!$  allocate(mA1Db%depth(mA1Db%nlayers,maxnnodes),mA1Db%rho(mA1Db%nlayers,maxnnodes), &
!!$           mA1Db%vp(mA1Db%nlayers,maxnnodes),mA1Db%vs(mA1Db%nlayers,maxnnodes),mA1Db%sprho(mA1Db%nlayers,maxnnodes), &
!!$           mA1Db%spvp(mA1Db%nlayers,maxnnodes),mA1Db%spvs(mA1Db%nlayers,maxnnodes))
!!$
!!$  mA1Db%depth(:,:) = 0.
!!$  mA1Db%rho(:,:) = 0.
!!$  mA1Db%vp(:,:) = 0.
!!$  mA1Db%vs(:,:) = 0.
!!$  mA1Db%sprho(:,:) = 0.
!!$  mA1Db%spvp(:,:) = 0.
!!$  mA1Db%spvs(:,:) = 0.
!!$
!!$  write(IMAIN,*) 'size(mA1Db%depth) == ',size(mA1Db%depth)
!!$  write(IMAIN,*) 'size(mA1Db%rho) == ',size(mA1Db%rho)
!!$  write(IMAIN,*) 'size(mA1Db%vp) == ',size(mA1Db%vp)
!!$  write(IMAIN,*) 'size(mA1Db%vs) == ',size(mA1Db%vs)
!!$
!!$  ! set parameters (read from files actually)
!!$  mA1Db%depth(1,1) =   0.0000 
!!$  mA1Db%depth(1,2) =   0.5329 
!!$  mA1Db%depth(1,3) =   1.0657 
!!$  mA1Db%depth(1,4) =   1.5986 
!!$  mA1Db%depth(1,5) =   2.1314 
!!$  mA1Db%depth(1,6) =   2.6643 
!!$  mA1Db%depth(2,1) =   2.6643 
!!$  mA1Db%depth(2,2) =   3.2607 
!!$  mA1Db%depth(2,3) =   3.8571 
!!$  mA1Db%depth(2,4) =   4.4535 
!!$  mA1Db%depth(2,5) =   5.0499 
!!$  mA1Db%depth(2,6) =   5.6463 
!!$  mA1Db%depth(2,7) =   6.2427 
!!$  mA1Db%depth(2,8) =   6.8391 
!!$  mA1Db%depth(2,9) =   7.4355 
!!$  mA1Db%depth(2,10) =  8.0319  
!!$  mA1Db%depth(2,11) =  8.6283  
!!$  mA1Db%depth(2,12) =  9.2247  
!!$  mA1Db%depth(2,13) =  9.8210  
!!$  mA1Db%depth(2,14) = 10.4174  
!!$  mA1Db%depth(2,15) = 11.0138  
!!$  mA1Db%depth(2,16) = 11.6102  
!!$  mA1Db%depth(2,17) = 12.2066  
!!$  mA1Db%depth(2,18) = 12.8030  
!!$  mA1Db%depth(2,19) = 13.3994  
!!$  mA1Db%depth(2,20) = 13.9958  
!!$  mA1Db%depth(2,21) = 14.5922  
!!$  mA1Db%depth(2,22) = 15.1886  
!!$  mA1Db%depth(2,23) = 15.7850  
!!$  mA1Db%depth(2,24) = 16.3814  
!!$  mA1Db%depth(3,1) =  16.3814  
!!$
!!$  mA1Db%rho(1,1) =  1600.
!!$  mA1Db%rho(1,2) =  1600.
!!$  mA1Db%rho(1,3) =  1600.
!!$  mA1Db%rho(1,4) =  1600.
!!$  mA1Db%rho(1,5) =  1600.
!!$  mA1Db%rho(1,6) =  1600.
!!$  mA1Db%rho(2,1) =  1600.
!!$  mA1Db%rho(2,2) =  1600.
!!$  mA1Db%rho(2,3) =  1600.
!!$  mA1Db%rho(2,4) =  1600.
!!$  mA1Db%rho(2,5) =  1600.
!!$  mA1Db%rho(2,6) =  1600.
!!$  mA1Db%rho(2,7) =  1600.
!!$  mA1Db%rho(2,8) =  1600.
!!$  mA1Db%rho(2,9) =  1600.
!!$  mA1Db%rho(2,10) = 1600.
!!$  mA1Db%rho(2,11) = 1600.
!!$  mA1Db%rho(2,12) = 1600
!!$  mA1Db%rho(2,13) = 1600. 
!!$  mA1Db%rho(2,14) = 1600. 
!!$  mA1Db%rho(2,15) = 1600. 
!!$  mA1Db%rho(2,16) = 1600. 
!!$  mA1Db%rho(2,17) = 1600. 
!!$  mA1Db%rho(2,18) = 1600. 
!!$  mA1Db%rho(2,19) = 1600. 
!!$  mA1Db%rho(2,20) = 1600.
!!$  mA1Db%rho(2,21) = 1600.
!!$  mA1Db%rho(2,22) = 1600.
!!$  mA1Db%rho(2,23) = 1600.
!!$  mA1Db%rho(2,24) = 1600.
!!$  mA1Db%rho(3,1) =  2300.
!!$                  
!!$  mA1Db%vp(1,1) =  131.2
!!$  mA1Db%vp(1,2) =  245.3
!!$  mA1Db%vp(1,3) =  337.8
!!$  mA1Db%vp(1,4) =  408.5
!!$  mA1Db%vp(1,5) =  457.5
!!$  mA1Db%vp(1,6) =  484.8
!!$  mA1Db%vp(2,1) =  483.0
!!$  mA1Db%vp(2,2) =  493.3
!!$  mA1Db%vp(2,3) =  503.6
!!$  mA1Db%vp(2,4) =  513.9
!!$  mA1Db%vp(2,5) =  524.2
!!$  mA1Db%vp(2,6) =  534.5
!!$  mA1Db%vp(2,7) =  544.8
!!$  mA1Db%vp(2,8) =  555.1
!!$  mA1Db%vp(2,9) =  565.4
!!$  mA1Db%vp(2,10) = 575.7 
!!$  mA1Db%vp(2,11) = 586.1 
!!$  mA1Db%vp(2,12) = 596.4 
!!$  mA1Db%vp(2,13) = 606.7 
!!$  mA1Db%vp(2,14) = 617.0 
!!$  mA1Db%vp(2,15) = 627.3 
!!$  mA1Db%vp(2,16) = 637.6 
!!$  mA1Db%vp(2,17) = 647.9 
!!$  mA1Db%vp(2,18) = 658.2 
!!$  mA1Db%vp(2,19) = 668.5 
!!$  mA1Db%vp(2,20) = 678.8 
!!$  mA1Db%vp(2,21) = 689.1 
!!$  mA1Db%vp(2,22) = 699.4 
!!$  mA1Db%vp(2,23) = 709.7 
!!$  mA1Db%vp(2,24) = 720.0 
!!$  mA1Db%vp(3,1) =  3702.4
!!$
!!$  mA1Db%vs(1,1) =  71.0
!!$  mA1Db%vs(1,2) =  156.3
!!$  mA1Db%vs(1,3) =  217.9
!!$  mA1Db%vs(1,4) =  255.8
!!$  mA1Db%vs(1,5) =  270.1
!!$  mA1Db%vs(1,6) =  260.8
!!$  mA1Db%vs(2,1) =  258.8
!!$  mA1Db%vs(2,2) =  267.6
!!$  mA1Db%vs(2,3) =  276.3
!!$  mA1Db%vs(2,4) =  285.0
!!$  mA1Db%vs(2,5) =  293.8
!!$  mA1Db%vs(2,6) =  302.5
!!$  mA1Db%vs(2,7) =  311.2
!!$  mA1Db%vs(2,8) =  320.0
!!$  mA1Db%vs(2,9) =  328.7
!!$  mA1Db%vs(2,10) = 337.4 
!!$  mA1Db%vs(2,11) = 346.2 
!!$  mA1Db%vs(2,12) = 354.9 
!!$  mA1Db%vs(2,13) = 363.6 
!!$  mA1Db%vs(2,14) = 372.4 
!!$  mA1Db%vs(2,15) = 381.1 
!!$  mA1Db%vs(2,16) = 389.8 
!!$  mA1Db%vs(2,17) = 398.6 
!!$  mA1Db%vs(2,18) = 407.3 
!!$  mA1Db%vs(2,19) = 416.1 
!!$  mA1Db%vs(2,20) = 424.8 
!!$  mA1Db%vs(2,21) = 433.5 
!!$  mA1Db%vs(2,22) = 442.3 
!!$  mA1Db%vs(2,23) = 451.0 
!!$  mA1Db%vs(2,24) = 459.7 
!!$  mA1Db%vs(3,1) =  2150.8

  ! now compute the splines for each layer (and each parameter) in form of values for sprho,spvp,spvs = p = s''(xj) 
  ! (procedure, as well as notation of p,s,xj as written in "Algorithms" by Robert Sedgewick, ADDISON-WESLEY 2002, Chapter 38)
  allocate(d(maxnnodes),u(maxnnodes),wrho(maxnnodes),wvp(maxnnodes),wvs(maxnnodes),wQmu(maxnnodes),wQkappa(maxnnodes))  ! local variables

  ! for each layer calculate the second derivative of the respective spline at all nodes
  do ilayer = 1,mA1Db%nlayers
  if(mA1Db%nnodes(ilayer) .ge. 3) then 
     ! in case of mA1Db%nnodes(ilayer) == 2, the spline interpolation (in that case linear interpolation) works, as sprho,spvp,spvs = 0. initially

     ! initiate temporary variables and calculate their values
     d(:) = 0.
     u(:) = 0.
     wrho(:) = 0.
     wvp(:) = 0.
     wvs(:) = 0.
     wQmu(:) = 0.
     wQkappa(:) = 0.
     do i = 2,mA1Db%nnodes(ilayer) - 1
        d(i) = 2.*(mA1Db%depth(ilayer,i+1)-mA1Db%depth(ilayer,i-1))
     end do

     do i = 1,mA1Db%nnodes(ilayer) - 1
        u(i) = mA1Db%depth(ilayer,i+1)-mA1Db%depth(ilayer,i)
     end do

     do i = 2,mA1Db%nnodes(ilayer) - 1
        wrho(i) = 6.*((mA1Db%rho(ilayer,i+1)-mA1Db%rho(ilayer,i))/u(i) - &
                      (mA1Db%rho(ilayer,i)-mA1Db%rho(ilayer,i-1))/u(i-1))
        wvp(i) = 6.*((mA1Db%vp(ilayer,i+1)-mA1Db%vp(ilayer,i))/u(i) - &
                      (mA1Db%vp(ilayer,i)-mA1Db%vp(ilayer,i-1))/u(i-1))
        wvs(i) = 6.*((mA1Db%vs(ilayer,i+1)-mA1Db%vs(ilayer,i))/u(i) - &
                      (mA1Db%vs(ilayer,i)-mA1Db%vs(ilayer,i-1))/u(i-1))
        wQmu(i) = 6.*((mA1Db%Qmu(ilayer,i+1)-mA1Db%Qmu(ilayer,i))/u(i) - &
                      (mA1Db%Qmu(ilayer,i)-mA1Db%Qmu(ilayer,i-1))/u(i-1))
        wQkappa(i) = 6.*((mA1Db%Qkappa(ilayer,i+1)-mA1Db%Qkappa(ilayer,i))/u(i) - &
                      (mA1Db%Qkappa(ilayer,i)-mA1Db%Qkappa(ilayer,i-1))/u(i-1))
     end do

     ! now calculate the second derivatives of the spline, assuming them being zero at the extremal nodes (natural boundary conditions)
     mA1Db%sprho(ilayer,1) = 0.; mA1Db%sprho(ilayer,mA1Db%nnodes(ilayer)) = 0.
     mA1Db%spvp(ilayer,1) = 0.; mA1Db%spvp(ilayer,mA1Db%nnodes(ilayer)) = 0.
     mA1Db%spvs(ilayer,1) = 0.; mA1Db%spvs(ilayer,mA1Db%nnodes(ilayer)) = 0.
     mA1Db%spQmu(ilayer,1) = 0.; mA1Db%spQmu(ilayer,mA1Db%nnodes(ilayer)) = 0.
     mA1Db%spQkappa(ilayer,1) = 0.; mA1Db%spQkappa(ilayer,mA1Db%nnodes(ilayer)) = 0.

     ! then calculate the others by solving a tridiagonal system of equations
     if(mA1Db%nnodes(ilayer) > 3) then
        do i = 2,mA1Db%nnodes(ilayer) - 2
           wrho(i+1) = wrho(i+1) - wrho(i)*u(i)/d(i)
           wvp(i+1) = wvp(i+1) - wvp(i)*u(i)/d(i)
           wvs(i+1) = wvs(i+1) - wvs(i)*u(i)/d(i)
           wQmu(i+1) = wQmu(i+1) - wQmu(i)*u(i)/d(i)
           wQkappa(i+1) = wQkappa(i+1) - wQkappa(i)*u(i)/d(i)
           d(i+1) = d(i+1) - (u(i)**2)/d(i)
        end do
     endif

     do i = mA1Db%nnodes(ilayer)-1,2,-1
        mA1Db%sprho(ilayer,i) = (wrho(i) - u(i)*mA1Db%sprho(ilayer,i+1))/d(i)
        mA1Db%spvp(ilayer,i) = (wvp(i) - u(i)*mA1Db%spvp(ilayer,i+1))/d(i)
        mA1Db%spvs(ilayer,i) = (wvs(i) - u(i)*mA1Db%spvs(ilayer,i+1))/d(i)
        mA1Db%spQmu(ilayer,i) = (wQmu(i) - u(i)*mA1Db%spQmu(ilayer,i+1))/d(i)
        mA1Db%spQkappa(ilayer,i) = (wQkappa(i) - u(i)*mA1Db%spQkappa(ilayer,i+1))/d(i)
     end do

   end if
   end do ! ilayer

   deallocate(d,u,wrho,wvp,wvs,wQmu,wQkappa)
   ! done calculating splines

  end subroutine read_ASKI_external_background_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_ASKI_external_inverted_model()
!---
!
! ADD YOUR MODEL HERE
!
!---

  integer :: ier
  character(len=800) :: error_message

  ier = -1

  ! read file "DATA/ASKI_inverted_model" constituting the text file produced by the exportKim ASKI programm
  call read_model_ASKI_kim_export(file_ASKI_inverted_model)

  ! write log about this model to OUTPUT_FILES/output_mesher.txt
  write(IMAIN,*) "successfully read model file '",'DATA/'//trim(file_ASKI_inverted_model),"' containing the ASKI inverted model"
  select case (model_ASKI_interpolation_type)
  case(1)
     write(IMAIN,*) "will interpolate by standard shepard method"
  case(2)
     write(IMAIN,*) "will interpolate by shepard method applying an additional radius factor of ",&
          model_ASKI_factor_shepard_radius
  case default
     write(error_message,*) "in read_ASKI_external_inverted_model: model_ASKI_interpolation_type = ",&
          model_ASKI_interpolation_type,"; must be either 1 or 2, should have been checked before. ",&
          "THIS ERROR SHOULD NOT OCCURR, HENCE THIS MODULE IS INCONSISTENT!"
     call stop_error_model_ASKI(error_message)
  end select

  write(IMAIN,*) "ncell = ",mAc%ncell
  write(IMAIN,*) "maximal number of neighbours = ",mAc%max_nnb
  select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLame )
        write(IMAIN,*) "parametrization is isoLame"
        write(IMAIN,*) "nval_rho,nval_lambda,nval_mu = ",mAisoL%nval_rho,mAisoL%nval_lambda,mAisoL%nval_mu
        write(IMAIN,*) "maxr_rho,maxr_lambda,maxr_mu = ",mAisoL%maxr_rho,mAisoL%maxr_lambda,mAisoL%maxr_mu
     case ( ipmtrz_isoVelocity )
        write(IMAIN,*) "parametrization is isoVelocity"
        write(IMAIN,*) "nval_rho,nval_vp,nval_vs = ",mAisoV%nval_rho,mAisoV%nval_vp,mAisoV%nval_vs
        write(IMAIN,*) "maxr_rho,maxr_vp,maxr_vs = ",mAisoV%maxr_rho,mAisoV%maxr_vp,mAisoV%maxr_vs
  end select
  write(IMAIN,*) ""

  end subroutine read_ASKI_external_inverted_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_ASKI_kim_export(modelfile)

  character(len=*) :: modelfile

  integer :: ier,IOASKI,iline
  character(len=800) :: error_message

!
!NOW READ KIM EXPORT FILE
!

  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file='DATA/'//trim(modelfile),form='formatted', &
       status='old',action='read',iostat=ier)
  if(ier .ne. 0) then
     close(IOASKI) ! ?? necessary? sensible?
     ! write error message to file and stop
     write(error_message,*) "in read_model_ASKI_kim_export: could not open file '"//'DATA/'//&
          trim(modelfile)//"' to read (was given on second line of file '"//'DATA/'//&
          'ASKI_inverted_model'//"')"
     call stop_error_model_ASKI(error_message)
  end if

  iline = 0

  call check_pmtrz_ASKI_kim_export(IOASKI,iline,'DATA/'//trim(modelfile))

  call read_model_ASKI_kim_export_cells(IOASKI,iline,'DATA/'//trim(modelfile))

  select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLame )
        call read_external_model_ASKI_kim_export_isoLame(IOASKI,iline,'DATA/'//trim(modelfile))

     case ( ipmtrz_isoVelocity )
        call read_external_model_ASKI_kim_export_isoVelocity(IOASKI,iline,'DATA/'//trim(modelfile))

     case default
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) 'in read_model_ASKI_kim_export: model_ASKI_pmtrz = '&
             ,model_ASKI_pmtrz,';  this parametrization index is not known: routines '//&
             'in model_external_values.f90 are inconsistent!'
        call stop_error_model_ASKI(error_message)
  end select

  close(IOASKI)

  end subroutine read_model_ASKI_kim_export

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_pmtrz_ASKI_kim_export(IOASKI,iline,filename)

  integer :: IOASKI,iline
  character(len=*) :: filename

  integer :: ier,nparam
  character(len=11) :: pmtrz
  character(len=400) :: line
  character(len=800) :: error_message
  character(len=6), dimension(:), allocatable :: param

  ier = -1

  read(IOASKI,"(a400)",iostat=ier) line  ; iline = iline + 1
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: could not read line ",iline,&
          "of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if
  read(line,*,iostat=ier) pmtrz
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: could not read model parametrization "//&
          "from line ",iline," of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if

  select case(pmtrz)
  case('isoLame','isoVelocity') ! ok, do nothing
  case default
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: model parametrization '"//&
          trim(pmtrz)//"' (on line ",iline," of file '"//trim(filename)//"')"//&
          " is not supported, only 'isoLame','isoVelocity' supported so far"
     call stop_error_model_ASKI(error_message)
  end select

  read(IOASKI,"(a400)",iostat=ier) line   ; iline = iline + 1
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: could not read line ",iline,&
          "of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if
  read(line,*,iostat=ier) nparam
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: could not read number of parameters "//&
          "as first entry from line ",iline," of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if
  if(nparam < 1) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: number of parameters '",nparam,&
          "(first entry on line ",iline," of file '"//trim(filename)//&
          "') should be positive"
     call stop_error_model_ASKI(error_message)
  end if
  allocate(param(nparam))
  read(line,*,iostat=ier) nparam,param
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in check_pmtrz_ASKI_kim_export: could not read ",nparam," parameter names "//&
          "as second to ",nparam+1,"-th entry from line ",iline," of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if

  select case(pmtrz)
  case('isoLame')
     model_ASKI_pmtrz = ipmtrz_isoLame
     if(nparam/=3 .or. .not.(any(param == 'rho') .and. any(param == 'lambda') .and. any(param == 'mu') )) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in check_pmtrz_ASKI_kim_export: model parametrization isoLame "//&
             "expects the 3 parameter names (on line ",iline," of file '"//trim(filename)//&
             "'), namely 'rho','lambda','mu' (in any order)"
        call stop_error_model_ASKI(error_message)
     end if
  case('isoVelocity')
     model_ASKI_pmtrz = ipmtrz_isoVelocity
     if(nparam/=3 .or. .not.(any(param == 'rho') .and. any(param == 'vp') .and. any(param == 'vs') )) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in check_pmtrz_ASKI_kim_export: model parametrization isoVelocity "//&
             "expects the 3 parameter names (on line ",iline," of file '"//trim(filename)//&
             "'), namely 'rho','vp','vs' (in any order)"
        call stop_error_model_ASKI(error_message)
     end if
  case default
  end select

  deallocate(param)

  end subroutine check_pmtrz_ASKI_kim_export

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_ASKI_kim_export_cells(IOASKI,iline,filename)

  integer :: IOASKI,iline
  character(len=*) :: filename

  integer :: ier,icell,nnb
  character(len=400) :: line
  character(len=800) :: error_message
  real :: c1,c2,c3,r
  integer, dimension(:,:), pointer :: tmp

  ier = -1

  read(IOASKI,*,iostat=ier) mAc%ncell   ; iline = iline + 1
  if(ier .ne. 0) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in read_model_ASKI_kim_export_cells: could not read number of invgrid cells "//&
          "from line ",iline," of file '"//trim(filename)//"'"
     call stop_error_model_ASKI(error_message)
  end if
  if(mAc%ncell<1) then
     close(IOASKI)
     ! write error message to file and stop
     write(error_message,*) "in read_model_ASKI_kim_export_cells: number of invgrid cells ",mAc%ncell,&
          "(on line ",iline," of file '"//trim(filename)//"') must be positive"
     call stop_error_model_ASKI(error_message)
  end if

  mAc%max_nnb = 0
  allocate(mAc%cc(3,mAc%ncell),mAc%r(mAc%ncell),mAc%nb(mAc%max_nnb+1,mAc%ncell))
  
  do icell = 1,mAc%ncell

     read(IOASKI,"(a400)",iostat=ier) line   ; iline = iline+1
     if(ier .ne. 0) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_model_ASKI_kim_export_cells: could not read "//&
             "line ",iline," of file '"//trim(filename)//"'"
        call stop_error_model_ASKI(error_message)
     end if
     read(line,*,iostat=ier) mAc%cc(:,icell),mAc%r(icell),nnb
     if(ier .ne. 0) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_model_ASKI_kim_export_cells: could not read cell center, "//&
             "radius and number of neighbours from line ",iline," of file '"//trim(filename)//"'"
        call stop_error_model_ASKI(error_message)
     end if

     if(nnb>0) then
        mAc%nb(1,icell) = nnb
        if(nnb > mAc%max_nnb) then
           allocate(tmp(nnb+1,mAc%ncell))
           tmp(1:mAc%max_nnb+1,:) = mAc%nb
           deallocate(mAc%nb)
           mAc%nb => tmp
           nullify(tmp)
           mAc%max_nnb = nnb
        end if
        read(line,*,iostat=ier) c1,c2,c3,r,nnb,mAc%nb(2:nnb+1,icell)
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_model_ASKI_kim_export_cells: could not read ",nnb,&
                " neighbour indices after number of neighbours from line ",iline," of file '"//&
                trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
     elseif(nnb==0) then
        mAc%nb(1,icell) = 0
     else
        write(error_message,*) "in read_model_ASKI_kim_export_cells: number of neighbours ",nnb,&
             "of cell ",icell," on line ",iline," of file '"//trim(filename)//&
             "' must not be negative: must be 0 if no neighbours"
        call stop_error_model_ASKI(error_message)
     end if

  end do ! icell

  end subroutine read_model_ASKI_kim_export_cells

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_external_model_ASKI_kim_export_isoLame(IOASKI,iline,filename)

  integer :: IOASKI,iline
  character(len=*) :: filename

  integer :: ier,iparam
  character(len=800) :: error_message
  character(len=6) :: one_param
  logical :: rho_found,lambda_found,mu_found

  ier = -1

  rho_found = .false.; lambda_found = .false.; mu_found = .false.

  ! read the three blocks containing model values for rho,lambda,mu
  do iparam = 1,3

     read(IOASKI,*,iostat=ier) one_param   ; iline = iline + 1
     if(ier .ne. 0) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read parameter name"//&
             "from line ",iline," of file '"//trim(filename)//"'"
        call stop_error_model_ASKI(error_message)
     end if

     select case(one_param)
     case('rho')
        if(rho_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: parameter name 'rho' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'rho' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        rho_found = .true.
        read(IOASKI,*,iostat=ier) mAisoL%nval_rho   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read number of rho "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoL%nval_rho < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: number of rho model values ",&
                mAisoL%nval_rho," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoL%nval_rho > 0) then
           allocate(mAisoL%idx_rho(mAisoL%nval_rho))
           read(IOASKI,*,iostat=ier) mAisoL%idx_rho   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_rho,&
                   " cell indices for rho model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoL%idx_rho < 1 .or. mAisoL%idx_rho > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: there are cell indices for "//&
                   "rho model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoL%rho(mAisoL%nval_rho))
           read(IOASKI,*,iostat=ier) mAisoL%rho   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_rho,&
                   " rho model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case('lambda')
        if(lambda_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: parameter name 'lambda' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'lambda' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        lambda_found = .true.
        read(IOASKI,*,iostat=ier) mAisoL%nval_lambda   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read number of lambda "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoL%nval_lambda < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: number of lambda model values ",&
                mAisoL%nval_lambda," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoL%nval_lambda > 0) then
           allocate(mAisoL%idx_lambda(mAisoL%nval_lambda))
           read(IOASKI,*,iostat=ier) mAisoL%idx_lambda   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_lambda,&
                   " cell indices for lambda model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoL%idx_lambda < 1 .or. mAisoL%idx_lambda > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: there are cell indices for "//&
                   "lambda model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoL%lambda(mAisoL%nval_lambda))
           read(IOASKI,*,iostat=ier) mAisoL%lambda   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_lambda,&
                   " lambda model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case('mu')
        if(mu_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: parameter name 'mu' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'mu' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        mu_found = .true.
        read(IOASKI,*,iostat=ier) mAisoL%nval_mu   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read number of mu "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoL%nval_mu < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: number of mu model values ",&
                mAisoL%nval_mu," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoL%nval_mu > 0) then
           allocate(mAisoL%idx_mu(mAisoL%nval_mu))
           read(IOASKI,*,iostat=ier) mAisoL%idx_mu   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_mu,&
                   " cell indices for mu model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoL%idx_mu < 1 .or. mAisoL%idx_mu > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: there are cell indices for "//&
                   "mu model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoL%mu(mAisoL%nval_mu))
           read(IOASKI,*,iostat=ier) mAisoL%mu   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: could not read ",mAisoL%nval_mu,&
                   " mu model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case default
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_external_model_ASKI_kim_export_isoLame: parameter name '"//trim(one_param)//&
             "' on line ",iline," of file '"//trim(filename)//"' is not valid in parametrization isoVelocity"
        call stop_error_model_ASKI(error_message)

     end select

  end do ! iparam

  close(IOASKI)

  ! compute maxr_rho,maxr_lambda,maxr_mu
  if(mAisoL%nval_rho > 0) mAisoL%maxr_rho = maxval(mAc%r(mAisoL%idx_rho))
  if(mAisoL%nval_lambda > 0) mAisoL%maxr_lambda = maxval(mAc%r(mAisoL%idx_lambda))
  if(mAisoL%nval_mu > 0) mAisoL%maxr_mu = maxval(mAc%r(mAisoL%idx_mu))

  end subroutine read_external_model_ASKI_kim_export_isoLame

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_external_model_ASKI_kim_export_isoVelocity(IOASKI,iline,filename)

  integer :: IOASKI,iline
  character(len=*) :: filename

  integer :: ier,iparam
  character(len=800) :: error_message
  character(len=6) :: one_param
  logical :: rho_found,vp_found,vs_found

  ier = -1

  rho_found = .false.; vp_found = .false.; vs_found = .false.

  ! read the three blocks containing model values for rho,vp,vs
  do iparam = 1,3

     read(IOASKI,*,iostat=ier) one_param   ; iline = iline + 1
     if(ier .ne. 0) then
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read parameter name"//&
             "from line ",iline," of file '"//trim(filename)//"'"
        call stop_error_model_ASKI(error_message)
     end if

     select case(one_param)
     case('rho')
        if(rho_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: parameter name 'rho' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'rho' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        rho_found = .true.
        read(IOASKI,*,iostat=ier) mAisoV%nval_rho   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read number of rho "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoV%nval_rho < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: number of rho model values ",&
                mAisoV%nval_rho," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoV%nval_rho > 0) then
           allocate(mAisoV%idx_rho(mAisoV%nval_rho))
           read(IOASKI,*,iostat=ier) mAisoV%idx_rho   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_rho,&
                   " cell indices for rho model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoV%idx_rho < 1 .or. mAisoV%idx_rho > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: there are cell indices for "//&
                   "rho model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoV%rho(mAisoV%nval_rho))
           read(IOASKI,*,iostat=ier) mAisoV%rho   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_rho,&
                   " rho model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case('vp')
        if(vp_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: parameter name 'vp' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'vp' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        vp_found = .true.
        read(IOASKI,*,iostat=ier) mAisoV%nval_vp   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read number of vp "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoV%nval_vp < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: number of vp model values ",&
                mAisoV%nval_vp," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoV%nval_vp > 0) then
           allocate(mAisoV%idx_vp(mAisoV%nval_vp))
           read(IOASKI,*,iostat=ier) mAisoV%idx_vp   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_vp,&
                   " cell indices for vp model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoV%idx_vp < 1 .or. mAisoV%idx_vp > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: there are cell indices for "//&
                   "vp model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoV%vp(mAisoV%nval_vp))
           read(IOASKI,*,iostat=ier) mAisoV%vp   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_vp,&
                   " vp model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case('vs')
        if(vs_found) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: parameter name 'vs' "//&
                "occurs for the second time on line ",iline," of file '"//trim(filename)//&
                "', only one block for 'vs' allowed"
           call stop_error_model_ASKI(error_message)
        end if
        vs_found = .true.
        read(IOASKI,*,iostat=ier) mAisoV%nval_vs   ; iline = iline + 1
        if(ier .ne. 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read number of vs "//&
                "model values from line ",iline," of file '"//trim(filename)//"'"
           call stop_error_model_ASKI(error_message)
        end if
        if(mAisoV%nval_vs < 0) then
           close(IOASKI)
           ! write error message to file and stop
           write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: number of vs model values ",&
                mAisoV%nval_vs," on line ",iline," of file '"//trim(filename)//&
                "' must not be negative: must be 0 if no values"
           call stop_error_model_ASKI(error_message)
        elseif(mAisoV%nval_vs > 0) then
           allocate(mAisoV%idx_vs(mAisoV%nval_vs))
           read(IOASKI,*,iostat=ier) mAisoV%idx_vs   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_vs,&
                   " cell indices for vs model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
           if(any(mAisoV%idx_vs < 1 .or. mAisoV%idx_vs > mAc%ncell)) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: there are cell indices for "//&
                   "vs model values on line ",iline," of file '"//trim(filename)//&
                   "' which exceed the lower limit 1 or the upper limit ncell = ",mAc%ncell
              call stop_error_model_ASKI(error_message)
           end if
           allocate(mAisoV%vs(mAisoV%nval_vs))
           read(IOASKI,*,iostat=ier) mAisoV%vs   ; iline = iline + 1
           if(ier .ne. 0) then
              close(IOASKI)
              ! write error message to file and stop
              write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: could not read ",mAisoV%nval_vs,&
                   " vs model values from line ",iline," of file '"//trim(filename)//"'"
              call stop_error_model_ASKI(error_message)
           end if
        end if

     case default
        close(IOASKI)
        ! write error message to file and stop
        write(error_message,*) "in read_external_model_ASKI_kim_export_isoVelocity: parameter name '"//trim(one_param)//&
             "' on line ",iline," of file '"//trim(filename)//"' is not valid in parametrization isoVelocity"
        call stop_error_model_ASKI(error_message)

     end select

  end do ! iparam

  close(IOASKI)

  ! compute maxr_rho,maxr_vp,maxr_vs
  if(mAisoV%nval_rho > 0) mAisoV%maxr_rho = maxval(mAc%r(mAisoV%idx_rho))
  if(mAisoV%nval_vp > 0) mAisoV%maxr_vp = maxval(mAc%r(mAisoV%idx_vp))
  if(mAisoV%nval_vs > 0) mAisoV%maxr_vs = maxval(mAc%r(mAisoV%idx_vs))

  end subroutine read_external_model_ASKI_kim_export_isoVelocity
!
!-------------------------------------------------------------------------------------------------
!
  subroutine values_ASKI_external_model(iregion_code,xmesh,ymesh,zmesh,r,&
       vpv,vph,vsv,vsh,rho,Qmu,Qkappa,eta_aniso,dvp, &
       c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
       c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
!! ASKI FS FS: THIS ROUTINE IS CALLED FROM OUTSIDE THE MODULE IN get_model()

! given a GLL point, returns super-imposed velocity model values

    
    integer :: iregion_code

    ! GLL point Cartesian coordinates
    double precision :: xmesh,ymesh,zmesh

    double precision :: r,scaleval
    double precision :: vpv,vph,vsv,vsh
    double precision :: Qkappa,Qmu
    double precision :: rho,eta_aniso,dvp

    ! the 21 coefficients for an anisotropic medium in reduced notation
    double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
         c34,c35,c36,c44,c45,c46,c55,c56,c66

    ! local parameters
    double precision :: rho_dim,vp_dim,vs_dim
    real :: vp_real,vs_real
    real :: rho_real,Qkappa_real,Qmu_real
    logical :: values_defined

!---
!
! ADD YOUR MODEL HERE
!
!---

    if(.not.(use_ASKI_background_model .or. impose_ASKI_inverted_model .or. impose_ASKI_checker_model)) return

    ! non-dimensionalization scaling factor
    scaleval=dsqrt(PI*GRAV*RHOAV) ! time scaling (s^{-1}) is done with scaleval


    ! NOW DIMENSIONALIZE INCOMING DENSITY AND VELOCITIES FOR FURTHER USE IN THIS SUBROUTINE
    ! for now treat no transverse anisotropy, i.e. average out incoming horizontal and vertical 
    ! velocities to one value as an input for imposing the inverted model
    !   dimensionalize the velocities by factor (R_EARTH*scaleval)/1000.d0 which yields km/s
    !   and the density by RHOAV/1000.0d0 which yields g/cm^3
    rho_dim = rho*RHOAV/1000.0d0
    vp_dim = 0.5d0*(vpv+vph)*(R_EARTH*scaleval)/1000.d0
    vs_dim = 0.5d0*(vsv+vsh)*(R_EARTH*scaleval)/1000.d0
    ! Qmu, Qkappa do not need to be dimensionalized


    if(use_ASKI_background_model) then
       ! get the background velocities vp_real,vs_real, and rho_real, Qmu_real, Qkappa_real
       call model_1Dbackground_values_ASKI(r,rho_real,vp_real,vs_real,Qmu_real,Qkappa_real,values_defined)

       if(values_defined) then
          ! overwrite the dimensionalized model values (used for further processing below) 
          ! with the ones returned by routine model_1Dbackground_values_ASKI
          rho_dim = dble(rho_real)
          vp_dim = dble(vp_real)
          vs_dim = dble(vs_real)
          Qkappa = dble(Qkappa_real)
          Qmu = dble(Qmu_real)
       end if
    end if ! use_ASKI_background_model


    if(impose_ASKI_inverted_model) then
       ! no attenuation treated so far for inverted models, so simply leave values Qkappa,Qmu untouched here

       !! properly smooth out to the outside (if you are outside the reach of a cell center, check if
       !! there is a cell center in reach for twice the shepard radius, if so use all in reach and smooth out 
       !! linearly to the background (where for distance = twice the radius, the background model has full 
       !! contribution, and for distance = 1 shepard radius, the background model has no contribution)
       select case (model_ASKI_pmtrz)
       case ( ipmtrz_isoLame )
          call model_external_values_ASKI_isoLame(xmesh*R_EARTH_KM,ymesh*R_EARTH_KM,zmesh*R_EARTH_KM,rho_dim,vp_dim,vs_dim)
       case ( ipmtrz_isoVelocity )
          call model_external_values_ASKI_isoVelocity(xmesh*R_EARTH_KM,ymesh*R_EARTH_KM,zmesh*R_EARTH_KM,rho_dim,vp_dim,vs_dim)
       end select
    end if ! impose_ASKI_inverted_model


    if(impose_ASKI_checker_model) then
       call model_checker_anomalies_ASKI(xmesh,ymesh,zmesh,rho_dim,vp_dim,vs_dim,Qmu,Qkappa)
    end if ! impose_ASKI_checker_model


    ! NON-DIMENSIONALIZE THE RETURN VALUES:
    !   density rho in  (Kg/m^3)/RHOAV
    !   velocities vpv,vph,vsv,vsh in  (m/s)/(R_EARTH_m*sqrt(PI*GRAV*RHOAV))
    !
    ! constants GRAV,RHOAV are defined in constants.h, hence contained in module constants (?!)

    ! non-dimensionalize the model values
    ! assumed units of variables: rho [g/cm^3], vp,vs [km/s]
    rho=rho_dim*1000.0d0/RHOAV
    vpv=vp_dim*1000.0d0/(R_EARTH*scaleval)
    vsv=vs_dim*1000.0d0/(R_EARTH*scaleval)

    ! return the same values for horizontal and vertical velocities, 
    ! this routine cannot treat transverse anisotropy yet
    vph = vpv
    vsh = vsv

  end subroutine values_ASKI_external_model

!
!--------------------------------------------------------------------------------------------
!
  subroutine model_checker_anomalies_ASKI(xmesh,ymesh,zmesh,rho,vp,vs,Qmu,Qkappa)
  ! point coordinates
  double precision :: xmesh,ymesh,zmesh
  ! density, Vp and Vs, damping Qmu Qkappa
  double precision :: vp,vs,rho,Qmu,Qkappa

  ! local variables
  double precision :: r,x_lat,y_lon,ztmp,sign_checker_anomaly
  integer :: j,ichecker_radius,ichecker_lat,ichecker_lon,sum_mod_ichecker

  ! FOR COMPUTATIONAL OPTIMIZATION, FIRST CHECK IF THE RADIUS IS WITHIN A CHECKER LAYER
  ! IF NOT, RETURN
  r = sqrt(xmesh*xmesh + ymesh*ymesh + zmesh*zmesh)

  ichecker_radius = 0
  do j = 1,mAchk%nchecker_radius
     ! the radial checker layers are sorted from surface of earth downwards
     ! so first check if the incoming point is ABOVE the current checker layer, i.e. if it has larger 
     ! radius than upper boundary (it cannot be in any possible checker layer above, since those
     ! layers were checked before) . 
     ! If so, just exit this loop: the initialization ichecker_radius = 0 will then cause this routine 
     ! to return without imposing the relative anomalies, since this point is not within a checker cell
     if(r > mAchk%bchecker_radius(j,2)) exit

     ! secondly check if the incoming point is ABOVE the LOWER boundary of the current checker layer
     ! if so, it must be inside this checker layer. Indicate so.
     if(r >= mAchk%bchecker_radius(j,1)) then
        ichecker_radius = j
        exit
     end if

     ! if this was not the case, just iterate the loop
  end do
  ! if inside the loop no index of a checker layer was found, the point is not in a checker layer. Hence, return
  if(ichecker_radius == 0) return


  ! NOW CHECK LATERALLY, IF THE INCOMING POINT IS INSIDE A CHECKER CELL


  ! TRANSFORM INCOMING MESH POINT TO REFERENCE FRAME (chunk on the north pole, x pointing south. 
  ! additionally project x and y to tangential plane (assumed radius of the sphere is 1.0 )
  ztmp = ASKI_Mrot_chunk1(3,1)*xmesh + ASKI_Mrot_chunk1(3,2)*ymesh + ASKI_Mrot_chunk1(3,3)*zmesh
  ! additionally to rotation, project rotated x and y onto tangential plane by dividing by rotated z coordinate
  x_lat = (ASKI_Mrot_chunk1(1,1)*xmesh + ASKI_Mrot_chunk1(1,2)*ymesh + ASKI_Mrot_chunk1(1,3)*zmesh)/ztmp
  y_lon = (ASKI_Mrot_chunk1(2,1)*xmesh + ASKI_Mrot_chunk1(2,2)*ymesh + ASKI_Mrot_chunk1(2,3)*zmesh)/ztmp

  ! loop on checker cell boundaries in LAT direction:
  ichecker_lat = 0
  do j = 1,mAchk%nchecker_lat
     ! the lateral checker cells are sorted from min to max coordinates in the tangential plane
     ! so first check if the incoming point is LEFT of the current checker cell, i.e. if it has a smaller 
     ! coordinate than the lower boundary (it cannot be in any possible checker cell to the left, since those
     ! cells were checked before) . 
     ! If so, just exit this loop: the initialization ichecker_lat = 0 will then cause this routine 
     ! to return without imposing the relative anomalies, since this point is not within a checker cell
     if(x_lat < mAchk%bchecker_lat(j,1)) exit

     ! secondly check if the incoming point is LEFT the UPPER boundary of the current checker cell
     ! if so, it must be inside this checker cell. Indicate so.
     if(x_lat <= mAchk%bchecker_lat(j,2)) then
        ichecker_lat = j
        exit
     end if

     ! if this was not the case, just iterate the loop
  end do
  ! if inside the loop no index of a checker cell was found, the point is not in a checker cell in lat direction. 
  ! Hence, return
  if(ichecker_lat == 0) return

  ! loop on checker cell boundaries in LON direction:
  ichecker_lon = 0
  do j = 1,mAchk%nchecker_lon
     ! the lateral checker cells are sorted from min to max coordinates in the tangential plane
     ! so first check if the incoming point is LEFT of the current checker cell, i.e. if it has a smaller 
     ! coordinate than the lower boundary (it cannot be in any possible checker cell to the left, since those
     ! cells were checked before) . 
     ! If so, just exit this loop: the initialization ichecker_lon = 0 will then cause this routine 
     ! to return without imposing the relative anomalies, since this point is not within a checker cell
     if(y_lon < mAchk%bchecker_lon(j,1)) exit

     ! secondly check if the incoming point is LEFT the UPPER boundary of the current checker cell
     ! if so, it must be inside this checker cell. Indicate so.
     if(y_lon <= mAchk%bchecker_lon(j,2)) then
        ichecker_lon = j
        exit
     end if

     ! if this was not the case, just iterate the loop
  end do
  ! if inside the loop no index of a checker cell was found, the point is not in a checker cell in lon direction. 
  ! Hence, return
  if(ichecker_lon == 0) return


  ! IF THE CODE COMES HERE, THERE WERE INDICES ichecker_radius,ichecker_lat,ichecker_lon FOUND
  ! WHICH DEFINE AN EXPLICIT CHECKER CELL INSIDE WHICH THE INCOMING POINT IS LOCATED!
  ! HENCE, IMPOSE RELATIVE MODEL ANOMALIES

  ! first find out if this checker cell has positive or negative model anomaly
  ! the following sum is an indicator of the sign of anomaly:
  sum_mod_ichecker = mod(ichecker_radius,2) + mod(ichecker_lat,2) + mod(ichecker_lon,2)
  ! this value being odd or even can decide about the sign. this way, you get a 3D checker pattern
  if( mod(sum_mod_ichecker,2) == 1) then
     sign_checker_anomaly = 1.d0
  else
     sign_checker_anomaly = -1.d0
  end if

  ! rho anomaly
  if(mAchk%anomaly(1) > 0.) rho = rho * (1.d0 + sign_checker_anomaly*mAchk%anomaly(1)*0.01d0)

  ! vp anomaly
  if(mAchk%anomaly(2) > 0.) vp = vp * (1.d0 + sign_checker_anomaly*mAchk%anomaly(2)*0.01d0)

  ! vs anomaly
  if(mAchk%anomaly(3) > 0.) vs = vs * (1.d0 + sign_checker_anomaly*mAchk%anomaly(3)*0.01d0)

  ! Qmu anomaly
  if(mAchk%anomaly(4) > 0.) Qmu = Qmu * (1.d0 + sign_checker_anomaly*mAchk%anomaly(4)*0.01d0)

  ! Qkappa anomaly
  if(mAchk%anomaly(5) > 0.) Qkappa = Qkappa * (1.d0 + sign_checker_anomaly*mAchk%anomaly(5)*0.01d0)

  end subroutine model_checker_anomalies_ASKI

!
!--------------------------------------------------------------------------------------------
!

  subroutine model_1Dbackground_values_ASKI(r,rho,vp,vs,qmu_atten,qkappa_atten,values_defined)

! given a GLL point, do spline interpolation at depth zmax - z

  ! density, Vp and Vs
  real :: vp,vs,rho

  ! attenuation
  real :: qmu_atten,qkappa_atten

  ! success flag
  logical :: values_defined

  ! imaterial_id: (something like REGION_CODE), associated material flag, could be used as a dummy variable
  ! to indicate on which side of a discontinuity this point is
!  integer :: imaterial_id

  double precision :: r

  integer :: inode,ilayer
  double precision :: depth,t

  character(len=400) :: error_message

  ! FOR TEST/DEBUGGING OUTPUT ONLY:
  !real :: depth2,vp2,vs2,rho2,t2
  !character(len=250) :: filename
  !integer :: iz,ier
  !real :: layerthickness



  depth = (1.d0 - r)*R_EARTH_KM  ! [km]

  ! CHECK IF THE REQUESTED DEPTH IS ABOVE THE VERY BOTTOM OF THIS MODEL DEFINITION
  if(depth > mA1Db%depth(mA1Db%nlayers,mA1Db%nnodes(mA1Db%nlayers))) then
     values_defined = .false.
     return
  end if

!!$  layerthickness = mA1Db%depth(imaterial_id,mA1Db%nnodes(imaterial_id)) - mA1Db%depth(imaterial_id,1)
!!$
!!$  if( (depth < mA1Db%depth(imaterial_id,1) - 0.001*layerthickness) .or. &
!!$      (depth > mA1Db%depth(imaterial_id,mA1Db%nnodes(imaterial_id)) + 0.001*layerthickness) ) then
!!$     write(*,*) 'in model_external_values_layered_spline_gradients: at this depth, imaterial_id is wrong!   x,y,z=',x,y,z, &
!!$                ';  depth=',depth,';  imaterial_id=',imaterial_id,';  layerthickness=',layerthickness
!!$     call MPI_ABORT(MPI_COMM_WORLD, 30, ier)
!!$     stop 'in model_external_values_layered_spline_gradients: depth and material_id do not fit the layered gradients model!'
!!$  endif

  ! FIRST FIND OUT THE LAYER IN WHICH THE CURRENT POINT IS, I.E. WITHIN WHICH SPLINE INTERPOLATION SHOULD BE CONDUCTED
  do ilayer = 1,mA1Db%nlayers
     if(depth <= mA1Db%depth(ilayer,mA1Db%nnodes(ilayer))) exit
  end do
  ! after this loop, ilayer should be the index of the layer which contains the current point (each layer contains its bottom depth, the first layer contains the surface of the earth)

  values_defined = .false.
  do inode = 2,mA1Db%nnodes(ilayer)
     if(depth <= mA1Db%depth(ilayer,inode)) then
        ! interpolate values at current depth
        t = (depth - mA1Db%depth(ilayer,inode-1)) / &
            (mA1Db%depth(ilayer,inode) - mA1Db%depth(ilayer,inode-1))
        rho = t*mA1Db%rho(ilayer,inode) + (1.-t)*mA1Db%rho(ilayer,inode-1) + &
              (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
              ((t**3-t)*mA1Db%sprho(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%sprho(ilayer,inode-1))/6.
        vp = t*mA1Db%vp(ilayer,inode) + (1.-t)*mA1Db%vp(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spvp(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spvp(ilayer,inode-1))/6.
        vs = t*mA1Db%vs(ilayer,inode) + (1.-t)*mA1Db%vs(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spvs(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spvs(ilayer,inode-1))/6.
        qmu_atten = t*mA1Db%Qmu(ilayer,inode) + (1.-t)*mA1Db%Qmu(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spQmu(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spQmu(ilayer,inode-1))/6.
        qkappa_atten = t*mA1Db%Qkappa(ilayer,inode) + (1.-t)*mA1Db%Qkappa(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spQkappa(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spQkappa(ilayer,inode-1))/6.
        ! position between the nodes is found, so leave the do loop
        values_defined = .true.
        exit
     end if
  end do ! inode

  ! after this loop, there should have been found an interpolation
  ! if not, raise an error
  if(.not.values_defined) then
     write(error_message,*) "in routine model_1Dbackground_values_ASKI: at depth ",depth,&
          " there was no model value defined, although it should have been. This routine is erroneous"
     call stop_error_model_ASKI(error_message)
  end if

!! FS FS ##################################  TEST INTERPOLATION OF MODEL  #####################################
!!$if(testlayers(1) == 0) then
!!$   write(filename,"('/data/Kernel/specfem3D/OUTPUT_FILES/test_model_spline_proc',i3.3,'.dat')") mA1Db%myrank
!!$   open(unit=30,file=trim(filename),status='unknown',action='write')
!!$   do iz=0,2800
!!$      depth2 = real(iz)*0.01
!!$      ! find out layer, in which testpoint lies
!!$      do ilayer = 1,mA1Db%nlayers
!!$         layerthickness = mA1Db%depth(ilayer,mA1Db%nnodes(ilayer)) - mA1Db%depth(ilayer,1)
!!$         if(depth2 < mA1Db%depth(ilayer,mA1Db%nnodes(ilayer)) + 0.001*layerthickness) exit
!!$      end do ! ilayer
!!$      do inode = 2,mA1Db%nnodes(ilayer)
!!$         if(depth2 < mA1Db%depth(ilayer,inode) + 0.001*layerthickness) then
!!$            ! interpolate values at current depth
!!$            t2 = (depth2 - mA1Db%depth(ilayer,inode-1)) / &
!!$                 (mA1Db%depth(ilayer,inode) - mA1Db%depth(ilayer,inode-1))
!!$            rho2 = t2*mA1Db%rho(ilayer,inode) + (1.-t2)*mA1Db%rho(ilayer,inode-1) + &
!!$                 (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
!!$                 ((t2**3-t2)*mA1Db%sprho(ilayer,inode) + ((1-t2)**3-(1-t2))*mA1Db%sprho(ilayer,inode-1))/6.
!!$            vp2 = t2*mA1Db%vp(ilayer,inode) + (1.-t2)*mA1Db%vp(ilayer,inode-1) + &
!!$                 (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
!!$                 ((t2**3-t2)*mA1Db%spvp(ilayer,inode) + ((1-t2)**3-(1-t2))*mA1Db%spvp(ilayer,inode-1))/6.
!!$            vs2 = t2*mA1Db%vs(ilayer,inode) + (1.-t2)*mA1Db%vs(ilayer,inode-1) + &
!!$                 (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
!!$                 ((t2**3-t2)*mA1Db%spvs(ilayer,inode) + ((1-t2)**3-(1-t2))*mA1Db%spvs(ilayer,inode-1))/6.
!!$
!!$            write(30,*) depth2,rho2,vp2,vs2
!!$            ! position between the nodes is found, so leave the do loop
!!$            exit
!!$         end if
!!$      end do ! inode
!!$   end do ! iz
!!$   close(30)
!!$   testlayers(1) = 1
!!$endif
!! FS FS ##################################  TEST INTERPOLATION OF MODEL  #####################################

!!$if(imaterial_id == 1) then
!!$   if(testlayers(1) == 0) then
!!$      if((depth < maxval(mA1Db%depth(1,:)) + 0.0001) .and. (depth > minval(mA1Db%depth(1,:)) - 0.0001)) then
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 1. depth is OK' 
!!$      else
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 1. WRONG! depth ==',depth
!!$      endif
!!$      testlayers(1) = 1
!!$   endif
!!$elseif(imaterial_id == 2) then
!!$   if(testlayers(2) == 0) then
!!$      if((depth < maxval(mA1Db%depth(2,:)) + 0.0001) .and. (depth > minval(mA1Db%depth(2,:)) - 0.0001)) then
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 2. depth is OK' 
!!$      else
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 2. WRONG! depth ==',depth
!!$      endif
!!$      testlayers(2) = 1
!!$   endif
!!$elseif(imaterial_id == 3) then
!!$   if(testlayers(3) == 0) then
!!$      if(depth > minval(mA1Db%depth(3,:)) - 0.0001) then
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 3. depth is OK' 
!!$      else
!!$         write(*,*) 'this is model_external_values, myrank == ',myrank,'. imaterial_id == 3. WRONG! depth ==',depth
!!$      endif
!!$      testlayers(3) = 1
!!$   endif
!!$endif

  end subroutine model_1Dbackground_values_ASKI
!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_values_ASKI_isoLame(x,y,z,rho,vp,vs)

! given a GLL point, do neighbour interpolation between cell centers defined by values in mAc,mAisoL

! dummy variables
  ! point coordinates
  double precision :: x,y,z
  ! density, Vp and Vs
  double precision :: vp,vs,rho

! local variables
  double precision :: rho_bg,lambda_bg,mu_bg
  double precision :: rho_interpolated,lambda_interpolated,mu_interpolated
  real :: factor_shepard_radius,enlarge_influence,f
  integer :: n
  integer, dimension(:), pointer :: interpolation_index
  real, dimension(:), pointer :: interpolation_weight

  ! compute lambda,mu,(rho) from the incoming background velocity model
  ! whenever values are missing below, use the background values
  mu_bg = vs*vs*rho
  lambda_bg = vp*vp*rho-2.d0*mu_bg
  rho_bg = rho

  nullify(interpolation_index,interpolation_weight)

  select case(model_ASKI_interpolation_type)
  case(1); factor_shepard_radius = 2.
  case(2); factor_shepard_radius = model_ASKI_factor_shepard_radius 
  end select

  ! deal with rho
  ! first assign the background value as the interpolated one (in case that no if-clause below is entered)
  rho_interpolated = rho_bg
  if(mAisoL%nval_rho > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_rho,mAisoL%idx_rho,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        rho_interpolated = sum( dble(interpolation_weight) * mAisoL%rho(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_rho,mAisoL%idx_rho,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) rho_interpolated &
                = (1 - (f-1) )*sum( dble(interpolation_weight) * mAisoL%rho(interpolation_index) ) &
                + (f-1)*rho_bg
           deallocate(interpolation_index,interpolation_weight)
        end if
     end if
  end if

  ! deal with lambda
  ! first assign the background value as the interpolated one (in case that no if-clause below is entered)
  lambda_interpolated = lambda_bg
  if(mAisoL%nval_lambda > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_lambda,mAisoL%idx_lambda,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        lambda_interpolated = sum( dble(interpolation_weight) * mAisoL%lambda(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_lambda,mAisoL%idx_lambda,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) lambda_interpolated &
                = (1 - (f-1) )*sum( dble(interpolation_weight) * mAisoL%lambda(interpolation_index) ) &
                + (f-1)*lambda_bg
           deallocate(interpolation_index,interpolation_weight)
        end if
     end if
  end if

  ! deal with mu
  ! first assign the background value as the interpolated one (in case that no if-clause below is entered)
  mu_interpolated = mu_bg
  if(mAisoL%nval_mu > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_mu,mAisoL%idx_mu,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        mu_interpolated = sum( dble(interpolation_weight) * mAisoL%mu(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoL%nval_mu,mAisoL%idx_mu,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) mu_interpolated = (1 - (f-1) )*sum( dble(interpolation_weight) * mAisoL%mu(interpolation_index) ) &
                + (f-1)*mu_bg
           deallocate(interpolation_index,interpolation_weight)
        end if
     end if
  end if

  ! compute vp,vs,(rho) from the above interpolated isoLame values
  vp = dsqrt( (lambda_interpolated+2.*mu_interpolated) / rho_interpolated )
  vs = dsqrt( (mu_interpolated) / rho_interpolated )
  rho = rho_interpolated

  end subroutine model_external_values_ASKI_isoLame

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_values_ASKI_isoVelocity(x,y,z,rho,vp,vs)

! given a GLL point, do neighbour interpolation between cell centers defined by values in mAc,mAisoV

! dummy variables
  ! point coordinates
  double precision :: x,y,z
  ! density, Vp and Vs
  double precision :: vp,vs,rho

! local variables
  real :: factor_shepard_radius,enlarge_influence,f
  integer :: n
  integer, dimension(:), pointer :: interpolation_index
  real, dimension(:), pointer :: interpolation_weight

  nullify(interpolation_index,interpolation_weight)

  select case(model_ASKI_interpolation_type)
  case(1); factor_shepard_radius = 2.
  case(2); factor_shepard_radius = model_ASKI_factor_shepard_radius 
  end select

  ! deal with rho
  if(mAisoV%nval_rho > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_rho,mAisoV%idx_rho,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        rho = sum( dble(interpolation_weight) * mAisoV%rho(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_rho,mAisoV%idx_rho,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) rho = (1 - (f-1) )*sum( dble(interpolation_weight) * mAisoV%rho(interpolation_index) ) &
                + (f-1)*rho
           deallocate(interpolation_index,interpolation_weight)
        else
           ! do nothing, leave incoming model values as they are (incoming values contain background model)
        end if
     end if
  end if

  ! deal with vp
  if(mAisoV%nval_vp > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_vp,mAisoV%idx_vp,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        vp = sum( dble(interpolation_weight) * mAisoV%vp(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_vp,mAisoV%idx_vp,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) vp = (1 - (f-1) )*sum( dble(interpolation_weight) * mAisoV%vp(interpolation_index) ) &
                + (f-1)*vp
           deallocate(interpolation_index,interpolation_weight)
        else
           ! do nothing, leave incoming model values as they are (incoming values contain background model)
        end if
     end if
  end if

  ! deal with vs
  if(mAisoV%nval_vs > 0) then
     call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_vs,mAisoV%idx_vs,&
          n,interpolation_index,interpolation_weight)
     if(n>0) then
        ! assume that in this case, interpolation_index,interpolation_weight are associated and have size n
        vs = sum( dble(interpolation_weight) * mAisoV%vs(interpolation_index) )
        deallocate(interpolation_index,interpolation_weight)
     else
        ! check if there is a cell center in reach for twice the radius. If so, linearly interpolate
        ! with the background model: 
        !   get the fractional distance f*shepard_radius of the current point to the closest cell ( 1 <= f <= 2)
        !   do a linear interpolation with the background model, using f-1  (0 <= f-1 <= 1)
        !   assign the return model values as:  (1 - (f-1) )*shepard_interpolation + (f-1)*background
        !   i.e. for distance = twice the radius (f=2), the background model has full 
        !   contribution, and for distance = 1 shepard radius (f=1), the background model has no contribution at all

        enlarge_influence = 2.0
        call shepard_interpolation_model_ASKI(x,y,z,factor_shepard_radius,mAisoV%nval_vs,mAisoV%idx_vs,&
             n,interpolation_index,interpolation_weight,&
             enlarge_influence=enlarge_influence,radius_factor_closest_cell=f)
        if(n>0) then
           if(f<2.0) vs = (1 - (f-1) )*sum( dble(interpolation_weight) * mAisoV%vs(interpolation_index) ) &
                + (f-1)*vs
           deallocate(interpolation_index,interpolation_weight)
        else
           ! do nothing, leave incoming model values as they are (incoming values contain background model)
        end if
     end if
  end if

  end subroutine model_external_values_ASKI_isoVelocity

!
!-------------------------------------------------------------------------------------------------
!
  subroutine shepard_interpolation_model_ASKI(x,y,z,factor_radius,nidx_cell,idx_cell,&
       n_C,C_prime,w,enlarge_influence,radius_factor_closest_cell)

    double precision :: x,y,z
    real :: factor_radius
    integer :: nidx_cell,n_C
    integer, dimension(nidx_cell) :: idx_cell
    integer, dimension(:), pointer :: C_prime
    real, dimension(:), pointer :: w
    real, optional :: enlarge_influence,radius_factor_closest_cell

    double precision, dimension(nidx_cell) :: d
    logical, dimension(nidx_cell) :: larray
    integer :: iclose,i,j,n_C_tmp
    real :: r,r_prime,h,sum_s,smax_2
    real, dimension(:), pointer :: s,s_tmp
    real, dimension(:), allocatable :: t
    logical, dimension(:), allocatable :: larray2
    integer, dimension(:), pointer :: C_prime_tmp

!
! BEWARE: THIS ROUTINE MIGHT WELL BE OF NOT VERY GOOD PERFORMANCE!
! PLEASE DON'T HESITATE TO IMPROVE!
!

! nomenclature and method as in paper
!    Donald Shepard, "A two-dimensional interpolation function for irregularly-spaced data", 
!    Proceedings-1968 ACM National Conference

    n_C = 0
    nullify(C_prime,w)

    ! for given point P=(x,y,z) compute the distances d_i to all cell centers D_i
    d = dsqrt( (mAc%cc(1,idx_cell)-x)**2 + (mAc%cc(2,idx_cell)-y)**2 + (mAc%cc(3,idx_cell)-z)**2 )

    ! select those cells, for which P is within their radius
    if(present(enlarge_influence)) then
       larray = d <= mAc%r(idx_cell)*enlarge_influence
    else
       larray = d <= mAc%r(idx_cell)
    end if
    if(count(larray) == 0) return

    ! among all cells, for which P is within their radius, select the one with cell center closest to P
    iclose = minloc(d,1,larray)

    ! if radius_factor_closest_cell is present, define it here
    if(present(radius_factor_closest_cell)) then
       radius_factor_closest_cell = d(iclose)/mAc%r(idx_cell(iclose))
    end if

    ! use the closest cell to P found above (has index iclose in incoming array idx_cell), 
    ! to define a total influence radius r= factor_radius*mAc%r(idx_cell(iclose)) about P 
    ! all cells with centers within radius r (n_C many) are chosen as collection
    ! C_prime, among which the interpolation will be done
    r = factor_radius*mAc%r(idx_cell(iclose))
    if(present(enlarge_influence)) r = r * enlarge_influence
    larray = d <= r
    n_C = count(larray)
    if(n_C == 0) return
    if(n_C == 1) then
       allocate(C_prime(1),w(1))
       C_prime = minloc(d,1,larray)
       w = 1.
       return
    end if

    ! if the minimum distance to a cell is exactly zero, we have to stop here, as we 
    ! cannot conduct the computations below and return the closest cell as only interpolation point
    if(minval(d,larray) == 0.) then
       n_C = 1
       allocate(C_prime(1),w(1))
       C_prime(1) = minloc(d,1,larray)
       w(1) = 1.
       return
    end if

    ! otherwise, pack preliminary collection C_prime_tmp
    n_C_tmp = n_C
    allocate(C_prime_tmp(n_C_tmp))
    C_prime_tmp = pack( (/ (i,i=1,nidx_cell) /) , larray)

    ! choose radius r_prime, beyond which the influence of a cell center will be absolutely zero, 
    ! as the shortest distance to P of a cell center which is not in collection C_prime_tmp, 
    ! or r if all cell centers are in the collection
    if(n_C_tmp == nidx_cell) then
       r_prime = r
    else
       r_prime = minval(d,.not.larray)
    end if

    ! compute preliminary values s_i
    allocate(s_tmp(n_C_tmp))
    do i = 1,n_C_tmp
       if(d(C_prime_tmp(i)) .le. r_prime/3.) then
          s_tmp(i) = 1./d(C_prime_tmp(i)) ! this should not be a problem, since it was checked above that all d(C_prime_tmp) > 0.
       else
          h = d(C_prime_tmp(i))/r_prime - 1.
          s_tmp(i) = 27.*h*h/(4.*r_prime)
       end if
    end do ! i

    ! for reasons of numerical stability, we need to assure that the values (s_i)^2 do not
    ! cause any significant roundoff when summing them up
    ! for this reason, we neglect all cell centers i in the collection for which (s_i)^2/(s_max)^2 < epsilon(1.0)
    ! (as those will behave as zero-weight points)
    smax_2 = maxval(s_tmp); smax_2 = smax_2*smax_2
    allocate(larray2(n_C_tmp))
    larray2 = s_tmp*s_tmp/smax_2 >= epsilon(1.0)
    n_C = count(larray2) ! note that n_C > 0 as s_tmp*s_tmp/smax_2 == 1 for maxval(s_tmp)
    if(n_C < n_C_tmp) then
       allocate(C_prime(n_C),s(n_C))
       C_prime = pack(C_prime_tmp,larray2)
       s = pack(s_tmp,larray2)
       deallocate(C_prime_tmp,s_tmp)
    else
       C_prime => C_prime_tmp
       s => s_tmp
       nullify(C_prime_tmp,s_tmp)
    end if
    deallocate(larray2)

    ! if there is only a single cell center left (very close to point P), the only weight is 1.
    if(n_C == 1) then
       allocate(w(1))
       w(1) = 1.
       deallocate(s)
       return
    end if

    ! otherwise define interpolation weights for final collection C_prime from s and 
    ! direction factors t (interpolation function f_3 in the paper)

    sum_s = sum(s)

    allocate(t(n_C))
    do i = 1,n_C
       t(i) = 0.
       do j = 1,n_C
          if(j==i) cycle
          h = (mAc%cc(1,C_prime(i)) - x)*(mAc%cc(1,C_prime(j)) - x) + &
              (mAc%cc(2,C_prime(i)) - y)*(mAc%cc(2,C_prime(j)) - y) + &
              (mAc%cc(3,C_prime(i)) - z)*(mAc%cc(3,C_prime(j)) - z)
          t(i) = t(i) + s(j)*(1.-h/(d(C_prime(i))*d(C_prime(j))))
       end do ! j
       t(i) = t(i)/sum_s
    end do ! i

    allocate(w(n_C))
    w = s*s*(1.+t)
    w = w/sum(w)

    deallocate(s,t)
  end subroutine shepard_interpolation_model_ASKI
!
!---------------------------------------------------------------------------------
!
subroutine read_Par_file_ASKI()

  implicit none

  character(len=500), dimension(:), allocatable :: val_parfile
  character(len=100), dimension(:), allocatable :: key_parfile
  character(len=601) :: line
  character(len=500) :: val,error_message
  integer :: npar,ios,IOASKI,eqindx

  ! open Par_file_ASKI and find number of valid lines
  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file='DATA/'//'Par_file_ASKI',&
       form='formatted',status='old',action='read',iostat=ios)
  if(ios/=0) then
     close(IOASKI)
     ! if there is no file Par_file_ASKI:
     if(model_ASKI_myrank==0) write(IMAIN,*) "could not open file '"//'DATA/'//&
          "Par_file_ASKI', so no ASKI external models are imposed" ! actually this subroutine is only called by rank 0 ...
     use_ASKI_background_model = .false.
     impose_ASKI_inverted_model = .false.
     impose_ASKI_checker_model = .false.
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
  allocate(key_parfile(npar),val_parfile(npar))

  ! now open again and store key,val pairs of valid lines
  call get_file_unit_model_ASKI(IOASKI)
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
                key_parfile(npar) = line(1:eqindx-1)
                val_parfile(npar) = adjustl(line(eqindx+1:))
           end if
        end if
     end if
  end do
  close(IOASKI)

  ! now set values of variables use_ASKI_background_model, impose_ASKI_inverted_model

  ! USE_ASKI_BACKGROUND_MODEL
  call get_value_Par_file_ASKI('USE_ASKI_BACKGROUND_MODEL',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) use_ASKI_background_model
  if(ios/=0) call exit_MPI_without_rank("invalid logical value for parameter 'USE_ASKI_BACKGROUND_MODEL' in '"&
       //'DATA/'//"Par_file_ASKI'")

  if(use_ASKI_background_model) then
     ! FILE_ASKI_BACKGROUND_MODEL
     call get_value_Par_file_ASKI('FILE_ASKI_BACKGROUND_MODEL',val,key_parfile,val_parfile,npar)
     read(val,*,iostat=ios) file_ASKI_background_model
     if(ios/=0) call exit_MPI_without_rank("invalid character value for parameter 'FILE_ASKI_BACKGROUND_MODEL' in '"&
          //'DATA/'//"Par_file_ASKI'")
  end if

  ! IMPOSE_ASKI_INVERTED_MODEL
  call get_value_Par_file_ASKI('IMPOSE_ASKI_INVERTED_MODEL',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) impose_ASKI_inverted_model
  if(ios/=0) call exit_MPI_without_rank("invalid logical value for parameter 'IMPOSE_ASKI_INVERTED_MODEL' in '"&
       //'DATA/'//"Par_file_ASKI'")

  if(impose_ASKI_inverted_model) then
     ! FILE_ASKI_INVERTED_MODEL
     call get_value_Par_file_ASKI('FILE_ASKI_INVERTED_MODEL',val,key_parfile,val_parfile,npar)
     read(val,*,iostat=ios) file_ASKI_inverted_model
     if(ios/=0) call exit_MPI_without_rank("invalid character value for parameter 'FILE_ASKI_INVERTED_MODEL' in '"&
          //'DATA/'//"Par_file_ASKI'")

     call get_value_Par_file_ASKI('ASKI_INVERTED_MODEL_INTERPOLATION_TYPE',val,key_parfile,val_parfile,npar)
     select case (val)
     case ('shepard_standard')
        model_aski_interpolation_type = 1
     case('shepard_factor_radius')
        model_ASKI_interpolation_type = 2
     case default
        call stop_error_model_ASKI("value '"//trim(val)//"' of parameter 'ASKI_INVERTED_MODEL_INTERPOLATION_TYPE' in '"&
       //'DATA/'//"Par_file_ASKI' not supported: must be one of 'shepard_standard', 'shepard_factor_radius'")
     end select
  end if
  if(model_ASKI_interpolation_type == 2) then
     call get_value_Par_file_ASKI('ASKI_INVERTED_MODEL_FACTOR_SHEPARD_RADIUS',val,key_parfile,val_parfile,npar)
     read(val,*,iostat=ios) model_ASKI_factor_shepard_radius
     if(ios/=0) call exit_MPI_without_rank("invalid real value for parameter 'ASKI_INVERTED_MODEL_FACTOR_SHEPARD_RADIUS'"//&
          " in '"//'DATA/'//"Par_file_ASKI'")
     if(model_ASKI_factor_shepard_radius <= 0.0) then
        write(error_message,*) "value ",model_ASKI_factor_shepard_radius,&
             " of parameter 'ASKI_INVERTED_MODEL_FACTOR_SHEPARD_RADIUS'"//&
             " in '"//'DATA/'//"Par_file_ASKI' is invalid: must be strictly positve"
        call stop_error_model_ASKI(error_message)
     end if
  end if

  ! IMPOSE_ASKI_CHECKER_MODEL
  call get_value_Par_file_ASKI('IMPOSE_ASKI_CHECKER_MODEL',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) impose_ASKI_checker_model
  if(ios/=0) call exit_MPI_without_rank("invalid logical value for parameter 'IMPOSE_ASKI_CHECKER_MODEL' in '"&
       //'DATA/'//"Par_file_ASKI'")

  if(impose_ASKI_checker_model) then
     ! FILE_ASKI_CHECKER_MODEL
     call get_value_Par_file_ASKI('FILE_ASKI_CHECKER_MODEL',val,key_parfile,val_parfile,npar)
     read(val,*,iostat=ios) file_ASKI_checker_model
     if(ios/=0) call exit_MPI_without_rank("invalid character value for parameter 'FILE_ASKI_CHECKER_MODEL' in '"&
          //'DATA/'//"Par_file_ASKI'")
  end if

  if(allocated(key_parfile)) deallocate(key_parfile)
  if(allocated(val_parfile)) deallocate(val_parfile)
end subroutine read_Par_file_ASKI
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine get_value_Par_file_ASKI(key,val,key_parfile,val_parfile,npar)
  character(len=*), intent(in) :: key
  integer, intent(in) :: npar
  character(len=*), dimension(npar), intent(in) :: key_parfile,val_parfile
  character(len=500), intent(out) :: val
  integer :: ipar
  logical :: found
  found = .false.
  do ipar = 1,size(key_parfile)
     if(key == key_parfile(ipar)) then
        val = val_parfile(ipar)
        found = .true.
        exit
     end if
  end do ! ipar
  if(.not.found) call exit_MPI_without_rank("definition of parameter '"//trim(key)//"' not found in '"&
          //'DATA/'//"Par_file_ASKI'")
end subroutine get_value_Par_file_ASKI
!
!---------------------------------------------------------------------------------
!
  subroutine stop_error_model_ASKI(error_message)
  character(len=*) :: error_message
  character(len=400) :: filename
  integer :: IOASKI

  write(*,*) "myrank = ",model_ASKI_myrank," , EXTERNAL MODEL ERROR: ",trim(error_message)

  write(filename,"(a,i6.6,a)") 'OUTPUT_FILES/'//'ERROR_model_external_ASKI_',model_ASKI_myrank,'.txt'

  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=filename,form='formatted',status='unknown',action='write')
  write(IOASKI,*) trim(error_message)
  close(IOASKI)
  call abort_mpi()

  end subroutine stop_error_model_ASKI
!
!---------------------------------------------------------------------------------
!
  subroutine get_file_unit_model_ASKI(unit_out)
   integer :: unit_out
   integer :: fu
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
 end subroutine get_file_unit_model_ASKI


  end module ASKI_external_model
