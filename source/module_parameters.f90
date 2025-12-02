!>
!! @authors Pedde, Zhaokun, Veraldi
!! @date 1-12-2025
!! @brief module_parameters.f90 contains all the parameters of the simulation


!>
!! @brief calculation_types is a collection of the type for all the code
!! @var wp   :: bit for the real number in the simulation
!! @var iowp :: type of the real numbers for the I/O interface with netCDF
!! @var ip   :: bit for the integer numbers in the simulation
module calculation_types
  use netcdf, only : nf90_real, nf90_double
  use iso_fortran_env
  use mpi
  implicit none
  public
  integer, parameter :: wp = real64
  integer, parameter :: iowp = nf90_double
  integer, parameter :: ip = int32
end module calculation_types

!>
!! @brief physical_constants is a collection of all the physical constants of the simulation
module physical_constants
  use calculation_types, only : wp
  implicit none
  public
  real(wp), parameter :: pi = 3.14159265358979323846264338327_wp !>pi greek
  real(wp), parameter :: hpi = 0.5_wp * pi                       !>half pi greek
  real(wp), parameter :: grav = 9.80665_wp                       !>gravity
  real(wp), parameter :: boltzk = 1.3806490e-23_wp               !>boltzmann
  real(wp), parameter :: navgdr = 6.02214076e23_wp               !>avogadro number
  real(wp), parameter :: rgasmol = navgdr*boltzk                 !>ideal gas const R
  real(wp), parameter :: amd = 28.96454_wp                       !>gr/mol dry air
  real(wp), parameter :: amw = 18.01528_wp                       !>gr/mol wet air
  real(wp), parameter :: rgas = (rgasmol/amd)*1000.0_wp          !>R for dry air
  real(wp), parameter :: cp = 3.5_wp*rgas                        !>specific heat const pressure
  real(wp), parameter :: cv = 2.5_wp*rgas                        !>specifc heat const volume
  real(wp), parameter :: rd = rgas                               !>ALIAS FOR R FOR DRY AIR
  real(wp), parameter :: p0 = 1.e5_wp                            !>REFERENCE PRESSURE
  real(wp), parameter :: t0 = 298.0_wp                           !>REFERENCE TEMPERATURE
  real(wp), parameter :: cdocv = 3.5_wp/2.5_wp                   !>adiabatic index gamma
  real(wp), parameter :: cvocd = 2.5_wp/3.5_wp                   !>1/gamma (inverse adiabatic)
  real(wp), parameter :: rdocp = rd/cp                           !>R_dry/cp
  real(wp), parameter :: rdocv = rd/cv                           !>R_dry/cv
  real(wp), parameter :: c0 = (rd**cdocv)*p0**(-rdocv)           !>normalization factor
  real(wp), parameter :: theta0 = 300.0_wp                       !>REFERENCE POTENTIAL TEMPERATURE (PERT BASE)
  real(wp), parameter :: exner0 = 1.0_wp                         !>REFERENCE ADIMENSIONAL PRESSURE (exner pr)
end module physical_constants


!>
!! @brief physical_parameters is a collection of all parameters related to the box and perturbation
module physical_parameters
  use calculation_types, only : wp
  implicit none
  public
  real(wp), parameter :: xlen = 20000.0_wp                       !>width of simulation
  real(wp), parameter :: zlen = 10000.0_wp                       !>height of simulation
  real(wp), parameter :: hxlen = 0.5_wp*xlen                     !>half x lenght
  real(wp), parameter :: p1 = xlen/10.0_wp                       !>10% OF DOMAIN PERTURBATION
  real(wp), parameter :: p2 = 4.0_wp*p1                          !>40% OF DOMAIN PERTURBATION
  real(wp), parameter :: p3 = 2.0_wp*p1                          !>20% OF DOMAIN PERTURBATION
  real(wp), parameter :: mnt_width = xlen/8.0_wp                 !>WIDTH OF PEAK OF GAUSSIAN
  real(wp), parameter :: hv_beta = 0.05_wp                       !>hyper viscosity coefficient
  real(wp), parameter :: cfl = 1.5_wp                            !>CFL stability number (dt comp)
  real(wp), parameter :: max_speed = 450.0_wp                    !>max velocity expected
  real(wp), parameter :: x0 = mnt_width                          !>X COORDINATE PERTURBATION
  real(wp), parameter :: z0 = 1000.0_wp                          !>Z COORDINATE PERTURBATION
  real(wp), parameter :: xrad = 500.0_wp                         !>XRADIUS OF PERTURBATION (ELLIPSE)
  real(wp), parameter :: zrad = 500.0_wp                         !>ZRADIUS OF PERTURBATION (ELLIPSE)
  real(wp), parameter :: amp = 0.01_wp                           !>AMPLITUDE OF PERTURBATION IN TEMPERATURE POT
end module physical_parameters

!>
!! @brief indexing is the memory mapping of the code
module indexing
  implicit none
  public
  integer, parameter :: NVARS = 4                                !>number of variables to resolve
  integer, parameter :: STEN_SIZE = 4                            !>accuracy of derivative (need 4 elements)
  integer, parameter :: I_DENS  = 1                              !>density indx
  integer, parameter :: I_UMOM  = 2                              !>linear moment along x (\rho u) indx
  integer, parameter :: I_WMOM  = 3                              !>linear moment along z (\rho w) indx
  integer, parameter :: I_RHOT  = 4                              !>potential temp. density (\rho \theta) indx
  integer, parameter :: DIR_X = 1                                !>x direction indx
  integer, parameter :: DIR_Z = 2                                !>z direction indx
end module indexing

!>
!! @brief iodir contains the variables for input/output/error redirection
module iodir
  use iso_fortran_env
  public
  integer, public, parameter :: stdin = input_unit
  integer, public, parameter :: stdout = output_unit
  integer, public, parameter :: stderr = error_unit
end module iodir

!>
!! @brief parallel_parameters contains the parameters for parallelization
!! @var i_beg :: is the beginning index of physical array (NOT THE HALO WILL BE IN NEGATIVE INDEX AND OUT IND)
!! @var k_beg :: is the beginning index of physical array ...
!! @var hs    :: is the dimension left/right of halo (column to exchange)-> two left two right
module parallel_parameters
  implicit none
  public
  integer, parameter :: i_beg = 1
  integer, parameter :: k_beg = 1
  integer, parameter :: hs = 2
  integer :: ierr, rank, size, comm, prev_rank, next_rank !< Initialized in model.f90
  integer :: z_local, z_global, nz_loc, base, rest !< Initialized in init
end module parallel_parameters

!>
!! @brief legendre_quadrature contains the parameters for the legendre operator
!! @var nqpoints :: number of point for the quadrature
!! @var qpoints  :: value of the quadrature points
!! @var qweights :: weights for the quadrature
module legendre_quadrature
  use calculation_types, only : wp
  implicit none
  public
  integer , parameter :: nqpoints = 3
  real(wp), parameter :: qpoints(nqpoints) =    &
      [ 0.112701665379258311482073460022E0_wp , &
        0.500000000000000000000000000000E0_wp , &
        0.887298334620741688517926539980E0_wp ]
  real(wp), parameter :: qweights(nqpoints) =   &
      [ 0.277777777777777777777777777779E0_wp , &
        0.444444444444444444444444444444E0_wp , &
        0.277777777777777777777777777779E0_wp ]
end module legendre_quadrature

!>
!! @brief dimensions contains the grid dimensions
!! @var nx          :: number of cells in x axis
!! @var nz          :: number of cells in the z axis
!! @var sim_time    :: the entire time of simulation
!! @var output_freq :: the frequency of output in time unit
module dimensions
  use calculation_types, only : wp
  use physical_parameters, only : zlen, xlen
  use indexing
  implicit none
  public
  integer , parameter :: nx = 100
  integer , parameter :: nz = int(nx * zlen/xlen)
  real(wp), parameter :: sim_time = 1000.0_wp
  real(wp), parameter :: output_freq = 10.0_wp
end module dimensions
