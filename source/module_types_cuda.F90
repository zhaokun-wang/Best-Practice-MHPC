
!> @brief Module containing object orientation to "satisfy the eye"
module module_types_cuda
  #ifdef _CUDA_KERN
  use physical_constants
  use physical_parameters
  use parallel_parameters
  use indexing
  use legendre_quadrature
  use dimensions
  use iodir
  use mpi
  use cudafor

  implicit none

  private

  public :: reference_state
  public :: atmospheric_state
  public :: atmospheric_flux
  public :: atmospheric_tendency

  public :: assignment(=)

  !> @brief Reference state (Initial + boundary conditions)
  type reference_state
    !> Density
    real(wp), device, allocatable, dimension(:) :: density
    !> Density * Theta
    real(wp), device, allocatable, dimension(:) :: denstheta
    !> Initial density
    real(wp), device, allocatable, dimension(:) :: idens
    !> Initial density * theta
    real(wp), device, allocatable, dimension(:) :: idenstheta
    !> Pressure
    real(wp), device, allocatable, dimension(:) :: pressure
    contains
    procedure, device, public :: new_ref
    procedure, device, public :: del_ref
  end type reference_state

  !> @brief Atmospheric state to be evolved
  type atmospheric_state
    !> Backing storage for data (first two entries are x and z directions, last one is the array containing dens, umom, wmom and rhot)
    real(wp), device, pointer, dimension(:,:,:) :: mem => null( )
    !> Density
    real(wp), device, pointer, dimension(:,:) :: dens
    !> Current horizontal velocity
    real(wp), device, pointer, dimension(:,:) :: umom
    !> Current vertical velocity
    real(wp), device, pointer, dimension(:,:) :: wmom
    !> Rho theta value
    real(wp), device, pointer, dimension(:,:) :: rhot
    !> same data array in GPU
    contains
    procedure, device, public :: new_state
    procedure, device, public :: set_state
    procedure, device, public :: del_state
    procedure, device, public :: update
    procedure, device, public :: exchange_halo_x
    procedure, device, public :: exchange_halo_z
  end type atmospheric_state

  !> @brief Flux computed at volume interfaces
  type atmospheric_flux
    !> Backing storage for flux data
    real(wp), device, pointer, dimension(:,:,:) :: mem => null( )
    !> Density
    real(wp), device, pointer, dimension(:,:) :: dens
    !> Horizontal velocity
    real(wp), device, pointer, dimension(:,:) :: umom
    !> Vertical velocity
    real(wp), device, pointer, dimension(:,:) :: wmom
    !> Rho theta value
    real(wp), device, pointer, dimension(:,:) :: rhot
    contains
    procedure, device, public :: new_flux
    procedure, device, public :: set_flux
    procedure, device, public :: del_flux
  end type atmospheric_flux
  !> @brief Tendency used to update the stat
  type atmospheric_tendency
    !> Backing storage data for tendency
    real(wp), device, pointer, dimension(:,:,:) :: mem => null( )
    !> Density
    real(wp), device, pointer, dimension(:,:) :: dens
    !> Horizontal velocity
    real(wp), device, pointer, dimension(:,:) :: umom
    !> Vertical velocity
    real(wp), device, pointer, dimension(:,:) :: wmom
    !> Rho theta value
    real(wp), device, pointer, dimension(:,:) :: rhot
    contains
    procedure, device, public :: new_tendency
    procedure, device, public :: set_tendency
    procedure, device, public :: del_tendency
    procedure, device, public :: xtend
    procedure, device, public :: ztend
  end type atmospheric_tendency
  !> Interface implementing the assignment operator between atmospheric states
  interface assignment(=)
    module procedure state_equal_to_state
  end interface assignment(=)

  contains

    !> @brief Instantiates a new atmospheric state
    !> @param[inout] atmo Atmospheric state
  subroutine new_state(atmo)
    implicit none
    class(atmospheric_state), intent(inout) :: atmo
    !> If memory is already allocated, deallocate
    if ( associated(atmo%mem) ) deallocate(atmo%mem)
    !> Allocate memory and set pointers
    !PARALLEL
    allocate(atmo%mem(1-hs:nx+hs, 1-hs:nz_loc+hs, NVARS))
    atmo%dens(1-hs:,1-hs:) => atmo%mem(:,:,I_DENS)
    atmo%umom(1-hs:,1-hs:) => atmo%mem(:,:,I_UMOM)
    atmo%wmom(1-hs:,1-hs:) => atmo%mem(:,:,I_WMOM)
    atmo%rhot(1-hs:,1-hs:) => atmo%mem(:,:,I_RHOT)
  end subroutine new_state

    !> @brief Sets an existing atmospheric state to a given value
    !> @param[inout] atmo Existing atmospheric state
    !> @param[in] xval New value to be assigned to atmo
    attribute(device) subroutine set_state(atmo, xval)
      implicit none
      class(atmospheric_state), device, intent(inout) :: atmo
      real(wp), intent(in) :: xval
       
      i = blockIdx%x * blockDim%x + threadIdx%x - hs
      k = blockIdx%y * blockDim%y + threadIdx%y - hs
      ll = blockIdx%z * blockDim%z + threadIdx%z - hs
      if ( i >= 1-hs .and. i <= nx+hs k >= 1-hs .and. k <= nz_loc+hs ) then
        atmo%mem(i, k, ll) = xval
      end if
    end subroutine set_state

    !> @brief Deletes existing atmospheric state
    !> @param[inout] atmo Existing atmospheric state to be deleted
  subroutine del_state(atmo)
    implicit none
    class(atmospheric_state), intent(inout) :: atmo
    if ( associated(atmo%mem) ) then
      deallocate(atmo%mem)
    end if
    nullify(atmo%dens)
    nullify(atmo%umom)
    nullify(atmo%wmom)
    nullify(atmo%rhot)
  end subroutine del_state


    !> @brief Evolve atmospheric state according to tendency and time step
    !> param[in] so Current state to be evolve
    !> param[inout] s2 Output evolved state
    !> param[in] tend Tendency of the evolution
    !> param[in] dt Evolution time step
  attribute(device) subroutine update(s2,s0,tend,dt)
    implicit none
    class(atmospheric_state), device, intent(inout) :: s2
    class(atmospheric_state), device, intent(in) :: s0
    class(atmospheric_tendency), device, intent(in) :: tend
    real(wp), intent(in) :: dt
    if ( i >= 1 .and i <= nx .and. k >= 1 .and. k <= nz_loc ) then
      s2%mem(i,k,ll) = s0%mem(i,k,ll) + dt * tend%mem(i,k,ll)
    end if
  end subroutine update

    !> @brief Computes the atmospheric tendency along x
    !> param[inout] tendency Atmospheric tendency
    !> param[inout] flux Atmospheric flux
    !> param[in] ref Reference atmospheric state
    !> param[inout] atmostat Atmospheric stategit c
    !> param[in] dx Horizontal cell size
    !> param[in] dt Time step
  attribute(device) subroutine xtend(tendency,flux,ref,atmostat,dx,dt)
    implicit none
    class(atmospheric_tendency), device, intent(inout) :: tendency
    class(atmospheric_flux), device, intent(inout) :: flux
    class(reference_state), device, intent(in) :: ref
    class(atmospheric_state), device, intent(inout) :: atmostat
    real(wp), intent(in) :: dx, dt
    integer :: s
    real(wp) :: r, u, w, t, p, hv_coef
    real(wp), dimension(STEN_SIZE) :: stencil
    real(wp), dimension(NVARS) :: d3_vals, vals

    call atmostat%exchange_halo_x( ) !< Load the interior values into halos in x

    hv_coef = -hv_beta * dx / (16.0_wp*dt) !< hyperviscosity coeff, normalized for 4th order stencil
    if ( i >= 1 .and. i <= nx+1 .and. k >= 1 .and. k <= nz_loc ) then
      do s = 1, STEN_SIZE
        stencil(s) = atmostat%mem(i-hs-1+s,k,ll) !< Collecting neighbor values in the stencil
      end do
      vals(ll) = - 1.0_wp * stencil(1)/12.0_wp & !< Compute fluxes with 4th-order scheme
                 + 7.0_wp * stencil(2)/12.0_wp &
                 + 7.0_wp * stencil(3)/12.0_wp &
                 - 1.0_wp * stencil(4)/12.0_wp
      d3_vals(ll) = - 1.0_wp * stencil(1) & !< Numerical dissipation prop to 3rd derivative
                  + 3.0_wp * stencil(2) &
                  - 3.0_wp * stencil(3) &
                  + 1.0_wp * stencil(4)

      r = vals(I_DENS) + ref%density(k)  !< Total density
      u = vals(I_UMOM) / r               !< Velocity in x
      w = vals(I_WMOM) / r               !< Velocity in z
      t = ( vals(I_RHOT) + ref%denstheta(k) ) / r   !< Temperature
      p = c0*(r*t)**cdocv !< Equation of state, pressure
      !> Physical fluxes + dissipation terms
      flux%dens(i,k) = r*u - hv_coef*d3_vals(I_DENS)
      flux%umom(i,k) = r*u*u+p - hv_coef*d3_vals(I_UMOM)
      flux%wmom(i,k) = r*u*w - hv_coef*d3_vals(I_WMOM)
      flux%rhot(i,k) = r*u*t - hv_coef*d3_vals(I_RHOT)
      !> Compute the tendency in x through flux differences
      if (i <= nx) then
        tendency%mem(i,k,ll) = -( flux%mem(i+1,k,ll) - flux%mem(i,k,ll) ) / dx
      end if
    end if
  end subroutine xtend

    !> @brief Computes the atmospheric tendency along z
    !> param[inout] tendency Atmospheric tendency
    !> param[inout] flux Atmospheric flux
    !> param[in] ref Reference atmospheric state
    !> param[inout] atmostat Atmospheric state
    !> param[in] dz Vertical cell size
    !> param[in] dt Time step
  attribute(device) subroutine ztend(tendency,flux,ref,atmostat,dz,dt)
    implicit none
    class(atmospheric_tendency), device, intent(inout) :: tendency
    class(atmospheric_flux), device, intent(inout) :: flux
    class(reference_state), device, intent(in) :: ref
    class(atmospheric_state), device, intent(inout) :: atmostat
    real(wp), intent(in) :: dz, dt
    integer :: i, k, ll, s
    real(wp) :: r, u, w, t, p, hv_coef
    real(wp), dimension(STEN_SIZE) :: stencil
    real(wp), dimension(NVARS) :: d3_vals, vals


    call atmostat%exchange_halo_z(ref) !< Load the fixed (given by ref) interior values into halos in z

    hv_coef = -hv_beta * dz / (16.0_wp*dt) !< hyperviscosity coeff, normalized for 4th order stencil
    if ( i >= 1 .and. i <= nx+1 .and. k >= 1 .and. k <= nz_loc) then
      do s = 1, STEN_SIZE
        stencil(s) = atmostat%mem(i,k-hs-1+s,ll)
      end do
      vals(ll) = - 1.0_wp * stencil(1)/12.0_wp &  !< Compute fluxes with 4th-order scheme
               + 7.0_wp * stencil(2)/12.0_wp &
               + 7.0_wp * stencil(3)/12.0_wp &
               - 1.0_wp * stencil(4)/12.0_wp
      d3_vals(ll) = - 1.0_wp * stencil(1) &   !< Numerical dissipation prop to 3rd derivative
                    + 3.0_wp * stencil(2) & 
                    - 3.0_wp * stencil(3) &
                    + 1.0_wp * stencil(4)
      r = vals(I_DENS) + ref%idens(k)  !< Total density
      u = vals(I_UMOM) / r             !< Total velocity in x
      w = vals(I_WMOM) / r             !< Total velocity in z
      t = ( vals(I_RHOT) + ref%idenstheta(k) ) / r   !< Temperature
      p = c0*(r*t)**cdocv - ref%pressure(k)          !< Equation of state, pressure
      if ((k == 1 .and. rank == 0) .or. (k == nz_loc+1 .and. rank == size - 1)) then
        w = 0.0_wp
        d3_vals(I_DENS) = 0.0_wp
      end if
      !> Physical fluxes + dissipation terms
      flux%dens(i,k) = r*w - hv_coef*d3_vals(I_DENS)
      flux%umom(i,k) = r*w*u - hv_coef*d3_vals(I_UMOM)
      flux%wmom(i,k) = r*w*w+p - hv_coef*d3_vals(I_WMOM)
      flux%rhot(i,k) = r*w*t - hv_coef*d3_vals(I_RHOT)

      !> Compute the tendency in z through flux differences
      do j = 1, 4
        do i = 1, nx
          k = merge( ks(j), ks(j) - 1 , j < 3)
          !> Compute the tendency in z through flux differences
          tendency%mem(i,k,ll) = &
              -( flux%mem(i,k+1,ll) - flux%mem(i,k,ll) ) / dz
          if (ll == I_WMOM) then
            tendency%wmom(i,k) = tendency%wmom(i,k) - atmostat%dens(i,k)*grav
          end if
        end do
      end do
      if ( k <= nz_loc ) then
        tendency%mem(i,k,ll) = -( flux%mem(i,k+1,ll) - flux%mem(i,k,ll) ) / dz
        if (ll == I_WMOM) then
            tendency%wmom(i,k) = tendency%wmom(i,k) - atmostat%dens(i,k)*grav
        end if
      end if
    end if
  end subroutine ztend


    !> @brief Implements periodic boundary conditions along x (2 rows of halos on each side)
    !> @param[inout] s Atmospheric state whose halos should be exchanged
  attribute(device) subroutine exchange_halo_x(s)
    implicit none
    class(atmospheric_state), device, intent(inout) :: s
    integer :: k, ll

    if ( k <= nz_loc ) then
        s%mem(-1,k,ll)   = s%mem(nx-1,k,ll)
        s%mem(0,k,ll)    = s%mem(nx,k,ll)
        s%mem(nx+1,k,ll) = s%mem(1,k,ll)
        s%mem(nx+2,k,ll) = s%mem(2,k,ll)
    end if
  end subroutine exchange_halo_x

    !> @brief Fixed boundary conditions along z (0 velocity at the boundary along z, and velocity given by ref along x)
    !> @param[inout] s Atmospheric state whose halos should be exchanged
    !> @param[in] ref Reference state
    attribute(device) subroutine exchange_halo_z(s,ref)
        implicit none
        class(atmospheric_state), device, intent(inout) :: s
        class(reference_state), device, intent(in) :: ref
        integer :: i, ll
        integer(8) :: rate
        integer :: send_count = 2 * (nx + 2 * hs)

        !PARALLEL COMMUNICATION DONE AT THE BEGINNING
        call system_clock(t_comm_start)
        ! SENDRECV DOWNWARDS
        call MPI_Sendrecv(s%mem(1-hs, 1, ll), send_count, MPI_DOUBLE_PRECISION, prev_rank, 0, &
        s%mem(1-hs, -1, ll), send_count, MPI_DOUBLE_PRECISION, prev_rank, 0, &
        comm, MPI_STATUS_IGNORE, ierr &
        )

        ! SENDRECV UPWARDS
        call MPI_Sendrecv(s%mem(1-hs, nz_loc-1, ll), send_count , MPI_DOUBLE_PRECISION, next_rank, 0, &
        s%mem(1-hs, nz_loc + 1, ll), send_count , MPI_DOUBLE_PRECISION, next_rank, 0, &
        =comm, MPI_STATUS_IGNORE, ierr &
        )
        call system_clock(t_comm_end,rate)
        T_communicate = T_communicate + dble(t_comm_end-t_comm_start)/dble(rate)


    !> FIRST BOUNDARY UPDATE
    if (rank == 0) then
        if (ll == I_WMOM) then
          !> Vertical velocities are set to 0
          s%mem(i,-1,ll) = 0.0_wp
          s%mem(i,0,ll) = 0.0_wp
        else if (ll == I_UMOM) then
          !> Horizontal velocities on z boundaries are scaled depending on ref values
          s%mem(i,-1,ll)   = s%mem(i,1,ll) /  &
              ref%density(1) * ref%density(-1)
          s%mem(i,0,ll)    = s%mem(i,1,ll) /  &
              ref%density(1) * ref%density(0)
        else
          !> Copying interior values into halos
          s%mem(i,-1,ll) = s%mem(i,1,ll)
          s%mem(i,0,ll) = s%mem(i,1,ll)
        end if
    end if

    !> LAST BOUNDARY UPDATE
    if (rank == size - 1) then
        if ( condition ) then
            
        end if
        do i = 1-hs,nx+hs
          if (ll == I_WMOM) then
            !> Vertical velocities are set to 0
            s%mem(i,nz_loc+1,ll) = 0.0_wp
            s%mem(i,nz_loc+2,ll) = 0.0_wp
          else if (ll == I_UMOM) then
            !> Horizontal velocities on z boundaries are scaled depending on ref values
            s%mem(i,nz_loc+1,ll) = s%mem(i,nz_loc,ll) / &
                    ref%density(nz_loc) * ref%density(nz_loc+1)
            s%mem(i,nz_loc+2,ll) = s%mem(i,nz_loc,ll) / &
                    ref%density(nz_loc) * ref%density(nz_loc+2)
          else
            !> Copying interior values into halos
            s%mem(i,nz_loc+1,ll) = s%mem(i,nz_loc,ll)
            s%mem(i,nz_loc+2,ll) = s%mem(i,nz_loc,ll)
          end if
        end do
      end do
      !$acc end parallel loop
      !$omp end parallel do
    end if


      !Sync of ranks
     ! call MPI_BARRIER(comm, ierr)

  end subroutine exchange_halo_z

    !> @brief Instantiates a new reference state
    !> @param[inout] ref New reference state to be allocated memory to
  subroutine new_ref(ref)
    implicit none
    class(reference_state), intent(inout) :: ref
    allocate(ref%density(1-hs:nz_loc+hs))
    allocate(ref%denstheta(1-hs:nz_loc+hs))
    allocate(ref%idens(nz_loc+1))
    allocate(ref%idenstheta(nz_loc+1))
    allocate(ref%pressure(nz_loc+1))
    !$acc enter data copyin(ref)
    !$acc enter data create(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
    !$acc enter data attach(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
  end subroutine new_ref

    !> @brief Delete existing reference state
    !> @param[inout] ref Reference state whose memory should be deallocated
  subroutine del_ref(ref)
    implicit none
    class(reference_state), intent(inout) :: ref
    !$acc exit data detach(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
    !$acc exit data delete(ref%density, ref%denstheta, ref%idens, ref%idenstheta, ref%pressure)
    !$acc exit data delete(ref)
    deallocate(ref%density)
    deallocate(ref%denstheta)
    deallocate(ref%idens)
    deallocate(ref%idenstheta)
    deallocate(ref%pressure)
  end subroutine del_ref

    !> @brief Instantiates a new flux object
    !> @param[inout] flux Flux object which should be initialized
  subroutine new_flux(flux)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
    if ( associated(flux%mem) ) deallocate(flux%mem)
    allocate(flux%mem(1:nx+1, 1:nz_loc+1,NVARS))
    flux%dens => flux%mem(:,:,I_DENS)
    flux%umom => flux%mem(:,:,I_UMOM)
    flux%wmom => flux%mem(:,:,I_WMOM)
    flux%rhot => flux%mem(:,:,I_RHOT)
    !$acc enter data copyin(flux)
    !$acc enter data create(flux%mem)
    !$acc enter data attach(flux%mem)
  end subroutine new_flux

    !> @brief Set an existing flux object to a given value
    !> @param[inout] flux Flux object whose value should be reassigned
    !> @param[int] xval New value to be assigned
    subroutine set_flux(flux, xval)
      implicit none
      class(atmospheric_flux), intent(inout) :: flux
      real(wp), intent(in) :: xval
      integer :: i, k, ll

      if ( .not. associated(flux%mem) ) then
        write(stderr,*) 'NOT ALLOCATED FLUX ERROR'
        stop
      end if

      !$acc parallel loop collapse(3) present(flux%mem)
      do ll = 1, NVARS
        do k = 1, nz_loc+1
          do i = 1, nx+1
            flux%mem(i,k,ll) = xval
          end do
        end do
      end do
      !$acc end parallel loop
    end subroutine set_flux

    !> @brief Deallocate an existing flux object
    !> @param[inout] flux Object which should be deallocated
  subroutine del_flux(flux)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
    if ( associated(flux%mem) ) then
      !$acc exit data detach(flux%mem)
      !$acc exit data delete(flux%mem)
      !$acc exit data delete(flux)
      deallocate(flux%mem)
    end if
    nullify(flux%dens)
    nullify(flux%umom)
    nullify(flux%wmom)
    nullify(flux%rhot)
  end subroutine del_flux

    !> @brief Allocate memory to a tendency object
    !> @param[inout] tend New tendency object which should be initialized
  subroutine new_tendency(tend)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
    if ( associated(tend%mem) ) deallocate(tend%mem)
    allocate(tend%mem(nx, nz_loc,NVARS))
    tend%dens => tend%mem(:,:,I_DENS)
    tend%umom => tend%mem(:,:,I_UMOM)
    tend%wmom => tend%mem(:,:,I_WMOM)
    tend%rhot => tend%mem(:,:,I_RHOT)
    !$acc enter data copyin(tend)
    !$acc enter data create(tend%mem)
    !$acc enter data attach(tend%mem)
  end subroutine new_tendency

    !> @brief Set an existing tendency object to a given value
    !> @param[inout] tend Tendency object whose value should be reassigned
    !> @param[int] xval New value to be assigned
    subroutine set_tendency(tend, xval)
      implicit none
      class(atmospheric_tendency), intent(inout) :: tend
      real(wp), intent(in) :: xval
      integer :: i, k, ll

      if ( .not. associated(tend%mem) ) then
        write(stderr,*) 'NOT ALLOCATED TENDENCY ERROR AT LINE ', __LINE__
        stop
      end if

      !$acc parallel loop collapse(3) present(tend%mem)
      do ll = 1, NVARS
        do k = 1, nz_loc
          do i = 1, nx
            tend%mem(i,k,ll) = xval
          end do
        end do
      end do
      !$acc end parallel loop
    end subroutine set_tendency

    !> @brief Deallocate an existing tendency object
    !> @param[inout] tend Object which should be deallocated
  subroutine del_tendency(tend)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
    if ( associated(tend%mem) ) then
      !$acc exit data detach(tend%mem)
      !$acc exit data delete(tend%mem)
      !$acc exit data delete(tend)
      deallocate(tend%mem)
    end if
    nullify(tend%dens)
    nullify(tend%umom)
    nullify(tend%wmom)
    nullify(tend%rhot)
  end subroutine del_tendency

    !> @brief Assigment operator for atmospheric states
    !> @param[inout] x Atmospheric state to be assigned a new value to
    !> @param[in] y Atmospheric state whose value is taken to be assigned
    subroutine state_equal_to_state(x,y)
      implicit none
      type(atmospheric_state), intent(inout) :: x
      type(atmospheric_state), intent(in) :: y
      integer :: i, k, ll

      !$acc parallel loop collapse(3) present(x%mem, y%mem)
      do ll = 1, NVARS
        do k = 1-hs, nz_loc+hs
          do i = 1-hs, nx+hs
            x%mem(i,k,ll) = y%mem(i,k,ll)
          end do
        end do
      end do
      !$acc end parallel loop
    end subroutine state_equal_to_state
#endif
end module module_types_cuda
