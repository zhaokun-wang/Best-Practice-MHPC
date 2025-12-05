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
  public :: atmospheric_state_host
  public :: atmospheric_flux
  public :: atmospheric_tendency
  public :: exchange_halo_z_mpi
  public :: exchange_halo_z_device

  public :: assignment(=)

  !> @brief Reference state (Initial + boundary conditions)
  type reference_state
    real(wp), device, allocatable, dimension(:) :: density
    real(wp), device, allocatable, dimension(:) :: denstheta
    real(wp), device, allocatable, dimension(:) :: idens
    real(wp), device, allocatable, dimension(:) :: idenstheta
    real(wp), device, allocatable, dimension(:) :: pressure
    contains
    procedure, public :: new_ref
    procedure, public :: del_ref
  end type reference_state

  !> @brief Atmospheric state to be evolved
  type atmospheric_state
    real(wp), device, pointer, dimension(:,:,:) :: mem => null()
    real(wp), device, pointer, dimension(:,:) :: dens
    real(wp), device, pointer, dimension(:,:) :: umom
    real(wp), device, pointer, dimension(:,:) :: wmom
    real(wp), device, pointer, dimension(:,:) :: rhot
    contains
    procedure, public :: new_state
    procedure, public :: set_state
    procedure, public :: del_state
    procedure, public :: update
    procedure, public :: exchange_halo_x
  end type atmospheric_state

  
  !type atmospheric_state_host
  !  real(wp), pointer, dimension(:,:,:) :: mem => null()
  !  real(wp), pointer, dimension(:,:) :: dens
  !  real(wp), pointer, dimension(:,:) :: umom
  !  real(wp), pointer, dimension(:,:) :: wmom
  !  real(wp), pointer, dimension(:,:) :: rhot
  !  contains
  !  procedure, public :: new_state
  !  procedure, public :: set_state
  ! procedure, public :: del_state
  !  procedure, public :: update
  !  procedure, public :: exchange_halo_x
  !end type atmospheric_state

  !> @brief Flux computed at volume interfaces
  type atmospheric_flux
    real(wp), device, pointer, dimension(:,:,:) :: mem => null()
    real(wp), device, pointer, dimension(:,:) :: dens
    real(wp), device, pointer, dimension(:,:) :: umom
    real(wp), device, pointer, dimension(:,:) :: wmom
    real(wp), device, pointer, dimension(:,:) :: rhot
    contains
    procedure, public :: new_flux
    procedure, public :: set_flux
    procedure, public :: del_flux
  end type atmospheric_flux

  !> @brief Tendency used to update the state
  type atmospheric_tendency
    real(wp), device, pointer, dimension(:,:,:) :: mem => null()
    real(wp), device, pointer, dimension(:,:) :: dens
    real(wp), device, pointer, dimension(:,:) :: umom
    real(wp), device, pointer, dimension(:,:) :: wmom
    real(wp), device, pointer, dimension(:,:) :: rhot
    contains
    procedure, public :: new_tendency
    procedure, public :: set_tendency
    procedure, public :: del_tendency
    procedure, public :: xtend
    procedure, public :: ztend
  end type atmospheric_tendency

  interface assignment(=)
    module procedure state_equal_to_state
  end interface

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
    attributes(device) subroutine set_state(atmo, xval)
      implicit none
      class(atmospheric_state), device, intent(inout) :: atmo
      real(wp), intent(in) :: xval
      integer :: i, k, ll

      i = blockIdx%x * blockDim%x + threadIdx%x - hs
      k = blockIdx%y * blockDim%y + threadIdx%y - hs
      ll = blockIdx%z * blockDim%z + threadIdx%z - hs
      if ( i >= 1-hs .and. i <= nx+hs .and. k >= 1-hs .and. k <= nz_loc+hs ) then
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
  attributes(device) subroutine update(s2,s0,tend,dt)
    implicit none
    class(atmospheric_state), device, intent(inout) :: s2
    class(atmospheric_state), device, intent(in) :: s0
    class(atmospheric_tendency), device, intent(in) :: tend
    real(wp), intent(in) :: dt
    integer :: i, k, ll
    i = blockIdx%x * blockDim%x + threadIdx%x
    k = blockIdx%y * blockDim%y + threadIdx%y
    ll = blockIdx%z * blockDim%z + threadIdx%z
    
    if ( i >= 1 .and. i <= nx .and. k >= 1 .and. k <= nz_loc ) then
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
  attributes(device) subroutine xtend(tendency,flux,ref,atmostat,dx,dt)
    implicit none
    class(atmospheric_tendency), device, intent(inout) :: tendency
    class(atmospheric_flux), device, intent(inout) :: flux
    class(reference_state), device, intent(in) :: ref
    class(atmospheric_state), device, intent(inout) :: atmostat
    real(wp), intent(in) :: dx, dt
    integer :: s, i, k, ll
    real(wp) :: r, u, w, t, p, hv_coef
    real(wp), dimension(STEN_SIZE) :: stencil
    real(wp), dimension(NVARS) :: d3_vals, vals
    
    i = blockIdx%x * blockDim%x + threadIdx%x - hs
    k = blockIdx%y * blockDim%y + threadIdx%y
    ll = blockIdx%z * blockDim%z + threadIdx%z
    
    call exchange_halo_x(atmostat)
    
    hv_coef = -hv_beta * dx / (16.0_wp*dt)
    if (i >= 1 .and. i <= nx+1 .and. k >= 1 .and. k <= nz_loc) then
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
      if (ll == 1) then
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
      end if

      if (i <= nx) then
        tendency%mem(i,k,ll) = -( flux%mem(i+1,k,ll) - flux%mem(i,k,ll) ) / dx
      end if
    endif
  end subroutine xtend

  attributes(device) subroutine ztend(tendency,flux,ref,atmostat,dz,dt)
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
    
    i = blockIdx%x * blockDim%x + threadIdx%x
    k = blockIdx%y * blockDim%y + threadIdx%y
    ll = blockIdx%z * blockDim%z + threadIdx%z
    
    call exchange_halo_z_device(atmostat, ref)
    
    hv_coef = -hv_beta * dz / (16.0_wp*dt)
    if (i >= 1 .and. i <= nx.and. k >= 1 .and. k <= nz_loc+1) then
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
      if (ll == 1) then
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
      end if
      tendency%mem(i,k,ll) = -( flux%mem(i,k+1,ll) - flux%mem(i,k,ll) ) / dz
      if (ll == I_WMOM) then
        tendency%wmom(i,k) = tendency%wmom(i,k) - atmostat%dens(i,k)*grav
      end if
    endif
  end subroutine ztend


    !> @brief Implements periodic boundary conditions along x (2 rows of halos on each side)
    !> @param[inout] s Atmospheric state whose halos should be exchanged
  attributes(device) subroutine exchange_halo_x(s)
    implicit none
    class(atmospheric_state), device, intent(inout) :: s
    integer :: i, k, ll
    
    i = blockIdx%x * blockDim%x + threadIdx%x - hs
    k = blockIdx%y * blockDim%y + threadIdx%y
    ll = blockIdx%z * blockDim%z + threadIdx%z
    
    if (k >= 1 .and. k <= nz_loc .and. ll >= 1 .and. ll <= NVARS) then
      s%mem(-1,k,ll)    = s%mem(nx-1,k,ll)
      s%mem(0,k,ll)     = s%mem(nx,k,ll)
      s%mem(nx+1,k,ll) = s%mem(1,k,ll)
      s%mem(nx+2,k,ll) = s%mem(2,k,ll)
    endif
  end subroutine exchange_halo_x

  attributes(device) subroutine exchange_halo_z_device(s, ref)
    implicit none
    type(atmospheric_state), device, intent(inout) :: s
    type(reference_state), device, intent(in) :: ref
    integer :: i, ll
    
    i = blockIdx%x * blockDim%x + threadIdx%x
    
    if (i >= 1-hs .and. i <= nx+hs) then
      if (rank == 0) then
        do ll = 1, NVARS
          if (ll == I_WMOM) then
            s%mem(i,-1,ll) = 0.0_wp
            s%mem(i,0,ll) = 0.0_wp
          else if (ll == I_UMOM) then
            s%mem(i,-1,ll) = s%mem(i,1,ll) / ref%density(1) * ref%density(-1)
            s%mem(i,0,ll) = s%mem(i,1,ll) / ref%density(1) * ref%density(0)
          else
            s%mem(i,-1,ll) = s%mem(i,1,ll)
            s%mem(i,0,ll) = s%mem(i,1,ll)
          end if
        end do
      endif
      
      if (rank == size - 1) then
        do ll = 1, NVARS
          if (ll == I_WMOM) then
            s%mem(i,nz_loc+1,ll) = 0.0_wp
            s%mem(i,nz_loc+2,ll) = 0.0_wp
          else if (ll == I_UMOM) then
            s%mem(i,nz_loc+1,ll) = s%mem(i,nz_loc,ll) / &
                                    ref%density(nz_loc) * ref%density(nz_loc+1)
            s%mem(i,nz_loc+2,ll) = s%mem(i,nz_loc,ll) / &
                                    ref%density(nz_loc) * ref%density(nz_loc+2)
          else
            s%mem(i,nz_loc+1,ll) = s%mem(i,nz_loc,ll)
            s%mem(i,nz_loc+2,ll) = s%mem(i,nz_loc,ll)
          end if
        end do
      endif
    endif
  end subroutine exchange_halo_z_device

  subroutine exchange_halo_z_mpi(s)
    use mpi
    use cudafor
    implicit none
    !type(atmospheric_state_host), intent(inout) :: s_host
    type(atmospheric_state), intent(inout) :: s
    integer :: send_count
    integer :: ierr
    integer :: prev_rank, next_rank
    !integer :: mpi_wp
    !integer :: comm_local
    !integer :: status(MPI_STATUS_SIZE)
    real(wp), allocatable, target, dimension(:) :: send_buffer_down, send_buffer_up
    real(wp), allocatable, target, dimension(:) :: recv_buffer_down, recv_buffer_up
    prev_rank = merge(rank - 1, MPI_PROC_NULL, rank /= 0 )
    next_rank = merge(rank + 1, MPI_PROC_NULL, rank /= size - 1)

    send_count = 2 * (nx + 2 * hs)
  ! Allocate host buffers (contiguous)
  allocate(send_buffer_down(send_count))
  allocate(send_buffer_up(send_count))
  allocate(recv_buffer_down(send_count))
  allocate(recv_buffer_up(send_count))

  do ll = 1, NVARS
    ierr = cudaMemcpy(send_buffer_down, s%mem(1-hs, 1, ll), send_count*8, cudaMemcpyDeviceToHost)
    ierr = cudaMemcpy(send_buffer_up, s%mem(1-hs, nz_loc, ll), send_count*8, cudaMemcpyDeviceToHost)

    ! SENDRECV DOWNWARDS
    call MPI_Sendrecv(send_buffer_down, send_count, MPI_DOUBLE_PRECISION, prev_rank, 0, &
    recv_buffer_down, send_count, MPI_DOUBLE_PRECISION, prev_rank, 0, &
    comm, MPI_STATUS_IGNORE, ierr &
    )

    ! SENDRECV UPWARDS
    call MPI_Sendrecv(send_buffer_up, send_count , MPI_DOUBLE_PRECISION, next_rank, 0, &
    recv_buffer_up, send_count , MPI_DOUBLE_PRECISION, next_rank, 0, &
    comm, MPI_STATUS_IGNORE, ierr &
    )
    ierr = cudaMemcpy(s%mem(1-hs, -1, ll), recv_buffer_down, send_count*8, cudaMemcpyHostToDevice)
    ierr = cudaMemcpy(s%mem(1-hs, -1, ll), recv_buffer_up, send_count*8, cudaMemcpyHostToDevice)

  end do
  deallocate(send_buffer_down, send_buffer_up)
  deallocate(recv_buffer_down, recv_buffer_up)

  end subroutine exchange_halo_z_mpi


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
  end subroutine new_ref

    !> @brief Delete existing reference state
    !> @param[inout] ref Reference state whose memory should be deallocated
  subroutine del_ref(ref)
    implicit none
    class(reference_state), intent(inout) :: ref
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
      ! FIX: Use array assignment instead of host loops for device memory
      flux%mem = xval
    end subroutine set_flux

    !> @brief Deallocate an existing flux object
    !> @param[inout] flux Object which should be deallocated
  subroutine del_flux(flux)
    implicit none
    class(atmospheric_flux), intent(inout) :: flux
    if ( associated(flux%mem) ) then
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

      ! FIX: Use array assignment instead of host loops for device memory
      tend%mem = xval
    end subroutine set_tendency

    !> @brief Deallocate an existing tendency object
    !> @param[inout] tend Object which should be deallocated
  subroutine del_tendency(tend)
    implicit none
    class(atmospheric_tendency), intent(inout) :: tend
    if ( associated(tend%mem) ) then
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
      
      ! FIX for NVFORTRAN-S-0519:
      ! Replaced explicit loops with array assignment.
      ! CUDA Fortran handles x%mem (device) = y%mem (device) automatically and efficiently.
      if (associated(x%mem) .and. associated(y%mem)) then
         x%mem = y%mem
      end if
      
    end subroutine state_equal_to_state
#endif
end module module_types_cuda