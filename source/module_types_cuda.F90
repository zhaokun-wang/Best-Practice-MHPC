module module_types_cuda
#ifdef _CUDA_KERN
  use physical_constants
  use physical_parameters
  use parallel_parameters
  use indexing
  use dimensions
  use mpi
  use cudafor

  implicit none

  private

  public :: reference_state
  public :: atmospheric_state
  public :: atmospheric_flux
  public :: atmospheric_tendency
  
  public :: new_state, del_state, new_ref, del_ref, new_flux, del_flux, new_tendency, del_tendency
  public :: xtend_kernel, ztend_kernel, exchange_halo_x_kernel, exchange_halo_z_host
  public :: assignment(=)

  !> @brief Reference state (SOLO DATI)
  type reference_state
    real(wp), device, allocatable, dimension(:) :: density
    real(wp), device, allocatable, dimension(:) :: denstheta
    real(wp), device, allocatable, dimension(:) :: idens
    real(wp), device, allocatable, dimension(:) :: idenstheta
    real(wp), device, allocatable, dimension(:) :: pressure
  end type reference_state

  !> @brief Atmospheric state (SOLO DATI)
  type atmospheric_state
    real(wp), device, pointer, dimension(:,:,:) :: mem => null()
    real(wp), device, pointer, dimension(:,:) :: dens
    real(wp), device, pointer, dimension(:,:) :: umom
    real(wp), device, pointer, dimension(:,:) :: wmom
    real(wp), device, pointer, dimension(:,:) :: rhot
  end type atmospheric_state

  !> @brief Flux (SOLO DATI)
  type atmospheric_flux
    real(wp), device, pointer, dimension(:,:,:) :: mem => null()
    real(wp), device, pointer, dimension(:,:) :: dens
    real(wp), device, pointer, dimension(:,:) :: umom
    real(wp), device, pointer, dimension(:,:) :: wmom
    real(wp), device, pointer, dimension(:,:) :: rhot
  end type atmospheric_flux

  !> @brief Tendency (SOLO DATI)
  type atmospheric_tendency
    real(wp), device, pointer, dimension(:,:,:) :: mem => null()
    real(wp), device, pointer, dimension(:,:) :: dens
    real(wp), device, pointer, dimension(:,:) :: umom
    real(wp), device, pointer, dimension(:,:) :: wmom
    real(wp), device, pointer, dimension(:,:) :: rhot
  end type atmospheric_tendency

  interface assignment(=)
    module procedure state_equal_to_state
  end interface assignment(=)

  contains

! ==============================================================================
! SEZIONE HOST (Gestione Memoria e Orchestrazione)
! ==============================================================================

  subroutine new_state(atmo)
    implicit none
    type(atmospheric_state), intent(inout) :: atmo
    if ( associated(atmo%mem) ) deallocate(atmo%mem)
    allocate(atmo%mem(1-hs:nx+hs, 1-hs:nz_loc+hs, NVARS))
    atmo%dens(1-hs:,1-hs:) => atmo%mem(:,:,I_DENS)
    atmo%umom(1-hs:,1-hs:) => atmo%mem(:,:,I_UMOM)
    atmo%wmom(1-hs:,1-hs:) => atmo%mem(:,:,I_WMOM)
    atmo%rhot(1-hs:,1-hs:) => atmo%mem(:,:,I_RHOT)
  end subroutine new_state

  subroutine del_state(atmo)
    implicit none
    type(atmospheric_state), intent(inout) :: atmo
    if ( associated(atmo%mem) ) deallocate(atmo%mem)
    nullify(atmo%dens, atmo%umom, atmo%wmom, atmo%rhot)
  end subroutine del_state

  subroutine new_ref(ref)
     implicit none
     type(reference_state), intent(inout) :: ref
     allocate(ref%density(1-hs:nz_loc+hs))
     allocate(ref%denstheta(1-hs:nz_loc+hs))
     allocate(ref%idens(nz_loc+1))
     allocate(ref%idenstheta(nz_loc+1))
     allocate(ref%pressure(nz_loc+1))
  end subroutine new_ref

  subroutine del_ref(ref)
     implicit none
     type(reference_state), intent(inout) :: ref
     if(allocated(ref%density)) deallocate(ref%density)
     if(allocated(ref%denstheta)) deallocate(ref%denstheta)
     if(allocated(ref%idens)) deallocate(ref%idens)
     if(allocated(ref%idenstheta)) deallocate(ref%idenstheta)
     if(allocated(ref%pressure)) deallocate(ref%pressure)
  end subroutine del_ref

  subroutine new_flux(flux)
    implicit none
    type(atmospheric_flux), intent(inout) :: flux
    if ( associated(flux%mem) ) deallocate(flux%mem)
    allocate(flux%mem(1:nx+1, 1:nz_loc+1,NVARS))
    flux%dens => flux%mem(:,:,I_DENS)
    flux%umom => flux%mem(:,:,I_UMOM)
    flux%wmom => flux%mem(:,:,I_WMOM)
    flux%rhot => flux%mem(:,:,I_RHOT)
  end subroutine new_flux

  subroutine del_flux(flux)
    implicit none
    type(atmospheric_flux), intent(inout) :: flux
    if ( associated(flux%mem) ) deallocate(flux%mem)
    nullify(flux%dens, flux%umom, flux%wmom, flux%rhot)
  end subroutine del_flux

  subroutine new_tendency(tend)
    implicit none
    type(atmospheric_tendency), intent(inout) :: tend
    if ( associated(tend%mem) ) deallocate(tend%mem)
    allocate(tend%mem(nx, nz_loc,NVARS))
    tend%dens => tend%mem(:,:,I_DENS)
    tend%umom => tend%mem(:,:,I_UMOM)
    tend%wmom => tend%mem(:,:,I_WMOM)
    tend%rhot => tend%mem(:,:,I_RHOT)
  end subroutine new_tendency

  subroutine del_tendency(tend)
    implicit none
    type(atmospheric_tendency), intent(inout) :: tend
    if ( associated(tend%mem) ) deallocate(tend%mem)
    nullify(tend%dens, tend%umom, tend%wmom, tend%rhot)
  end subroutine del_tendency


  !> @brief Scambio Halo Z gestito su HOST (Risolve errore MPI Generic Procedure)
  subroutine exchange_halo_z_host(s, ref, comm_id)
    implicit none
    type(atmospheric_state), intent(inout) :: s
    type(reference_state), intent(in) :: ref
    integer, intent(in) :: comm_id
    
    integer :: ierr, send_count, ll
    integer(8) :: t_comm_start, t_comm_end, rate
    
    ! Buffer HOST temporanei per MPI
    ! Assumiamo profondità halo = 2 (hs=2), quindi scambiamo 2 righe
    real(wp), allocatable, dimension(:,:) :: send_buf_down, recv_buf_down
    real(wp), allocatable, dimension(:,:) :: send_buf_up, recv_buf_up

    send_count = 2 * (nx + 2 * hs)
    
    ! Allochiamo buffer HOST della dimensione corretta
    ! Dim 1: Lunghezza riga X (inclusi halo X) = 1-hs : nx+hs
    ! Dim 2: Profondità halo Z = 2
    allocate(send_buf_down(1-hs:nx+hs, 2))
    allocate(recv_buf_down(1-hs:nx+hs, 2))
    allocate(send_buf_up(1-hs:nx+hs, 2))
    allocate(recv_buf_up(1-hs:nx+hs, 2))

    call system_clock(t_comm_start)

    ! Loop su tutte le variabili (NVARS)
    ! Necessario perché s%mem ha 3 dimensioni (x, z, var)
    do ll = 1, NVARS
        
        ! 1. DEVICE -> HOST (Copia implicita CUDA Fortran)
        ! Prendi le prime 2 righe interne (k=1,2) per mandarle GIÙ
        send_buf_down = s%mem(1-hs:nx+hs, 1:2, ll)
        
        ! Prendi le ultime 2 righe interne (k=nz_loc-1, nz_loc) per mandarle SU
        send_buf_up = s%mem(1-hs:nx+hs, nz_loc-1:nz_loc, ll)

        ! 2. MPI COMMUNICATION (Host to Host)
        ! Scambio verso il basso
        call MPI_Sendrecv(send_buf_down, send_count, MPI_DOUBLE_PRECISION, prev_rank, 0, &
                          recv_buf_down, send_count, MPI_DOUBLE_PRECISION, prev_rank, 0, &
                          comm_id, MPI_STATUS_IGNORE, ierr)

        ! Scambio verso l'alto
        call MPI_Sendrecv(send_buf_up, send_count, MPI_DOUBLE_PRECISION, next_rank, 0, &
                          recv_buf_up, send_count, MPI_DOUBLE_PRECISION, next_rank, 0, &
                          comm_id, MPI_STATUS_IGNORE, ierr)
        
        ! 3. HOST -> DEVICE (Copia implicita)
        ! I dati ricevuti dal basso vanno nell'halo inferiore (k=-1,0)
        ! Nota: recv_buf_down contiene k=-1 (riga 1) e k=0 (riga 2)
        s%mem(1-hs:nx+hs, -1:0, ll) = recv_buf_down
        
        ! I dati ricevuti dall'alto vanno nell'halo superiore (k=nz_loc+1, nz_loc+2)
        s%mem(1-hs:nx+hs, nz_loc+1:nz_loc+2, ll) = recv_buf_up

    end do
    
    call system_clock(t_comm_end, rate)
    ! T_communicate = T_communicate + dble(t_comm_end-t_comm_start)/dble(rate)

    deallocate(send_buf_down, recv_buf_down, send_buf_up, recv_buf_up)

    ! 4. Applicazione Boundary Conditions fisiche (Kernel GPU)
    ! Lanciamo il kernel per fissare le velocità a 0 o scalarle sui bordi
    call update_halo_z_kernel<<<dim3((nx+127)/128, 1, 1), dim3(128, 1, 1)>>>(s, ref)

  end subroutine exchange_halo_z_host


! ==============================================================================
! SEZIONE DEVICE (Kernels)
! ==============================================================================

  !> Kernel per calcolare XTEND.
  attributes(global) subroutine xtend_kernel(tendency, flux, ref, atmostat, dx, dt)
    implicit none
    type(atmospheric_tendency), device, intent(inout) :: tendency
    type(atmospheric_flux), device, intent(inout) :: flux
    type(reference_state), device, intent(in) :: ref
    type(atmospheric_state), device, intent(inout) :: atmostat
    real(wp), value :: dx, dt
    
    integer :: i, k, ll, s
    real(wp) :: r, u, w, t, p, hv_coef
    real(wp), dimension(STEN_SIZE) :: stencil
    real(wp), dimension(NVARS) :: d3_vals, vals

    i = blockIdx%x * blockDim%x + threadIdx%x
    k = blockIdx%y * blockDim%y + threadIdx%y
    ll = blockIdx%z * blockDim%z + threadIdx%z 

    hv_coef = -hv_beta * dx / (16.0_wp*dt)

    if ( i >= 1 .and. i <= nx .and. k >= 1 .and. k <= nz_loc .and. ll >= 1 .and. ll <= NVARS ) then
      do s = 1, STEN_SIZE
        stencil(s) = atmostat%mem(i-hs-1+s,k,ll)
      end do
      vals(ll) = - 1.0_wp * stencil(1)/12.0_wp + 7.0_wp * stencil(2)/12.0_wp + &
                 7.0_wp * stencil(3)/12.0_wp - 1.0_wp * stencil(4)/12.0_wp
      d3_vals(ll) = - 1.0_wp * stencil(1) + 3.0_wp * stencil(2) - &
                    3.0_wp * stencil(3) + 1.0_wp * stencil(4)

      r = vals(I_DENS) + ref%density(k)
      u = vals(I_UMOM) / r
      w = vals(I_WMOM) / r
      t = ( vals(I_RHOT) + ref%denstheta(k) ) / r
      p = c0*(r*t)**cdocv

      flux%dens(i,k) = r*u - hv_coef*d3_vals(I_DENS)
      flux%umom(i,k) = r*u*u+p - hv_coef*d3_vals(I_UMOM)
      flux%wmom(i,k) = r*u*w - hv_coef*d3_vals(I_WMOM)
      flux%rhot(i,k) = r*u*t - hv_coef*d3_vals(I_RHOT)
      
      if (i <= nx) then
         tendency%mem(i,k,ll) = -( flux%mem(i+1,k,ll) - flux%mem(i,k,ll) ) / dx
      end if
    end if
  end subroutine xtend_kernel

  attributes(global) subroutine ztend_kernel(tendency, flux, ref, atmostat, dz, dt)
      implicit none
      type(atmospheric_tendency), device, intent(inout) :: tendency
      type(atmospheric_flux), device, intent(inout) :: flux
      type(reference_state), device, intent(in) :: ref
      type(atmospheric_state), device, intent(inout) :: atmostat
      real(wp), value :: dz, dt
      integer :: i, k, ll, s, j, k_flux
      real(wp) :: r, u, w, t, p, hv_coef
      real(wp), dimension(STEN_SIZE) :: stencil
      real(wp), dimension(NVARS) :: d3_vals, vals
  
      i = blockIdx%x * blockDim%x + threadIdx%x
      k = blockIdx%y * blockDim%y + threadIdx%y
      ll = blockIdx%z * blockDim%z + threadIdx%z 
  
      hv_coef = -hv_beta * dz / (16.0_wp*dt)
      if ( i >= 1 .and. i <= nx .and. k >= 1 .and. k <= nz_loc .and. ll >= 1 .and. ll <= NVARS) then
        do s = 1, STEN_SIZE
          stencil(s) = atmostat%mem(i,k-hs-1+s,ll)
        end do
        vals(ll) = - 1.0_wp * stencil(1)/12.0_wp + 7.0_wp * stencil(2)/12.0_wp + &
                   7.0_wp * stencil(3)/12.0_wp - 1.0_wp * stencil(4)/12.0_wp
        d3_vals(ll) = - 1.0_wp * stencil(1) + 3.0_wp * stencil(2) - &
                      3.0_wp * stencil(3) + 1.0_wp * stencil(4)
        r = vals(I_DENS) + ref%idens(k)
        u = vals(I_UMOM) / r
        w = vals(I_WMOM) / r
        t = ( vals(I_RHOT) + ref%idenstheta(k) ) / r
        p = c0*(r*t)**cdocv - ref%pressure(k)
        
        if ((k == 1 .and. rank == 0) .or. (k == nz_loc+1 .and. rank == size - 1)) then
          w = 0.0_wp
          d3_vals(I_DENS) = 0.0_wp
        end if
        
        flux%dens(i,k) = r*w - hv_coef*d3_vals(I_DENS)
        flux%umom(i,k) = r*w*u - hv_coef*d3_vals(I_UMOM)
        flux%wmom(i,k) = r*w*w+p - hv_coef*d3_vals(I_WMOM)
        flux%rhot(i,k) = r*w*t - hv_coef*d3_vals(I_RHOT)
  
        ! Tendency calculation (Simplified for GPU context - avoiding inner loop dependencies if possible)
        ! Assuming flux is computed, we calculate tendency for THIS k
        if ( k <= nz_loc ) then
          tendency%mem(i,k,ll) = -( flux%mem(i,k+1,ll) - flux%mem(i,k,ll) ) / dz
          if (ll == I_WMOM) then
              tendency%wmom(i,k) = tendency%wmom(i,k) - atmostat%dens(i,k)*grav
          end if
        end if
      end if
  end subroutine ztend_kernel

  !> Kernel per aggiornare i bordi Z (da chiamare DOPO MPI)
  attributes(global) subroutine update_halo_z_kernel(s, ref)
    implicit none
    type(atmospheric_state), device, intent(inout) :: s
    type(reference_state), device, intent(in) :: ref
    integer :: i, ll
    
    ! Grid: X dimension (blockIdx%x), Y dimension = 1 (usiamo loop interno per ll)
    i = blockIdx%x * blockDim%x + threadIdx%x - hs
    
    if (i >= 1-hs .and. i <= nx+hs) then
      do ll = 1, NVARS
            if (rank == 0) then
                if (ll == I_WMOM) then
                    s%mem(i,-1,ll) = 0.0_wp
                    s%mem(i,0,ll) = 0.0_wp
                else if (ll == I_UMOM) then
                    s%mem(i,-1,ll) = s%mem(i,1,ll) / ref%density(1) * ref%density(-1)
                    s%mem(i,0,ll)  = s%mem(i,1,ll) / ref%density(1) * ref%density(0)
                else
                    s%mem(i,-1,ll) = s%mem(i,1,ll)
                    s%mem(i,0,ll)  = s%mem(i,1,ll)
                end if
            end if
            
            if (rank == size - 1) then
                 if (ll == I_WMOM) then
                    s%mem(i,nz_loc+1,ll) = 0.0_wp
                    s%mem(i,nz_loc+2,ll) = 0.0_wp
                  else if (ll == I_UMOM) then
                    s%mem(i,nz_loc+1,ll) = s%mem(i,nz_loc,ll) / ref%density(nz_loc) * ref%density(nz_loc+1)
                    s%mem(i,nz_loc+2,ll) = s%mem(i,nz_loc,ll) / ref%density(nz_loc) * ref%density(nz_loc+2)
                  else
                    s%mem(i,nz_loc+1,ll) = s%mem(i,nz_loc,ll)
                    s%mem(i,nz_loc+2,ll) = s%mem(i,nz_loc,ll)
                  end if
            end if
      end do
    end if
  end subroutine update_halo_z_kernel

  attributes(global) subroutine exchange_halo_x_kernel(s)
    implicit none
    type(atmospheric_state), device, intent(inout) :: s
    integer :: k, ll
    
    k = blockIdx%x * blockDim%x + threadIdx%x
    ll = blockIdx%y * blockDim%y + threadIdx%y
    
    if (k >= 1 .and. k <= nz_loc .and. ll >= 1 .and. ll <= NVARS) then
       s%mem(-1,k,ll)   = s%mem(nx-1,k,ll)
       s%mem(0,k,ll)    = s%mem(nx,k,ll)
       s%mem(nx+1,k,ll) = s%mem(1,k,ll)
       s%mem(nx+2,k,ll) = s%mem(2,k,ll)
    end if
  end subroutine exchange_halo_x_kernel

  subroutine state_equal_to_state(x,y)
      implicit none
      type(atmospheric_state), intent(inout) :: x
      type(atmospheric_state), intent(in) :: y
      !$acc parallel loop collapse(3) present(x%mem, y%mem)
      do ll = 1, NVARS
        do k = 1-hs, nz_loc+hs
          do i = 1-hs, nx+hs
            x%mem(i,k,ll) = y%mem(i,k,ll)
          end do
        end do
      end do
  end subroutine state_equal_to_state

#endif
end module module_types_cuda
