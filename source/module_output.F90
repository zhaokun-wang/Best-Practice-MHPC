!>
!! @authors Pedde, Zhaokun, Veraldi
!! @date 1-12-2025
!! @brief model_output.f90 is the module containing all related to I/O of the simulation

module module_output
  use calculation_types, only : wp, iowp
  use parallel_parameters, only : i_beg, k_beg, nz_loc, comm, rank, size, rest, base !**** PARALLEL ****
  use dimensions, only : nx, nz
  use module_types, only : atmospheric_state, reference_state
  use iodir, only : stderr
  use netcdf
  use mpi

  implicit none

  private

  public :: create_output
  public :: write_record
  public :: close_output

  real(wp), allocatable, dimension(:,:) :: dens     !>The density grid
  real(wp), allocatable, dimension(:,:) :: uwnd     !>the velocity on x grid
  real(wp), allocatable, dimension(:,:) :: wwnd     !>the velocity on z grid
  real(wp), allocatable, dimension(:,:) :: theta    !>the temperature potential grid

  integer :: ncid                                                       !>the id of the file
  integer :: dens_varid, uwnd_varid, wwnd_varid, theta_varid, t_varid   !>the id of the variables
  integer :: rec_out

  contains

    !>
    !!@brief creation of the output system with netcdf
    subroutine create_output
      implicit none

      !the dimension id fot the grid
      integer :: t_dimid, x_dimid, z_dimid
      integer :: status_nc
      logical :: is_parallel

      !allocate the arrayy of the 4 variables
      !**** PARALLEL **** array now have the local dimension for the MPI rank
      allocate(dens(nx,nz_loc))
      allocate(uwnd(nx,nz_loc))
      allocate(wwnd(nx,nz_loc))
      allocate(theta(nx,nz_loc))

      !create the gile output.nc, overwrite it (_clobber) if exist, from now id file is in ncid
      !! call ncwrap(nf90_create('output.nc',nf90_clobber,ncid), __LINE__)
      !**** PARALLEL ****
      !create the file to be write in parallel

      call ncwrap(nf90_create_par('output.nc', ior(nf90_clobber, nf90_netcdf4), comm, MPI_INFO_NULL, ncid), __LINE__)

      !definitions of dimensions of the grid per step, ATTENTION: time _unlimited because variable in sim
      !nf90_def_dim(id_file, <label>, dimension of grid for that axes, id_grid_direction)
      call ncwrap(nf90_def_dim(ncid,'time',nf90_unlimited,t_dimid), __LINE__)
      call ncwrap(nf90_def_dim(ncid,'x',nx,x_dimid), __LINE__)
      call ncwrap(nf90_def_dim(ncid,'z',nz,z_dimid), __LINE__)

      !definitions of variables in the grid (prepare to receive that variable)
      !nf90_def_var(id_file, <label>, type, id_time in which put var, id_var)
      call ncwrap(nf90_def_var(ncid,'time',iowp,[t_dimid],t_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'rho',iowp, &
          [x_dimid,z_dimid,t_dimid],dens_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'u',iowp, &
          [x_dimid,z_dimid,t_dimid],uwnd_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'w',iowp, &
          [x_dimid,z_dimid,t_dimid],wwnd_varid), __LINE__)
      call ncwrap(nf90_def_var(ncid,'theta',iowp, &
          [x_dimid,z_dimid,t_dimid],theta_varid), __LINE__)

      !close definitions area, set rec_out for the number of printing file
      call ncwrap(nf90_enddef(ncid), __LINE__)

      !*** PARALLEL ***
      !say to netCDF that all the variables have a collective acces, not indipendent
      !to write in parallel in all EXCEPT TIME
      call ncwrap(nf90_var_par_access(ncid, t_varid, nf90_collective), __LINE__)
      call ncwrap(nf90_var_par_access(ncid, dens_varid, nf90_collective), __LINE__)
      call ncwrap(nf90_var_par_access(ncid, uwnd_varid, nf90_collective), __LINE__)
      call ncwrap(nf90_var_par_access(ncid, wwnd_varid, nf90_collective), __LINE__)
      call ncwrap(nf90_var_par_access(ncid, theta_varid, nf90_collective), __LINE__)

      rec_out = 1
    end subroutine create_output

    !>
    !! @brief write an entry in the netCDF file
    !! @param[in] atmostat
    !! @param[in] ref
    !! @param[in] etime (actual time in the simulation that we are going to print)
    subroutine write_record(atmostat,ref,etime)
      implicit none
      type(atmospheric_state), intent(in) :: atmostat
      type(reference_state), intent(in) :: ref
      real(wp), intent(in) :: etime
      integer :: i, k
      integer, dimension(1) :: st1, ct1
      integer, dimension(3) :: st3, ct3
      real(wp), dimension(1) :: etimearr
      integer(8) :: t1, t2, rate

      !put the variables from atmostat in the right array of I/O
      !*** PARALLEL ***
      !now we fill the variables with nz_loc as dimension in z


      !$acc parallel loop collapse(2) copyin(atmostat%dens, atmostat%umom, atmostat%wmom, atmostat%rhot, ref%density, ref%denstheta) &
      !$acc copy(dens) copyout(uwnd, wwnd, theta)
      !$omp parallel do collapse(2) private(k,i)
      do k = 1, nz_loc
        do i = 1, nx
          dens(i,k) = atmostat%dens(i,k)
          uwnd(i,k) = atmostat%umom(i,k)/(ref%density(k)+dens(i,k))
          wwnd(i,k) = atmostat%wmom(i,k)/(ref%density(k)+dens(i,k))
          theta(i,k) = (atmostat%rhot(i,k) + ref%denstheta(k)) / &
              (ref%density(k) + dens(i,k)) - ref%denstheta(k)/ref%density(k)
        end do
      end do
      !$omp end parallel do
      !$acc end parallel loop



      !Writing part
      !*** PARALLEL ***
      !k_beg now I suppose is the point where is global starting of rank
      st3 = [ i_beg, k_beg, rec_out ]   !>cursor coordinate where starting writing
      ct3 = [ nx, nz_loc, 1 ]               !>define the dimension of the block where to write


      call ncwrap(nf90_put_var(ncid,dens_varid,dens,st3,ct3), __LINE__)
      call ncwrap(nf90_put_var(ncid,uwnd_varid,uwnd,st3,ct3), __LINE__)
      call ncwrap(nf90_put_var(ncid,wwnd_varid,wwnd,st3,ct3), __LINE__)
      call ncwrap(nf90_put_var(ncid,theta_varid,theta,st3,ct3), __LINE__)


      !writing the time
      !*** PARALLEL ***
      !only the rank 0 can do it (need only one time writing the time)

      st1 = [ rec_out ]
      ct1 = [ 1 ]
      etimearr(1) = etime
      call ncwrap(nf90_put_var(ncid, t_varid, etimearr, start=st1, count=ct1), __LINE__)

      !update the receive file number
      rec_out = rec_out + 1

    end subroutine write_record

    !>
    !! @brief close the otuput file and deallocate the array for the output
    subroutine close_output
      implicit none
      if ( allocated(dens) ) then
        deallocate(dens)
        deallocate(uwnd)
        deallocate(wwnd)
        deallocate(theta)
      end if
      call ncwrap(nf90_close(ncid), __LINE__)
    end subroutine close_output

    !>
    !! @brief check if there is an error in NetCDF I/O procedure
    subroutine ncwrap(ierr,line)
      implicit none
      integer, intent(in) :: ierr
      integer, intent(in) :: line
      if (ierr /= nf90_noerr) then
        write(stderr,*) 'NetCDF Error at line: ', line
        write(stderr,*) nf90_strerror(ierr)
        stop
      end if
    end subroutine ncwrap

end module module_output
