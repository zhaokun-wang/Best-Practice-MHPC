!>
!! @authors Pedde, Zhaokun, Veraldi
!! @date 1-12-2025
!! @brief model.f90 is the main of the flluid dynamics simulation

program atmosphere_model
  !**** Module load area ****
  use calculation_types, only : wp
  use module_physics, only : dt, oldstat, newstat, flux, tend, ref
  use module_physics, only : init, finalize
  use module_physics, only : rungekutta, total_mass_energy
  use module_output, only : create_output, write_record, close_output
  use dimensions , only : sim_time, output_freq
  use iodir, only : stdout
  use mpi
  use parallel_parameters
  implicit none

  !**** Variables declaration area ****
  real(wp) :: etime
  real(wp) :: ptime
  real(wp) :: output_counter
  real(wp) :: pctime
  real(wp) :: mass0, te0
  real(wp) :: mass1, te1
  integer(8) :: t1, t2, rate

  !Parallel init
  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  !Rank and Size
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, size, ierr)

  !Prev and Next ranks
  prev_rank = merge(rank - 1, MPI_PROC_NULL, rank /= 0 )
  next_rank = merge(rank + 1, MPI_PROC_NULL, rank /= size - 1)

  !**** Initialization region ****
  write(stdout, *) 'SIMPLE ATMOSPHERIC MODEL STARTING.'
  call init(etime,output_counter,dt)                    !>initialize old state and new state
  call total_mass_energy(mass0,te0)                     !>initalize mass and temperature at start
  call create_output( )                                 !>create the .nc for the output storing
  call write_record(oldstat,ref,etime)                  !>write the first record

  !*** timing ***
  call system_clock(t1)

  !**************************** SIMULATION CYCLE BLOCK ***************************
  ptime = int(sim_time/10.0)

  do while (etime < sim_time)
    !check case in which the last step to do to end is smaller thend set dt
    if (etime + dt > sim_time) dt = sim_time - etime

    !**************** EVOLUTION *****************
    call rungekutta(oldstat,newstat,flux,tend,dt)
    !********************************************

    !execution percentage write
    if ( mod(etime,ptime) < dt ) then
      if (rank == 0) then
        pctime = (etime/sim_time)*100.0_wp
        write(stdout,'(1x,a,i2,a)') 'TIME PERCENT : ', int(pctime), '%'
      end if
    end if

    !updating the actual time and the otput counter
    etime = etime + dt
    output_counter = output_counter + dt

    !printing area
    if (output_counter >= output_freq) then
      output_counter = output_counter - output_freq
      call write_record(oldstat,ref,etime)
    end if
  end do

  !******************************************************************************

  !**** final printing for checking and timings results ****
  call total_mass_energy(mass1,te1)
  call close_output( )

  if (rank == 0) then
    write(stdout,*) "----------------- Atmosphere check ----------------"
    write(stdout,*) "Fractional Delta Mass  : ", (mass1-mass0)/mass0
    write(stdout,*) "Fractional Delta Energy: ", (te1-te0)/te0
    write(stdout,*) "---------------------------------------------------"
  end if
  call finalize()
  call system_clock(t2,rate)

  if (rank == 0) then
    write(stdout,*) "SIMPLE ATMOSPHERIC MODEL RUN COMPLETED."
    write(stdout,*) "USED CPU TIME: ", dble(t2-t1)/dble(rate)
  endif

  call MPI_Finalize(ierr)

end program atmosphere_model
