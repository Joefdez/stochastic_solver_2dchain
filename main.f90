program main

  ! This is a skeleton code for the numerical solution of systems of N stochastic differential equations
  ! The basic and essential variables and tasks are included. The specalization to a specific problem should be easy following the instructions
  ! The parallelization is achieved using MPI and is isolated to the main.f90 file. Therefore all that is requiered from the user is the implementation
  !       of the equations to be solved.
  ! The solver included follows the Platen 2.0 weak scheme for the case in which there are no stocahstic couplings between the different unknown functions
  !       ( no off-diagonal terms in the B matrix containing the stochastic terms)

  ! This particular file is the main body of the code. It should be reduced as much as possible to simply a collection of function calls.

  ! Created by Joseph John Fernández, ULL (2017)

  use mpi
  use constants
  use initialization
  implicit none

  include "declarations.f90"

  integer :: rank, procs, status(MPI_STATUS_SIZE), source, ierr                 ! MPI variables
  call MPI_INIT(ierr)                                                           ! Neccesary MPI initialization calls
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)                                ! rank -> identifies processes
  call MPI_COMM_SIZE(MPI_COMM_WORLD, procs, ierr)                               ! procs -> number of processes

  ! Initialize the solution parameters. Other parameters can be added if needed for a specific application.
  ! The source code is included in "init/initialization.f90".
  ! The data is read from the file "data/solve_values.f90" by process 0 and distributed to the other processes.
  if(rank .eq. 0) then
    print*, "Process 0 Reading solution parameters"
    call initialize_system(neqs, tt, dt, traj, save_freq, fin, share_freq)
    print*, "Sharing data with other processes"
  end if

! After initializing the solution data, broadcast to the rest of processes following the example below:
  ! (repeat as many lines as neccesary, be careful to use the correct parameter in sned_type!!!!)
  ! call mpi_bcast(data_to_be_sent, send_count, send_type, broadcasting_process_ID, comm, ierr)

  call mpi_bcast(neqs, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(tt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(dt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(save_freq, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(fin, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(share_freq, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)

  nsteps        =  int(tt/dt)                                                   ! Number of steps in numerical integration
  nsample_steps =  int(nsteps/save_freq)                                        ! Number of steps to be saved (sample)
  final_steps   =  int(fin*nsteps)                                              ! Final number of steps for stead-state calculations
  half          =  neqs/2
  quarter       =  neqs/4

  ! Now that the basic parameters controlling the solution of the problem are loaded, the vector and matrix allocations can be performed
  include "allocations.f90"

  ! Distribute the stochastic realizations between the different processes
  local_traj = traj/procs
  rem = mod(traj, procs)
  if (rank .lt. rem) local_traj = local_traj + 1

! Now it is time to call the initialization call specific to the problem. This routine should be written for the specific
  ! application of the code and stored in "init/initialization.f90". Again, it is read by proc. 0 and distributed subsequently.
  ! Fill in as needed.
  if(rank .eq. 0) then
    print*, "Process 0 reading problem-specific parameters"
    call initialize_problem(nbath, nparticles, dim, mass, charge, omega0, omega1, omega2, Gam, delL, delR, IL, IR)
    print*, "Sharing data with other processes"
  end if
  ! A broadcast call for each value returned in initialize_system
  ! call mpi_bcast(data_to_be_sent, send_count, send_type, broadcasting_process_ID, comm, ierr)
  call mpi_bcast(nbath, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(nparticles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(dim, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(mass, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(charge, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(omega, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Gam, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(delL, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(delR, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(IL, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(IR, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  ! Intermediate calculations
  omegaL = delL + omega0
  omegaR = delR + omega0
  kL =  omegaL / cc   ! Not really wavelengths, dimensionally they are wavenumbers
  kR =  omegaR / cc
 ! Calculate diffusion and friction coefficients
  call doppler_values(kL, Gam, delL, IL, etaL, DL)
  call doppler_values(kR, Gam, delR, IR, etaR, DR)
  ! Calculate dimensionless Doppler cooling parameters
  call dimensionless_doppler_values(etaL, DL, mass, long_freq, char_length, aetaL, aDL)
  call dimensionless_doppler_values(etaR, DR, mass, long_freq, char_length, aetaR, aDR)
  ! Initialization of the stochastic terms of the equations. This is specific to the particular application.
  ! Several standard skeletons are included. See "functions/stochastic.f90"

  !call stochastic inizializer
  BB((2*nparticles+1):(2*nparticles+bath)) =
  BB(3)
  BB() =
  BB() =

  !call initial conditions functions
  if(rank .eq. 0) then
    print*, "Process 0 reading problem-specific parameters"
    call initial_conditions()
    print*, "Sharing data with other processes"
  end if
  ! A broadcast call for each value returned in initial_conditions. As many as neccesary
  ! call mpi_bcast(data_to_be_sent, send_count, send_type, broadcasting_process_ID, comm, ierr)
  call mpi_bcast(,0, MPI_COMM_WORLD, ierr)

  ! The following loop executes the calculation of the stochastic realizaitions
  do kk=1, local_traj, 1
    print*, "Proc ", rank, " on stoch. realization ", kk
    ! Call the time solver
    call platen_realization(YY0, YY, YYf, BB, nbath, dt, dst, nsteps, nsample_steps, final_steps)
    PP2 = PP2 + 0.5d0*(YY(:,half:(half+quarter))*YY(:,half:(half+quarter)) + YY(:,(half+quarter+1):)*YY(:,(half+quarter+1):))
    if( ( mod(kk,share_freq) .eq. 0) .and. (kk .lt. local_traj) ) then
      ! Partial sharing of results to proc. 0 and calculate intermediate averages. Use reduction functions.
      !call MPI_REDUCE(data_to_be_reduced, recieve_buffer, count, type, op, recieve_id, comm, ierr)
      call mpi_reduce(kk, partial_t, 1, mpi_integer, sum, 0, MPI_COMM_WORLD, ierr)
      call mpi_reduce(PP2, PP2_av, nelems, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
      if(rank .eq. 0) then
        ! calculate the averages
        PP2_av = PP2_av/partial_t
        ! save the results
        open(unit=12,file='PP2')
        do jj=1, nparticles, 1
          write(12,*) PP2_av(jj,:)
        end do
        close(unit=12)
        PP2_av    = 0.0d0
        partial_t = 0
      end if
    end if
  end do

  ! Final sharing of results to proc. 0 and calculate averages. Use reduction functions.
  !call MPI_REDUCE(data_to_be_reduced, recieve_buffer, count, type, op, recieve_id, comm, ierr)
  call mpi_reduce(kk, traj, 1, mpi_integer, sum, 0, MPI_COMM_WORLD, ierr)
  call mpi_reduce(PP2, PP2_av, nelems, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  if(rank .eq. 0) then
    ! calculate the averages
    PP2_av = PP2_av/traj
    ! save the results
    open(unit=12,file='PP2')
    do jj=1, nparticles, 1
      write(12,*) PP2_av(jj,:)
    end do
    close(unit=12)
  end if


  call mpi_finalize(ierr)                                                       ! Finalize the MPI related processes for correct end-of-run time
end program main
