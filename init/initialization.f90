module initialization
  use constants
  implicit none

contains
  subroutine initialize_system(neqs, tt, dt, traj, save_freq, fin, share_freq)
    ! based on http://jblevins.org/log/control-file
    implicit none

    integer, intent(inout)            :: neqs, traj, save_freq, share_freq
    real(kind=8), intent(inout)       :: dt, tt, fin

    character(len=100) :: buffer, local
    character(len=100) :: label
    integer :: pos
    integer, parameter :: fh=15
    integer :: ios=0, line=0


    open(unit=fh, file='values.dat',action='read')

    do while(ios==0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios==0) then
        line = line + 1

        ! find the first instance of whitespace. Split and label data
        pos=scan(buffer,' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)
        select case(label)
        case('neqs')
              read(buffer,*,iostat=ios) neqs
          case('dt')
              read(buffer,*,iostat=ios) dt
          case('tt')
              read(buffer,*,iostat=ios) tt
          case('traj')
              read(buffer,*,iostat=ios) traj
          case('save_freq')
              read(buffer,*,iostat=ios) traj
          case('fin')
              read(buffer,*,iostat=ios) traj
          case('share_freq')
              read(buffer,*,iostat=ios) share_freq
          case default
              print*, "Skipping invalid value."
        end select
      end if
    end do

    close(unit=fh)

  end subroutine initialize_system

  subroutine initialize_problem(nbath, nparticles, dim, mass, charge, omega0, Gam, delL, delR, IL, IR)
    ! based on http://jblevins.org/log/control-file
    implicit none

    integer, intent(inout)            :: nbath, nparticles, dim
    real(kind=8), intent(inout)       :: mass, charge
    real(kind=8), intent(inout)       :: omega, Gam
    real(kind=8), intent(inout)       :: delL, delR, IL, IR



    character(len=100) :: buffer, local
    character(len=100) :: label
    integer :: pos
    integer, parameter :: fh=15
    integer :: ios=0, line=0


    open(unit=fh, file='problem_data.dat',action='read')

    do while(ios==0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios==0) then
        line = line + 1

        ! find the first instance of whitespace. Split and label data
        pos=scan(buffer,' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)
        select case(label)
        case('ion_mass')
              read(buffer,*,iostat=ios) mass
        case('ion_charge')
              read(buffer,*,iostat=ios) charge
        case('line_freq')
              read(buffer,*,iostat=ios) omega
        case('line_width')
              read(buffer,*,iostat=ios) Gam
        case('detuning_left')
              read(buffer,*,iostat=ios) delL
        case('detuning_right')
              read(buffer,*,iostat=ios) delR
        case('intensity_left')
            read(buffer,*,iostat=ios) IL
        case('intensity_left')
            read(buffer,*,iostat=ios) IR
        case('nbath')
            read(buffer,*,iostat=ios) nbath
        case('nparticles')
            read(buffer,*,iostat=ios) nparticles
        case('dim')
            read(buffer,*,iostat=ios) dim
        case default
              print*, "Skipping invalid value."
        end select
      end if
    end do

    Gam  = 2.0d0 * pi * Gam
    del1  = del1 * Gam
    del2  = del2 * Gam
    omega0 = 2.0d0 * pi * omega0

    close(unit=fh)

  end subroutine initialize_problem

end module initialization
