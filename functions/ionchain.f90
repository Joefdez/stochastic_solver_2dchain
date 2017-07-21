module ionchain
  use constants
  implicit none

  contains

  subroutine coulombM(nparticles, xx, yy, fx, fy, invD)
    implicit none
    real(kind=8), dimension(:), intent(in)              :: xx, yy
    real(kind=8), dimension(:,:), intent(inout)         :: fx, fy, invD
    integer, intent(in)                                 :: nparticles
    integer                                             :: ii, jj
    real(kind=8)                                        :: dist
    invD = 0.0d0
    fx = 0.0d0
    fy = 0.0d0
    do ii=1, nparticles, 1
      do jj=ii+1, nparticles, 1
        dist = 0.0d0
        dist = sqrt( (xx(ii)-xx(jj))**2 + (yy(ii)-yy(jj))**2 )
        invD(ii,jj) = 1.0d0/dist
        invD(jj,ii) = 1.0d0/dist
        dist = dist*dist*dist
        fx(ii,jj) = (xx(ii)-xx(jj)) / dist
        fx(jj,ii) = -1.0d0*(xx(ii)-xx(jj)) / dist
        fy(ii,jj) = (yy(ii)-yy(jj)) / dist
        fy(jj,ii) = -1.0d0*(yy(ii)-yy(jj)) / dist
      end do
    end do
  end subroutine coulombM

  subroutine local_energy(nparticles, alpha, xx, yy, invD, ppx, ppy, energy)
    implicit none
    integer, intent(in)                                   :: nparticles
    real(kind=8), intent(in)                              :: alpha
    real(kind=8), dimension(:), intent(in)                :: xx, yy, ppx, ppy
    real(kind=8), dimension(:,:), intent(in)              :: invD
    real(kind=8), dimension(:), intent(inout)             :: energy
    real(kind=8)                                          :: kin_term
    integer                                               :: ii
    energy = 0.0d0

    ! Add kinetic energy contribution and harmonic potential energy contributions
    energy(1:nparticles) = 0.5*ppx(1:nparticles)*ppx(1:nparticles) + 0.5d0*ppy(1:nparticles)*ppy(1:nparticles) &
              + 0.5d0*xx(1:nparticles)*xx(1:nparticles) + 0.5d0*alpha*alpha*yy(1:nparticles)*yy(1:nparticles)
    do ii=1, nparticles, 1
      energy(ii) = energy(ii) + 0.5d0*sum(invD(ii,:), 1)          ! Add the coulomb potential contribution
    end do

  end subroutine local_energy

  subroutine heat_current(nparticles, cfx, cfy, ppx, ppy, hc)
    implicit none
    integer, intent(in)                                   :: nparticles
    real(kind=8), dimension(:,:), intent(in)              :: cfx, cfy
    real(kind=8), dimension(:), intent(in)                :: ppx, ppy
    real(kind=8), dimension(:,:), intent(inout)              :: hc
    integer                                               :: ii, jj
    hc=0.0d0
    do ii=1, nparticles, 1
      do jj=ii+1, nparticles, 1
        hc(ii,jj) = 0.5d0*(ppx(ii)+ppx(jj))*cfx(ii,jj) + 0.5d0*(ppy(ii)+ppy(jj))*cfy(ii,jj)
        hc(jj,ii) = -1.0d0*hc(ii,jj)
      end do
    end do

  end subroutine heat_current

  subroutine current_Flux(hc, energy, xx, yy, ppx, ppy, nparticles, JJintx, JJinty)
    implicit none
    real(kind=8), dimension(:,:), intent(in)              :: hc !cfx, cfy
    real(kind=8), dimension(:), intent(in)                :: xx, yy, ppx, ppy, energy
    integer, intent(in)                                   :: nparticles
    real(kind=8), intent(inout)                           :: JJintx, JJinty
    real(kind=8)                                          :: JJ0
    integer                                               :: nn, ll

    JJintx = 0.0d0
    JJinty = 0.0d0

    do nn=1, nparticles, 1
      JJintx = JJintx + energy(nn)*ppx(nn)
      JJinty = JJinty + energy(nn)*ppy(nn)
      do ll=1, nparticles, 1
        if(nn .eq. ll) cycle
        !JJ0 = 0.5d0*(ppx(nn+1)+ppx(ll))*cfx(nn+1,ll) + 0.5d0*(ppy(nn+1)+ppy(ll))*cfy(nn+1,ll)
        JJintx = JJintx + 0.5d0*(xx(nn)-xx(ll))*hc(nn,ll)
        JJinty = JJinty + 0.5d0*(yy(nn)-yy(ll))*hc(nn,ll)
      end do
    end do

  end subroutine current_Flux

  subroutine doppler_values(kk, gam, del, II, eta, DD)
    implicit none
    real(kind=8), intent(in)    :: kk, gam, del, II
    real(kind=8), intent(inout) :: eta, DD

    eta = -4.0d0*hbar*kk*kk*II*(2.0d0*del/Gam)/( (1 + 4.0d0*del*del/(Gam*Gam)) * (1 + 4.0d0*del*del/(Gam*Gam)) )
    DD   =  hbar*hbar*kk*kk*II*(Gam)/(1.0d0 + 4.0d0*del*del/(Gam*Gam))

  end subroutine doppler_values


  subroutine dimensionless_doppler_values(eta, D, mass, long_freq, char_length, aeta, aD)
    implicit none
    real(kind=8), intent(in)              :: eta, D, mass, long_freq, char_length
    real(kind=8), intent(inout)           :: aeta, aD
    ! calculate scaled friction and difussion coefficients
    aeta = eta / (mass * long_freq)
    aD   = D / (char_length*char_length*mass*mass*long_freq*long_freq*long_freq)

  end subroutine dimensionless_doppler_values


end module ionchain
