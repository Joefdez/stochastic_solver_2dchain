
  subroutine platen_realization_simpleB(YY0, YY, YYf, BB, nbath, dt, dst, nsteps, nsample_steps, final_steps)
    implicit none
    ! Calculates a realization of the stochastic process using the weak Platen scheme
    ! Simplified version using a vector for BB. 
    real(kind=8), intent(in), dimension(:)          :: YY0
    real(kind=8), intent(inout), dimension(:,:)     :: YY, YYf
    real(kind=8), dimension(:), intent(in)          :: BB
    integer, intent(in)                             :: nbath, neqs
    real(kind=8), intent(in)                        :: dt, dst, nsteps, nsample_steps, final_steps

    integer                                         :: ii, jj, kk
    real(kind=8), dimension(1:neqs)                 :: YYold, YYnew, YYi, AA, AAi, dOm
    real(kind=8), dimension(1:(neqs/2), 1:(neqs/2)) :: fx1, fx2, invD
    !allocate(YYi(1:neqs))
    !allocate(AA(1:neqs))
    !allocate(AAi(1:neqs))
    !allocate(dOm(1:neqs))
    !allocate(BB(1:neqs, 1:neqs))
    half    = neqs/2
    quarter = neqs/4
    YYold = YY0
    jj = 1
    kk = 1
    do ii=1, nsteps, 1
      call coulombM(npart, YYold(ii,1:quarter), YYold(ii,(quarter+1):half), fx1, fy2, invD1)
      call deterministic_terms(YYold, AA, neqs, nbath, cfx, cfy, alpha, eta1, eta2)
      call stochastic_step(dOm, dst)
      YYi = YYold + AA*dt + BB*dOm
      call coulombM(npart, YYi(1:quarter), YYi((quarter+1):half), fx2, fy2, invD2)
      call deterministic_terms(YYi, AAi, neqs, nbath, cfx, cfy, alpha, eta1, eta2)
      YYnew = YYold + 0.5d0*(AA+AAi)*dt + BB*dOm
      YYold = YYnew
      if(mod(ii,nsample_steps) .eq. 0) then                                     ! Save selected sample
        jj = jj + 1
        YY(jj) = YYnew
      end if
      if(ii .gt. final_steps) then
        kk = kk + 1
        YYf(kk) = YYnew
      end if
    end do

  end subroutine platen_realization
