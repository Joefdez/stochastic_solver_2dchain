module deterministic
  use constants
  implicit none

  contains

  subroutine deterministic_terms(YYn, AA, neqs, nbath, cfx, cfy, alpha, eta1, eta2)
    implicit none
    ! 2D ion chain with nbath ions on each edge connected to heat baths. Adimensional units. See references.
    real(kind=8), dimension(:), intent(in)                      :: YYn, cfx, cfy
    real(kind=8), dimension(:), intent(inout)                   :: AA
    integer, intent(in)                                         :: neqs, nbath
    real(kind=8), intent(in)                                    :: alpha, eta1, eta2
    integer                                                     :: ii, half, quarter

    half    = neqs/2
    quarter = neqs/4
    AA = 0.0d0

    ! RHS of position equations
    do ii=1, half, 1
      AA(ii) = YYn(half+ii)
    end do

    ! RHS of the momentum equations
    do ii=1, quarter, 1                           ! Contribution from harmonic and Coulomb terms in x
      AA(half+ii) = -1.0d0*YYn(ii) + cfx(ii)
    end do
    do ii=1, quarter, 1                           ! Contribution from harmonic and Coulomb terms in y
      AA(3*quarter + ii) = -1.0d0*YYn(quarter+ii)*alpha*alpha + cfy(ii)
    end do

    do ii=1, nbath, 1                             ! Contribution from left laser on x
      AA(half+ii) = AA(half+ii) - eta1*YYn(half+ii)
    end do
    do ii=1, nbath, 1                             ! Contribution from left laser on y
      AA(half+quarter+ii) = AA(half+quarter+ii) - eta1*YYn(half+quarter+ii)
    end do
    do ii=0, nbath-1, 1                          ! Contribution from right laser on x
      AA(neqs-quarter-ii) = AA(neqs-quarter-ii) - eta2*YYn(neqs-quarter-ii)
    end do
    do ii=0, nbath-1, 1                          ! Contribution from right laser on y
      AA(neqs-quarter-ii) = AA(neqs-ii) - eta2*YYn(neqs-ii)
    end do

  end subroutine deterministic_terms

end module deterministic
