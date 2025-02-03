!------------------------------------------------------------------------------!
!> This function defines the van leer limiter
  FUNCTION van_leer (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: van_leer

    REAL(KIND=8) :: tmp
 
    tmp = ABS(r)
    van_leer = (r + tmp)/(1.d0 + tmp)

  END FUNCTION van_leer
!------------------------------------------------------------------------------!
!> This function defines the minimod limiter
  FUNCTION minimod (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: minimod

    minimod = MAX(0.d0,MIN(1.d0,r))

  END FUNCTION minimod
!------------------------------------------------------------------------------!
!? This function defines the superbee limiter
  FUNCTION superbee (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: superbee

    superbee = MAX(0.d0,MAX(MIN(2.d0*r,1.d0),MIN(r,2.d0)))

  END FUNCTION superbee
!------------------------------------------------------------------------------!
!> This function defines the koren limiter
  FUNCTION koren (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: koren

    koren = MAX(0.d0,MIN(MIN(2.d0*r,(1.D0 + 2.d0*r)/3.d0),2.d0))

  END FUNCTION koren
!------------------------------------------------------------------------------!
!> This function defines the mc limiter
  FUNCTION mc (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: mc

    mc = MAX(0.d0,MIN(MIN(2.D0*r,(1.d0 + r)/2.d0),2.d0))

  END FUNCTION mc
!------------------------------------------------------------------------------!
!> This function defines the charm limiter
  FUNCTION charm (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: charm

    charm = r*(3.d0*r + 1.d0)/(r + 1.d0)**2.d0

  END FUNCTION charm
!------------------------------------------------------------------------------!
!> This function defines the hcus limiter
  FUNCTION hcus (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: hcus

    hcus = 1.5D0*(r + ABS(r))/(r + 2.d0)

  END FUNCTION hcus
!------------------------------------------------------------------------------!
!> This function defines the hquick limiter
  FUNCTION hquick (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: hquick

    hquick = 2.D0*(r + ABS(r))/(r + 3.d0)

  END FUNCTION hquick 
!------------------------------------------------------------------------------!
!> This function defines the first version of van Albada limiter
  FUNCTION van_albada1 (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: van_albada1

    REAL(KIND=8) :: tmp

    tmp = r*r
    van_albada1 = (tmp + r)/(tmp + 1.d0)

  END FUNCTION van_albada1 
!------------------------------------------------------------------------------!
!> This function defines the second version of van Albada limiter
  FUNCTION van_albada2 (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: van_albada2

    van_albada2 = 2.d0*r/(r*r + 1.d0)

  END FUNCTION van_albada2
!------------------------------------------------------------------------------!
!> This function defines the ospre limiter
  FUNCTION ospre (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: ospre

    REAL(KIND=8) :: tmp
  
    tmp = r*r
    ospre = 1.5d0*(tmp + r)/(tmp + r + 1.d0)

  END FUNCTION ospre
!------------------------------------------------------------------------------!
!> This function defines the no limiter (gradients are not limited) 
  FUNCTION no_limiter (r)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r
    REAL(KIND=8) :: no_limiter

    no_limiter = 1.d0

  END FUNCTION no_limiter
!------------------------------------------------------------------------------!
