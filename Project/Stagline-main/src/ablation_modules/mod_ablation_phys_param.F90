!------------------------------------------------------------------------------!
!> This modules stores the physical parameters used when the ablative
!! boundary condition is used.
MODULE mod_ablation_phys_param

#include "../config.h"

#ifdef CARBONABLA

        IMPLICIT NONE 

        ! Parameters for the Surface Energy Balance (SEB)
        REAL(KIND=8), PARAMETER :: SF_const   = 5.670373d-8        !< Stefan–Boltzmann constant

#endif

END MODULE mod_ablation_phys_param
