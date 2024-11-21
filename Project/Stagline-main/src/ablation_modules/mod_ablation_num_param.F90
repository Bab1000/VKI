!------------------------------------------------------------------------------!
!> This modules stores the numerical parameters used when the ablative
!! boundary condition is used.
MODULE mod_ablation_num_param

#include "../config.h"

#ifdef CARBONABLA

        IMPLICIT NONE 

        INTEGER, PARAMETER :: j_max  = 3                                                          !< Newton procedure step number
        REAL(KIND=8), PARAMETER :: eps = 1.d-10                                                   !< Newton procedure perturbation value

#endif

END MODULE mod_ablation_num_param
