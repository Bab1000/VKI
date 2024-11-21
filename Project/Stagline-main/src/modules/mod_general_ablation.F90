!------------------------------------------------------!
! This module contains data regarding ablation
!------------------------------------------------------!
!> @author
!> Bruno Dias
!>
!>@brief
!>Contains general data:w
!>
!>\b DESCRIPTION: 
!>
!> This module contains data regarding ablation
!------------------------------------------------------!
  MODULE mod_general_ablation

  IMPLICIT NONE


  REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: mdot_i          !< wall production terms \f$ mdot_i \f$ (working vector)

  REAL(KIND=8), SAVE :: mdot
  REAL(KIND=8), SAVE :: uwall
  LOGICAL, SAVE :: Flag_ablation      = .FALSE.
  
  REAL(KIND=8), SAVE :: Twall_translation
  REAL(KIND=8), SAVE :: Pwall_ablation
  


 END MODULE mod_general_ablation

 MODULE mod_surface_properties

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: rho_i_surface
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: T_surface
  REAL(KIND=8), SAVE :: u_surface
  REAL(KIND=8), SAVE :: v_surface
  

 END MODULE mod_surface_properties
