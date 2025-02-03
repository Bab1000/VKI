!------------------------------------------------------------------------------!
!> This subroutine associates the pointer to the subroutine computing transport 
!! coefficients (calorically perfect gas flows)
  SUBROUTINE set_transp_coeff ()

    USE mod_general_data,       ONLY: viscosity_law, flag_diss
    USE mod_function_pointer 

    IMPLICIT NONE

    IF (flag_diss.EQV..TRUE.) THEN
 
       SELECT CASE(viscosity_law)

         ! Sutherland's law for viscosity
         CASE('Sutherland') 
           get_transpCoeff => transpCoeff_sutherland    

         ! Inverse power (interaction potential)
         CASE('Inv_power')
           get_transpCoeff => transpCoeff_inv_power    

         CASE DEFAULT
           WRITE(*,10)'In set_transp_coeff.F90, viscosity law not implemented yet...'
           WRITE(*,10)viscosity_law
           PRINT*
           STOP 

       END SELECT

    ENDIF 

10 FORMAT(A)

  END SUBROUTINE set_transp_coeff
!------------------------------------------------------------------------------!
