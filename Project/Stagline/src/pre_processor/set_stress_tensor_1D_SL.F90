!------------------------------------------------------------------------------!
!> This subroutine associates the pointer for the stress tensor subroutine.
   SUBROUTINE set_stress_tensor_1D_SL ()

     USE mod_numerics_data,             ONLY: stag_line_geom 
     USE mod_function_pointer

     IMPLICIT NONE

     SELECT CASE(stag_line_geom)

       ! Cylinder
       CASE(0)
         get_stress_tensor_1D_SL => stress_tensor_1D_SL_cyl

       ! Sphere
       CASE(1)
         get_stress_tensor_1D_SL => stress_tensor_1D_SL_sph

     END SELECT

   END SUBROUTINE set_stress_tensor_1D_SL
!------------------------------------------------------------------------------!
