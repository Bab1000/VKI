!------------------------------------------------------------------------------!
!> This subroutine associates the pointers for the following subroutines (1D calorically perfect gas flows):
!! - convective flux 
!! - diffusive flux 
!! - convective flux Jacobian approximation
!! - convective (and eventually) diffusive flux Jacobian spectral radius
  SUBROUTINE set_flux_1D ()

    USE mod_general_data,      ONLY: flag_diss
    USE mod_numerics_data,     ONLY: flux_splitter, flux_Jac_approx, flag_full_Impl, diff_flux_Jac
    USE mod_function_pointer

    IMPLICIT NONE
 
    ! Convective flux 
    SELECT CASE(flux_splitter)

      ! van Leer flux vector splitting 
      CASE('van_leer','van_Leer')
        conv_flux_1D => van_Leer_1D             

      ! Steger-Warming flux vector splitting
      CASE('steg_war','Steg_War')
        conv_flux_1D => Steger_Warming_1D  

      ! Roe's approximate Riemann solver
      CASE('roe','Roe')
        conv_flux_1D => roe_1D   

      ! HLLE's approximate Riemann solver
      CASE('hlle','HLLE')
        conv_flux_1D => hlle_1D   

      CASE DEFAULT 
        WRITE(*,10)'In set_numerical_flux_1D.F90, error in numerical flux function selection ...'
        PRINT*
        STOP 

      END SELECT 

      ! Numerical flux Jacobian (only in case of a fully implicit time-integration method) 
      IF (flag_full_Impl.EQV..TRUE.) THEN

         SELECT CASE(flux_Jac_approx) 

           CASE('Numerical','numerical')
             conv_flux_Jac_1D => conv_flux_num_Jac_1D
 
           CASE('Yoon_Jameson','yoon_jameson')
             conv_flux_Jac_1D => Yoon_Jameson_Jac_1D

           CASE('pos_neg_split','Pos_Neg_split')
             conv_flux_Jac_1D => Pos_Neg_split_Jac_1D

           CASE DEFAULT
             WRITE(*,10)'In "set_numerical_flux_1D.F90", error in flux Jacobian approximation selection...'
             PRINT*
             STOP 

         END SELECT

      ENDIF

      ! Diffusive flux (and related Jacobian)
      IF (flag_diss.EQV..TRUE.) THEN

         diff_flux_1D => ns_flux_1D
         IF (diff_flux_Jac=='numerical') THEN
            diff_flux_1D_Jac => diff_flux_num_Jac_1D
         ELSEIF (diff_flux_Jac=='analytical') THEN
            diff_flux_1D_Jac => ns_flux_1D_Jac
         ELSE 
            WRITE(*,10)'In "set_num_flux_1D.F90", error in diffusive flux Jacobian approximation selection ...'
            PRINT*
            STOP 
         ENDIF

      ELSE 

         diff_flux_1D => null_ns_flux_1D
         diff_flux_1D_Jac => null_ns_flux_1D_Jac

      ENDIF

      ! Spectral radius 
      IF (flag_diss.EQV..TRUE.) THEN
         
         get_inv_spectral_radius_1D => inv_spectral_radius_1D 
         get_visc_spectral_radius_1D => visc_spectral_radius_1D 

      ELSE

         get_inv_spectral_radius_1D => inv_spectral_radius_1D 
         get_visc_spectral_radius_1D => null_visc_spectral_radius_1D 

      ENDIF

10 FORMAT(A)

  END SUBROUTINE set_flux_1D
!------------------------------------------------------------------------------!
