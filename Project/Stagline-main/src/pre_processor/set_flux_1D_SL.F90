!------------------------------------------------------------------------------!
!> This subroutine associates the pointers for the following subroutines (1D calorically perfect gas flows):
!! - convective flux 
!! - diffusive flux 
!! - convective flux Jacobian approximation
!! - convective (and eventually) diffusive flux Jacobian spectral radius
  SUBROUTINE set_flux_1D_SL ()

    USE mod_general_data,      ONLY: flag_diss
    USE mod_numerics_data,     ONLY: flux_splitter, flux_Jac_approx, flag_full_Impl, diff_flux_Jac, & 
                                   & stag_line_geom, flag_metrics
    USE mod_function_pointer

    IMPLICIT NONE
 
    ! Convective flux 
    SELECT CASE(flux_splitter)

      ! van Leer flux vector splitting 
      CASE('van_leer','van_Leer')
        conv_flux_1D_SL => van_Leer_1D_SL            

      ! Steger-Warming flux vector splitting
      CASE('steg_war','Steg_War')
        conv_flux_1D_SL => Steger_Warming_1D_SL  

      ! Roe's approximate Riemann solver
      CASE('roe','Roe')
        conv_flux_1D_SL => roe_1D_SL   

      ! Preconditioned AUSM 
      CASE('ausmPC','AUSMPC')
        conv_flux_1D_SL => ausmPC_1D_SL  

      ! HLLE's approximate Riemann solver
      !CASE('hlle','HLLE')
      !  conv_flux_1D_SL => hlle_1D_SL   

      CASE DEFAULT 
        WRITE(*,10)'In "set_numerical_flux_1D_SL.F90", error in numerical flux function selection...'
        PRINT*
        STOP 

    END SELECT 

      ! Numerical flux Jacobian (only in case of a fully implicit time-integration method) 
    IF (flag_full_Impl.EQV..TRUE.) THEN

       SELECT CASE(flux_Jac_approx) 

         CASE('Numerical','numerical')
           conv_flux_Jac_1D_SL => conv_flux_num_Jac_1D_SL
 
         CASE('Yoon_Jameson','yoon_jameson')
           conv_flux_Jac_1D_SL => Yoon_Jameson_Jac_1D_SL

         CASE('pos_neg_split','Pos_Neg_split')
           conv_flux_Jac_1D_SL => Pos_Neg_split_Jac_1D_SL
        
         CASE DEFAULT
           WRITE(*,10)'In "set_numerical_flux_1D_SL.F90", error in flux Jacobian approximation selection...'
           PRINT*
           STOP 

       END SELECT

    ENDIF

    ! Diffusive flux (diffusive flux Jacobians to be added)
    IF (flag_diss.EQV..TRUE.) THEN

       SELECT CASE(stag_line_geom)

         ! Cylinder 
         CASE(0)
           diff_flux_1D_SL => ns_flux_1D_SL_cyl
           IF (diff_flux_Jac.EQ.'numerical') THEN
              diff_flux_Jac_1D_SL => diff_flux_num_Jac_1D_SL
           ELSEIF (diff_flux_Jac.EQ.'analytical') THEN
              diff_flux_Jac_1D_SL => ns_flux_1D_SL_cyl_Jac
           ELSE 
              WRITE(*,10)'In "set_numerical_flux_1D_SL.F90", error in diffusive flux Jacobian approximation selection...'
              PRINT*
              STOP 
           ENDIF

         ! Sphere
         CASE(1)
           IF (flag_metrics.EQV..FALSE.) THEN 
              diff_flux_1D_SL => ns_flux_1D_SL_sph
           ELSE
              diff_flux_1D_SL => ns_flux_1D_SL_sph_metr
           ENDIF
           IF (diff_flux_Jac.EQ.'numerical') THEN
              diff_flux_Jac_1D_SL => diff_flux_num_Jac_1D_SL
           ELSEIF (diff_flux_Jac.EQ.'analytical') THEN
              IF (flag_metrics.EQV..FALSE.) THEN 
                 diff_flux_Jac_1D_SL => ns_flux_1D_SL_sph_Jac
              ELSE
                 diff_flux_Jac_1D_SL => ns_flux_1D_SL_sph_Jac_metr
              ENDIF
           ELSE 
              WRITE(*,10)'In "set_numerical_flux_1D_SL.F90", error in diffusive flux Jacobian approximation selection...'
              PRINT*
              STOP  
           ENDIF

       END SELECT

    ELSE

       diff_flux_1D_SL => null_ns_flux_1D_SL 
       diff_flux_Jac_1D_SL => null_ns_flux_1D_SL_Jac

    ENDIF

    ! Spectral radius 
    IF (flag_diss.EQV..TRUE.) THEN
       
       get_inv_spectral_radius_1D_SL => inv_spectral_radius_1D_SL 
       get_visc_spectral_radius_1D_SL => visc_spectral_radius_1D_SL

    ELSE

       get_inv_spectral_radius_1D_SL => inv_spectral_radius_1D_SL 
       get_visc_spectral_radius_1D_SL => null_visc_spectral_radius_1D_SL 

    ENDIF

10 FORMAT(A)

  END SUBROUTINE set_flux_1D_SL
!------------------------------------------------------------------------------!
