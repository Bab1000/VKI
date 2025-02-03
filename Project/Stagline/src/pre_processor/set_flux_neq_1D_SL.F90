!------------------------------------------------------------------------------!
!> This subroutine associates the pointers for the following subroutines (1D stagnation line nonequilibrium flows):
!! - convective flux 
!! - diffusive flux 
!! - convective flux Jacobian approximation
!! - convective (and eventually) diffusive flux Jacobian spectral radius
!! - eigensystem
  SUBROUTINE set_flux_neq_1D_SL ()

    USE mod_general_data,      ONLY: nb_temp, nb_int_temp, nb_te, flag_diss
    USE mod_numerics_data,     ONLY: flux_splitter, flux_Jac_approx, flag_full_Impl, diff_flux_Jac, & 
                                   & stag_line_geom, flag_metrics
    USE mod_function_pointer
 
    IMPLICIT NONE

    ! Numerical flux 
    SELECT CASE(flux_splitter)

      ! van Leer's flux vector splitting 
      CASE('van_leer','van_Leer')
        conv_flux_1D_SL => van_Leer_neq_1D_SL            
 
      ! Steger-Warming's flux vector splitting 
      CASE ('steg_war','Steg_War')
        conv_flux_1D_SL => StegerWarming_neq_1D_SL 
        
        ! Eigensystem  
        SELECT CASE(nb_temp)

          ! 1 Temperature 
          CASE(1) 
            eigsys_neq_1D_SL  => eigensystem_neq1T_1D_SL 

          ! N temperatures and eventually one for free electrons (Te)
          CASE DEFAULT 
            
            IF (nb_te.EQ.0) THEN
               eigsys_neq_1D_SL  => eigensystem_neqNT_1D_SL
            ELSE
               eigsys_neq_1D_SL  => eigensystem_neqNT_Te_1D_SL
            ENDIF

        END SELECT 

      ! AUSM+up All Speed Scheme
      CASE('AUSMP_UP_AS', 'ausmP_up_AS', "ausmp_up_as")
        conv_flux_1D_SL => ausmP_up_as_neq_1D_SL

      ! AUSM+up All Speed Scheme
      CASE('AUSMP_UP_AS_CD', 'ausmP_up_AS_CD', "ausmp_up_as_cd")
        conv_flux_1D_SL => ausmP_up_as_cd_neq_1D_SL

      ! AUSM+up_2 All Speed Scheme
      CASE('AUSMP_UP_AS_2', 'ausmP_up_AS_2', "ausmp_up_as_2")
        conv_flux_1D_SL => ausmP_up2_neq_1D_SL


     ! ! Preconditioned AUSM 
     ! CASE('AUSMPC', 'ausmPC')
     !   conv_flux_1D_SL => ausmPC_neq_1D_SL
        
     ! CASE('LDFSS(2)', 'ldfss(2)')
     !   conv_flux_1D_SL => LDFSS2_neq_1D_SL

      ! Roe approximate Riemann solver
      CASE('roe','Roe')
        conv_flux_1D_SL => roe_neq_1D_SL

        ! Eigensystem  
        SELECT CASE(nb_temp)

          ! 1 Temperature 
          CASE(1) 
            eigsys_neq_1D_SL  => eigensystem_neq1T_1D_SL
            roe_avg_neq_1D_SL => roe_avgState_neq1T_1D_SL

          ! N temperatures and eventually one for free electrons (Te)
          CASE DEFAULT 
            
            IF (nb_te.EQ.0) THEN
               eigsys_neq_1D_SL  => eigensystem_neqNT_1D_SL
               roe_avg_neq_1D_SL => roe_avgState_neqNT_1D_SL
            ELSE
               eigsys_neq_1D_SL  => eigensystem_neqNT_Te_1D_SL
               roe_avg_neq_1D_SL => roe_avgState_neqNT_Te_1D_SL
            ENDIF

        END SELECT 
      ! Roe approximate Riemann solver
      CASE('modified_SW', 'MODIFIED_SW')
        conv_flux_1D_SL => modified_SW_neq_1D_SL

        ! Eigensystem  
        SELECT CASE(nb_temp)

          ! 1 Temperature 
          CASE(1)

               eigsys_neq_1D_SL  => eigensystem_neq1T_1D_SL                    
            
          ! N temperatures and eventually one for free electrons (Te)
          CASE DEFAULT 
            
            IF (nb_te.EQ.0) THEN
               eigsys_neq_1D_SL  => eigensystem_neqNT_1D_SL
            ELSE
               eigsys_neq_1D_SL  => eigensystem_neqNT_Te_1D_SL
            ENDIF

        END SELECT   

      ! HLLE's approximate Riemann solver 
      CASE('hlle','HLLE')
        !conv_flux_1D_SL => hlle_neq_1D_SL

        ! Eigensystem  
        SELECT CASE(nb_temp)

          ! 1 Temperature 
          CASE(1) 
            roe_avg_neq_1D_SL => roe_avgState_neq1T_1D_SL

          ! N temperatures and eventually one for free electrons (Te)
          CASE DEFAULT 
            
            IF (nb_te.EQ.0) THEN
               roe_avg_neq_1D_SL => roe_avgState_neqNT_1D_SL
            ELSE
               roe_avg_neq_1D_SL => roe_avgState_neqNT_Te_1D_SL
            ENDIF

        END SELECT  

      CASE DEFAULT 

         WRITE(*,10)'In "set_flux_neq_1D_SL.F90", error in numerical flux function selection ...'
         PRINT*
         STOP 

      END SELECT 

      ! Inviscid flux Jacobian
      SELECT CASE(nb_temp)

        ! 1 Temperature 
        CASE(1) 
          inv_flux_Jac_neq_1D_SL => inviscid_flux_Jac_neq1T_1D_SL

        ! N temperatures and eventually one for free electrons (Te)
        CASE DEFAULT 
            
          IF (nb_te.EQ.0) THEN
             inv_flux_Jac_neq_1D_SL => inviscid_flux_Jac_neqNT_1D_SL
          ELSE
             inv_flux_Jac_neq_1D_SL => inviscid_flux_Jac_neqNT_Te_1D_SL
          ENDIF

      END SELECT 

      ! Numerical flux Jacobian (only in case of a fully implicit time-integration method)
      IF (flag_full_Impl.EQV..TRUE.) THEN

        SELECT CASE(flux_Jac_approx) 

          CASE('Numerical','numerical')
            conv_flux_Jac_1D_SL => conv_flux_num_Jac_1D_SL

          CASE('Yoon_Jameson')
            conv_flux_Jac_1D_SL => Yoon_Jameson_Jac_neq_1D_SL

          CASE('pos_neg_split','Pos_Neg_split')
            SELECT CASE(nb_temp)

              ! 1 Temperature 
              CASE(1) 
                conv_flux_Jac_1D_SL => Pos_Neg_split_Jac_neq1T_1D_SL

              ! N temperatures and eventually one for free electrons (Te)
              CASE DEFAULT 

                IF (nb_te.EQ.0) THEN
                   conv_flux_Jac_1D_SL => Pos_Neg_split_Jac_neqNT_1D_SL
                ELSE
                   conv_flux_Jac_1D_SL => Pos_Neg_split_Jac_neqNT_Te_1D_SL 
                ENDIF

            END SELECT
 
            CASE('pos_neg_general_split','Pos_Neg_General_split')
              conv_flux_Jac_1D_SL=> Pos_Neg_split_Jac_general_1D_SL
              ! Eigensystem  
            SELECT CASE(nb_temp)

              ! 1 Temperature 
              CASE(1)

                   eigsys_neq_1D_SL  => eigensystem_neq1T_1D_SL                    
                
              ! N temperatures and eventually one for free electrons (Te)
              CASE DEFAULT 
                
                IF (nb_te.EQ.0) THEN
                   eigsys_neq_1D_SL  => eigensystem_neqNT_1D_SL
                ELSE
                   eigsys_neq_1D_SL  => eigensystem_neqNT_Te_1D_SL
                ENDIF

            END SELECT 
  
          CASE DEFAULT
            WRITE(*,10)'In "set_num_flux_neq_1D_SL.F90", error in flux Jacobian approximation selection ...'
            PRINT*
            STOP 

        END SELECT

      ENDIF

      ! Diffusive flux (and related Jacobians)
      IF (flag_diss.EQV..TRUE.) THEN

         SELECT CASE(stag_line_geom)

           ! Cylinder 
           CASE(0)

             ! 1 temperature
             IF (nb_temp.EQ.1) THEN

                diff_flux_1D_SL => ns_flux_neq1T_1D_SL_cyl

             ! N temperatures and eventually one for free electrons (Te)
             ELSE 

                IF (nb_te.EQ.0) THEN
                   diff_flux_1D_SL => ns_flux_neqNT_1D_SL_cyl
                ELSE 
                   diff_flux_1D_SL => ns_flux_neqNT_Te_1D_SL_cyl
                ENDIF

             ENDIF

             ! Numerical diffusive flux Jacobian
             IF (diff_flux_Jac.EQ.'numerical') THEN

                diff_flux_Jac_1D_SL => diff_flux_num_Jac_1D_SL

             ! Analytical diffusive flux Jacobian
             ELSEIF (diff_flux_Jac.EQ.'analytical') THEN
 
                ! 1 temperature
                IF (nb_temp.EQ.1) THEN

                   diff_flux_Jac_1D_SL => ns_flux_neq1T_1D_SL_cyl_Jac

                ! N temperatures and eventually one for free electrons (Te)
                ELSE 

                   IF (nb_te.EQ.0) THEN
                      diff_flux_Jac_1D_SL => ns_flux_neqNT_1D_SL_cyl_Jac
                   ELSE 
                      diff_flux_Jac_1D_SL => ns_flux_neqNT_Te_1D_SL_cyl_Jac
                   ENDIF

                ENDIF

             ELSE 
                WRITE(*,10)'In "set_numerical_flux_neq_1D_SL.F90", error in diffusive flux Jacobian approximation selection...'
                PRINT*
                STOP 
             ENDIF

           ! Sphere 
           CASE(1)

             ! 1 temperature
             IF (nb_temp.EQ.1) THEN
                IF (flag_metrics.EQV..FALSE.) THEN 
                   diff_flux_1D_SL => ns_flux_neq1T_1D_SL_sph
                ELSE
                   diff_flux_1D_SL => ns_flux_neq1T_1D_SL_sph_metr
                ENDIF


             ! N temperatures and eventually one for free electrons (Te)
             ELSE 
                IF (flag_metrics.EQV..FALSE.) THEN 
                   IF (nb_te.EQ.0) THEN
                      diff_flux_1D_SL => ns_flux_neqNT_1D_SL_sph 
                   ELSE 
                      diff_flux_1D_SL => ns_flux_neqNT_Te_1D_SL_sph 
                   ENDIF
                ELSE 
                   IF (nb_te.EQ.0) THEN
                      diff_flux_1D_SL => ns_flux_neqNT_1D_SL_sph_metr 
                   ELSE 
                      diff_flux_1D_SL => ns_flux_neqNT_Te_1D_SL_sph_metr 
                   ENDIF

                ENDIF

             ENDIF

             ! Numerical diffusive flux Jacobian
             IF (diff_flux_Jac.EQ.'numerical') THEN
              
                diff_flux_Jac_1D_SL => diff_flux_num_Jac_1D_SL

             ! Analytical diffusive flux Jacobian
             ELSEIF (diff_flux_Jac.EQ.'analytical') THEN
 
                ! 1 temperature
                IF (nb_temp.EQ.1) THEN
                   IF (flag_metrics.EQV..FALSE.) THEN 
                      diff_flux_Jac_1D_SL => ns_flux_neq1T_1D_SL_sph_Jac
                   ELSE
                      diff_flux_Jac_1D_SL => ns_flux_neq1T_1D_SL_sph_Jac_metr
                   ENDIF

                ! N temperatures and eventually one for free electrons (Te)
                ELSE 

                   IF (nb_te.EQ.0) THEN
                      diff_flux_Jac_1D_SL => ns_flux_neqNT_1D_SL_sph_Jac
                   ELSE 
                      diff_flux_Jac_1D_SL => ns_flux_neqNT_Te_1D_SL_sph_Jac
                   ENDIF

                ENDIF

             ELSE 
                WRITE(*,10)'In "set_numerical_flux_neq_1D_SL.F90", error in diffusive flux Jacobian approximation selection...'
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
         get_visc_spectral_radius_1D_SL => visc_spectral_radius_neq_1D_SL

      ELSE

         get_inv_spectral_radius_1D_SL => inv_spectral_radius_1D_SL 
         get_visc_spectral_radius_1D_SL => null_visc_spectral_radius_1D_SL 

      ENDIF

10 FORMAT(A)

  END SUBROUTINE set_flux_neq_1D_SL
!------------------------------------------------------------------------------!
