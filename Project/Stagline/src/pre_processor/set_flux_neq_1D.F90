!------------------------------------------------------------------------------!
!> This subroutine associates the pointers for the following subroutines (1D nonequilibrium flows):
!! - convective flux 
!! - diffusive flux 
!! - convective flux Jacobian approximation
!! - convective (and eventually) diffusive flux Jacobian spectral radius
!! - eigensystem
  SUBROUTINE set_flux_neq_1D ()

    USE mod_general_data,      ONLY: nb_temp, nb_int_temp, nb_te, flag_diss
    USE mod_numerics_data,     ONLY: flux_splitter, flux_Jac_approx, flag_full_Impl, diff_flux_Jac
    USE mod_function_pointer

    IMPLICIT NONE

    ! Numerical flux 
    SELECT CASE(flux_splitter)

      ! van Leer's flux vector splitting 
      CASE('van_leer','van_Leer')
        conv_flux_1D => van_Leer_neq_1D             
 
      ! Steger-Warming's flux vector splitting 
      CASE ('steg_war','Steg_War')
        
        ! Eigensystem  
        SELECT CASE(nb_temp)

          ! 1 Temperature 
          CASE(1) 
            conv_flux_1D => Steger_Warming_neq1T_1D 

          ! N temperatures and eventually one for free electrons (Te)
          CASE DEFAULT 
            
            IF (nb_te.EQ.0) THEN
               conv_flux_1D => Steger_Warming_neqNT_1D 
            ELSE
               conv_flux_1D => Steger_Warming_neqNT_Te_1D 
            ENDIF

        END SELECT 

      ! Roe approximate Riemann solver
      CASE('roe','Roe')
        conv_flux_1D => roe_neq_1D

        ! Eigensystem  
        SELECT CASE(nb_temp)

          ! 1 Temperature 
          CASE(1) 
            eigsys_neq_1D  => eigensystem_neq1T_1D
            roe_avg_neq_1D => roe_avgState_neq1T_1D

          ! N temperatures and eventually one for free electrons (Te)
          CASE DEFAULT 
            
            IF (nb_te.EQ.0) THEN
               eigsys_neq_1D  => eigensystem_neqNT_1D
               roe_avg_neq_1D => roe_avgState_neqNT_1D
            ELSE
               eigsys_neq_1D  => eigensystem_neqNT_Te_1D
               roe_avg_neq_1D => roe_avgState_neqNT_Te_1D
            ENDIF

        END SELECT  

      ! HLLE's approximate Riemann solver 
      CASE('hlle','HLLE')
        conv_flux_1D => hlle_neq_1D

        ! Eigensystem  
        SELECT CASE(nb_temp)

          ! 1 Temperature 
          CASE(1) 
            roe_avg_neq_1D => roe_avgState_neq1T_1D

          ! N temperatures and eventually one for free electrons (Te)
          CASE DEFAULT 
            
            IF (nb_te.EQ.0) THEN
                roe_avg_neq_1D => roe_avgState_neqNT_1D
            ELSE
                roe_avg_neq_1D => roe_avgState_neqNT_Te_1D
            ENDIF

        END SELECT  

      CASE DEFAULT 

         WRITE(*,10)'In "set_flux_neq_1D.F90", error in numerical flux function selection ...'
         PRINT*
         STOP 

      END SELECT 

      ! Inviscid flux Jacobian
      SELECT CASE(nb_temp)

        ! 1 Temperature 
        CASE(1) 
          inv_flux_Jac_neq_1D => inviscid_flux_Jac_neq1T_1D

        ! N temperatures and eventually one for free electrons (Te)
        CASE DEFAULT 
            
          IF (nb_te.EQ.0) THEN
             inv_flux_Jac_neq_1D => inviscid_flux_Jac_neqNT_1D
          ELSE
             inv_flux_Jac_neq_1D => inviscid_flux_Jac_neqNT_Te_1D
          ENDIF

      END SELECT 

      ! Numerical flux Jacobian (only in case of a fully implicit time-integration method)
      IF (flag_full_Impl.EQV..TRUE.) THEN

        SELECT CASE(flux_Jac_approx) 

          CASE('Numerical','numerical')
             conv_flux_Jac_1D => conv_flux_num_Jac_1D

          CASE('Yoon_Jameson')
            conv_flux_Jac_1D => Yoon_Jameson_Jac_neq_1D

          CASE('pos_neg_split','Pos_Neg_split')
            SELECT CASE(nb_temp)

              ! 1 Temperature 
              CASE(1) 
              conv_flux_Jac_1D => Pos_Neg_split_Jac_neq1T_1D

              ! N temperatures and eventually one for free electrons (Te)
              CASE DEFAULT 

                IF (nb_te.EQ.0) THEN
                  conv_flux_Jac_1D => Pos_Neg_split_Jac_neqNT_1D
                ELSE
                  conv_flux_Jac_1D => Pos_Neg_split_Jac_neqNT_Te_1D 
                ENDIF

            END SELECT

          CASE DEFAULT
            WRITE(*,10)'In "set_num_flux_1D.F90", error in flux Jacobian approximation selection ...'
            PRINT*
            STOP 

        END SELECT

      ENDIF

      ! Diffusive flux (and related Jacobians)
      IF (flag_diss.EQV..TRUE.) THEN

         SELECT CASE(nb_temp)

           ! 1 Temperature 
           CASE(1) 
             diff_flux_1D => ns_flux_neq1T_1D
             IF (diff_flux_Jac=='numerical') THEN 
                diff_flux_1D_Jac => diff_flux_num_Jac_1D
             ELSE
                diff_flux_1D_Jac => ns_flux_neq1T_1D_Jac
             ENDIF

           ! N temperatures and eventually one for free electrons (Te)
           CASE DEFAULT 

             IF (nb_te.EQ.0) THEN
                diff_flux_1D => ns_flux_neqNT_1D
                IF (diff_flux_Jac=='numerical') THEN 
                   diff_flux_1D_Jac => diff_flux_num_Jac_1D
                ELSE
                   diff_flux_1D_Jac => ns_flux_neqNT_1D_Jac
                ENDIF
             ELSE
                diff_flux_1D => ns_flux_neqNT_Te_1D
                IF (diff_flux_Jac=='numerical') THEN 
                   diff_flux_1D_Jac => diff_flux_num_Jac_1D
                ELSE
                   diff_flux_1D_Jac => ns_flux_neqNT_Te_1D_Jac
                ENDIF
             ENDIF

         END SELECT

      ELSE 

         diff_flux_1D => null_ns_flux_1D
         diff_flux_1D_Jac => null_ns_flux_1D_Jac

      ENDIF

      ! Spectral radius 
      IF (flag_diss.EQV..TRUE.) THEN
         
         !get_spectral_radius_1D => visc_spectral_radius_1D 
         get_inv_spectral_radius_1D => inv_spectral_radius_1D 
         get_visc_spectral_radius_1D => visc_spectral_radius_neq_1D 

      ELSE

         !get_spectral_radius_1D => inv_spectral_radius_1D 
         get_inv_spectral_radius_1D => inv_spectral_radius_1D 
         get_visc_spectral_radius_1D => null_visc_spectral_radius_1D 

      ENDIF

10 FORMAT(A)

  END SUBROUTINE set_flux_neq_1D
!------------------------------------------------------------------------------!
