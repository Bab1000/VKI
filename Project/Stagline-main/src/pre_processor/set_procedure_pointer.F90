!------------------------------------------------------------------------------!
!> This subroutine selects the subroutine names from the specifications
!! provided in the input file. It is used in order to speed up computations.
!! In order to perform this task pointers for functions/subroutines are used. 
  SUBROUTINE set_procedure_pointer ()

    USE mod_general_data,          ONLY: model_name, nb_dim, flag_diss 
    USE mod_numerics_data,         ONLY: flag_stag_line

    IMPLICIT NONE

    ! Physical model selection
    SELECT CASE (model_name)

      ! Calorically perfect gas flow
      CASE ('pg')

        SELECT CASE (nb_dim)

          ! 1D flow
          CASE (1)
            CALL set_pointer_1D()

          ! 2D flow
          CASE (2) 

            ! Stagnation line flow
            IF (flag_stag_line.EQV..TRUE.) THEN
               CALL set_pointer_1D_SL()
            ELSE 
               WRITE(*,10)'In "set_procedure_pointer.F90", error in physical model selection ...'
               PRINT*
               STOP
            ENDIF 

        END SELECT 

      ! Nonequilibrium flow
      CASE ('neq')

        SELECT CASE (nb_dim)

          ! 1D flow
          CASE (1)
            CALL set_pointer_neq_1D()
            
          ! 2D flow
          CASE (2) 
            IF (flag_stag_line.EQV..TRUE.) THEN
               CALL set_pointer_neq_1D_SL()
            ELSE 
               WRITE(*,10)'In "set_procedure_pointer.F90", error in physical model selection ...'
               PRINT*
               STOP
            ENDIF 

        END SELECT 

      CASE DEFAULT 
        WRITE(*,10)'In "set_procedure_pointer.F90", error in physical model selection...'
        PRINT*
        STOP 

    END SELECT

10  FORMAT(A)

    ! Subroutines for the various cases.
    CONTAINS 

      !----------------------------------------------------!
      !> This subroutine associates function/subroutine pointers for 1D calorically perfect gas flows.
      SUBROUTINE set_pointer_1D ()

        USE mod_numerics_data,        ONLY: poly_rec 
        USE mod_function_pointer

        IMPLICIT NONE

        EXTERNAL set_bc_1D, set_flux_1D, set_time_step_1D
        EXTERNAL set_time_integration_1D
        EXTERNAL set_limiter, set_source_term

        ! Compute physical properties from conservative variables 
        IF (flag_diss.EQV..FALSE.) THEN
           get_phys_from_cons => cons_to_phys_1D_Eu
        ELSE 
           get_phys_from_cons => cons_to_phys_1D_Ns 
        ENDIF

        ! Compute physical properties and conservative variables from primitive
        ! variables
        IF (flag_diss.EQV..FALSE.) THEN
           get_cons_phys_from_prim => prim_to_cons_phys_1D_Eu
        ELSE 
           get_cons_phys_from_prim => prim_to_cons_phys_1D_Ns  
        ENDIF

        ! Compute pre-conditioning velocity 
        IF (flag_diss.EQV..FALSE.) THEN
           get_prec_vel_1D => inv_prec_vel_1D
        ELSE 
           get_prec_vel_1D => visc_prec_vel_1D  
        ENDIF

        ! Conversion from primitive to conservative variables
        get_cons_from_prim => prim_to_cons_1D  

        ! Conversion from conservative to primitive variables
        get_prim_from_cons => cons_to_prim_1D  

        ! Convective flux, diffusive flux (and related Jacobians)  
        CALL set_flux_1D ()

        IF (poly_rec.EQ.'linear') THEN

           ! Polynomial re-construction
           rec_1D => poly_rec_1D
        
           ! Limiter for polynomial re-construction 
           CALL set_limiter ()

        ELSEIF (poly_rec.EQ.'constant') THEN

           ! Null polynomial re-construction
           rec_1D => null_poly_rec_1D

        ELSE 

           WRITE(*,10)'in "set_pointer_1D", error in polynomial reconstruction selection...'
           PRINT*
           STOP

        ENDIF

        ! Boundary conditions
        CALL set_bc_1D ()

        ! Source term(s) and Jacobian(s)
        CALL set_source_term ()
        
        ! Time-step
        CALL set_time_step_1D () 

        ! Transport coefficients
        CALL set_transp_coeff ()

        ! Set time-integration method
        CALL set_time_integration_1D()

10    FORMAT(A)

      END SUBROUTINE set_pointer_1D

      !----------------------------------------------------!
      !> This subroutine associates function/subroutine pointers for 1D stagnation line calorically perfect gas flows.
      SUBROUTINE set_pointer_1D_SL ()

        USE mod_numerics_data,        ONLY: poly_rec, flag_metrics
        USE mod_function_pointer

        IMPLICIT NONE

        EXTERNAL set_bc_1D_SL, set_flux_1D_SL, set_time_step_1D_SL
        EXTERNAL set_limiter, set_source_term_SL
        EXTERNAL set_stress_tensor_1D_SL

        ! Compute physical properties from conservative variables 
        IF (flag_diss.EQV..FALSE.) THEN
           get_phys_from_cons => cons_to_phys_1D_SL_Eu
        ELSE 
           get_phys_from_cons => cons_to_phys_1D_SL_Ns 
        ENDIF

        ! Compute physical properties and conservative variables from primitive
        ! variables
        IF (flag_diss.EQV..FALSE.) THEN
           get_cons_phys_from_prim => prim_to_cons_phys_1D_SL_Eu
        ELSE 
           get_cons_phys_from_prim => prim_to_cons_phys_1D_SL_Ns  
        ENDIF

        ! Compute pre-conditioning velocity 
        IF (flag_diss.EQV..FALSE.) THEN
           get_prec_vel_1D_SL => inv_prec_vel_1D_SL
        ELSE 
           get_prec_vel_1D_SL => visc_prec_vel_1D_SL  
        ENDIF 

        ! Conversion from primitive to conservative variables
        get_cons_from_prim => prim_to_cons_1D_SL  

        ! Conversion from conservative to primitive variables
        get_prim_from_cons => cons_to_prim_1D_SL  

        ! Convective flux, diffusive flux (and related Jacobians)  
        CALL set_flux_1D_SL ()

        IF (poly_rec.EQ.'linear') THEN
           IF (flag_metrics.EQV..FALSE.) THEN
              ! Polynomial re-construction w/o metrics
              rec_1D_SL => poly_rec_1D_SL
           ELSE
              ! Polynomial re-construction w/ metrics (A.Turchi)
              rec_1D_SL => poly_rec_1D_SL_metr
           ENDIF
        
           ! Limiter for polynomial re-construction 
           CALL set_limiter ()

        ELSEIF (poly_rec.EQ.'constant') THEN

           ! Null polynomial re-construction
           rec_1D_SL => null_poly_rec_1D_SL

        ELSE 

           WRITE(*,10)'in "set_pointer_1D_SL", error in polynomial reconstruction selection...'
           PRINT*
           STOP

        ENDIF

        ! Boundary conditions
        CALL set_bc_1D_SL ()

        ! Source term(s) and Jacobian(s)
        CALL set_source_term_SL ()
        
        ! Time-step
        CALL set_time_step_1D_SL () 

        ! Transport coefficients
        CALL set_transp_coeff ()

        ! Stress tensor
        CALL set_stress_tensor_1D_SL ()

        ! Set time-integration method
        CALL set_time_integration_1D_SL()

10    FORMAT(A)

      END SUBROUTINE set_pointer_1D_SL

      !----------------------------------------------------!
      !> This subroutine defines function/subroutine pointers for 1D nonequilibrium flows.
      SUBROUTINE set_pointer_neq_1D ()

        USE mod_numerics_data,        ONLY: poly_rec 
        USE mod_function_pointer 

        IMPLICIT NONE
       
        EXTERNAL set_bc_neq_1D, set_flux_neq_1D, set_time_step_1D
        EXTERNAL set_time_integration_1D
        EXTERNAL set_limiter, set_source_term_neq

        ! Compute physical properties from conservative variables 
        IF (flag_diss.EQV..FALSE.) THEN
           get_phys_from_cons => cons_to_phys_neq_1D_Eu
        ELSE 
           get_phys_from_cons => cons_to_phys_neq_1D_Ns 
        ENDIF

        ! Compute physical properties and conservative variables from primitive
        ! variables
        IF (flag_diss.EQV..FALSE.) THEN
           get_cons_phys_from_prim => prim_to_cons_phys_neq_1D_Eu
        ELSE 
           get_cons_phys_from_prim => prim_to_cons_phys_neq_1D_Ns  
        ENDIF

        ! Compute pre-conditioning velocity 
        IF (flag_diss.EQV..FALSE.) THEN
           get_prec_vel_1D => inv_prec_vel_1D
        ELSE 
           get_prec_vel_1D => visc_prec_vel_1D  
        ENDIF

        ! Conversion from primitive to conservative variables
        get_cons_from_prim => prim_to_cons_neq_1D  

        ! Conversion from conservative to primitive variables
        get_prim_from_cons => cons_to_prim_neq_1D  

        ! Convective flux, diffusive flux (and related Jacobians) and eigensystem 
        CALL set_flux_neq_1D ()

        IF (poly_rec.EQ.'linear') THEN

           ! Polynomial re-construction
           rec_1D => poly_rec_neq_1D
        
           ! Limiter for polynomial re-construction 
           CALL set_limiter ()

        ELSEIF (poly_rec.EQ.'constant') THEN
          
           ! Null polynomial re-construction
           rec_1D => null_poly_rec_1D
          
        ELSE 

           WRITE(*,10)'in "set_pointer_neq_1D", error in polynomial reconstruction selection...'
           PRINT*
           STOP

        ENDIF

        ! Boundary conditions
        CALL set_bc_neq_1D ()

        ! Source term(s) and Jacobian(s)
        CALL set_source_term_neq ()
        
        ! Time-step
        CALL set_time_step_1D () 

        ! Set time-integration method
        CALL set_time_integration_1D()

10    FORMAT(A)

      END SUBROUTINE set_pointer_neq_1D 

      !----------------------------------------------------!
      !> This subroutine defines function/subroutine pointers for 1D stagnation line nonequilibrium flows.
      SUBROUTINE set_pointer_neq_1D_SL ()

        USE mod_numerics_data,        ONLY: poly_rec, flag_metrics
        USE mod_function_pointer 

        IMPLICIT NONE
       
        EXTERNAL set_bc_neq_1D_SL, set_flux_neq_1D_SL
        EXTERNAL set_time_step_1D_SL
        EXTERNAL set_time_integration_1D_SL
        EXTERNAL set_limiter
        EXTERNAL set_source_term_neq_SL
        EXTERNAL set_stress_tensor_1D_SL

        ! Compute physical properties from conservative variables 
        IF (flag_diss.EQV..FALSE.) THEN
           get_phys_from_cons => cons_to_phys_neq_1D_SL_Eu
        ELSE 
           get_phys_from_cons => cons_to_phys_neq_1D_SL_Ns 
        ENDIF

        ! Compute physical properties and conservative variables from primitive
        ! variables
        IF (flag_diss.EQV..FALSE.) THEN
           get_cons_phys_from_prim => prim_to_cons_phys_neq_1D_SL_Eu
        ELSE 
           get_cons_phys_from_prim => prim_to_cons_phys_neq_1D_SL_Ns  
        ENDIF

        ! Compute pre-conditioning velocity 
        IF (flag_diss.EQV..FALSE.) THEN
           get_prec_vel_1D_SL => inv_prec_vel_1D_SL
        ELSE 
           get_prec_vel_1D_SL => visc_prec_vel_1D_SL  
        ENDIF

        ! Conversion from primitive to conservative variables
        get_cons_from_prim => prim_to_cons_neq_1D_SL  

        ! Conversion from conservative to primitive variables
        get_prim_from_cons => cons_to_prim_neq_1D_SL

        ! Convective flux, diffusive flux (and related Jacobians) and eigensystem 
        CALL set_flux_neq_1D_SL ()

        IF (poly_rec.EQ.'linear') THEN
           IF (flag_metrics.EQV..FALSE.) THEN
              ! Polynomial re-construction w/o metrics
              rec_1D_SL => poly_rec_neq_1D_SL 
           ELSE
              ! Polynomial re-construction w/ metrics (A.Turchi)
              rec_1D_SL => poly_rec_neq_1D_SL_metr
           ENDIF

           ! Limiter for polynomial re-construction 
           CALL set_limiter ()

        ELSEIF (poly_rec.EQ.'constant') THEN
          
           ! Null polynomial re-construction
           rec_1D_SL => null_poly_rec_1D_SL
          
        ELSE 

           WRITE(*,10)'in "set_pointer_neq_1D", error in polynomial reconstruction selection...'
           PRINT*
           STOP

        ENDIF

        ! Boundary conditions
        CALL set_bc_neq_1D_SL ()

        ! Source term(s) and Jacobian(s)
        CALL set_source_term_neq_SL ()
        
        ! Time-step
        CALL set_time_step_1D_SL () 

        ! Stress tensor
        CALL set_stress_tensor_1D_SL ()

        ! Set time-integration method
        CALL set_time_integration_1D_SL()

10    FORMAT(A)  

      END SUBROUTINE set_pointer_neq_1D_SL

  END SUBROUTINE set_procedure_pointer
!------------------------------------------------------------------------------!
