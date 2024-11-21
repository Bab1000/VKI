!------------------------------------------------------------------------------!
!> This subroutine writes the solution flowfield output files.   
  SUBROUTINE write_sol ()

    USE mod_general_data,         ONLY: nb_dim, model_name
    USE mod_numerics_data,        ONLY: flag_stag_line
    USE mod_general_ablation,     ONLY: Flag_ablation

    IMPLICIT NONE

    EXTERNAL write_sol_1D, write_sol_1D_SL
    EXTERNAL write_sol_neq_1D, write_sol_neq_1D_SL
    EXTERNAL write_sol_2D, write_sol_neq_2D
    EXTERNAL write_surface_data_1D_SL
    EXTERNAL write_transp_flux_1D, write_transp_flux_1D_SL
    EXTERNAL write_transp_flux_neq_1D, write_transp_flux_neq_1D_SL
   
   

    ! Physical model selection 
    SELECT CASE(model_name)

      ! Calorically perfect gas flow
      CASE('pg') 

        SELECT CASE (nb_dim)

          ! 1D flow
          CASE (1)

            CALL write_sol_1D ()       
            CALL write_transp_flux_1D ()
          ! 2D flow
          CASE (2)

            ! Stagnation line flow
            IF (flag_stag_line.EQV..TRUE.) THEN

               CALL write_sol_1D_SL () 
               CALL write_transp_flux_1D_SL ()
               CALL write_surface_data_1D_SL()

            ELSE
 
               CALL write_sol_2D () 

            ENDIF

        END SELECT
      
      ! Nonequilibrium flow
      CASE('neq')

        SELECT CASE (nb_dim)

          ! 1D flow
          CASE (1)
           
            CALL write_sol_neq_1D ()      
            CALL write_transp_flux_neq_1D () 

          ! 2D flow
          CASE (2)

            ! Stagnation line flow
            IF (flag_stag_line.EQV..TRUE.) THEN

               CALL write_surface_data_1D_SL()
               CALL write_sol_neq_1D_SL () 
               CALL write_transp_flux_neq_1D_SL ()
               
              IF (Flag_ablation.EQV..TRUE.) THEN
               
!               CALL write_ablation_solution ()
               

              ENDIF
            ELSE

               CALL write_sol_neq_2D () 

            ENDIF

        END SELECT

    END SELECT 

  END SUBROUTINE write_sol
!------------------------------------------------------------------------------!  
