!------------------------------------------------------------------------------!
! This subroutine writes the transport fluxes for viscous calculations
  SUBROUTINE write_transport_fluxes()

    USE mod_general_data,            ONLY: nb_dim, flag_diss, model_name
    USE mod_numerics_data,           ONLY: flag_stag_line

    IMPLICIT NONE

    EXTERNAL write_transp_flux_1D, write_transp_flux_1D_SL
    EXTERNAL write_transp_flux_neq_1D, write_transp_flux_neq_1D_SL
    EXTERNAL write_transp_flux_2D, write_transp_flux_neq_2D

    IF (flag_diss.EQV..TRUE.) THEN

       ! Physical model selection 
       SELECT CASE(model_name)

         ! Calorically perfect gas flow
         CASE('pg') 

           SELECT CASE (nb_dim)

             ! 1D flow
             CASE (1)

               CALL write_transp_flux_1D ()       

             ! 2D flow
             CASE (2)

               ! Stagnation line flow
               IF (flag_stag_line.EQV..TRUE.) THEN

                 CALL write_transp_flux_1D_SL () 

               ELSE
 
                 CALL write_transp_flux_2D () 

              ENDIF

           END SELECT
      
         ! Nonequilibrium flow
         CASE('neq')

            SELECT CASE (nb_dim)

              ! 1D flow
              CASE (1)
           
                CALL write_transp_flux_neq_1D ()       

              ! 2D flow
              CASE (2)

                ! Stagnation line flow
                IF (flag_stag_line.EQV..TRUE.) THEN

                   CALL write_transp_flux_neq_1D_SL () 

                ELSE

                  CALL write_transp_flux_neq_2D () 

                ENDIF

            END SELECT

       END SELECT 

    ENDIF

  END SUBROUTINE write_transport_fluxes
!------------------------------------------------------------------------------!
