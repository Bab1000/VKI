!------------------------------------------------------------------------------!
!> This subroutine associates the pointer for the source term (and relative Jacobian) 
!! when solving the 1D stagnation line equation for nonequilibrium flows.
  SUBROUTINE set_source_term_neq_SL ()

    USE mod_general_data,       ONLY: nb_temp, nb_te
    USE mod_numerics_data,      ONLY: nb_source_terms, nb_inv_source_terms, source_terms, & 
                                    & flag_Jac, stag_line_geom, compute_source_term_Jac
    USE mod_general_data,       ONLY: flag_diss
    USE mod_function_pointer 

    IMPLICIT NONE

    INTEGER :: i

    ! Number of inviscid source terms
    IF (flag_diss.EQV..TRUE.) THEN
       nb_inv_source_terms = nb_source_terms - 1
    ELSE
       nb_inv_source_terms = nb_source_terms
    ENDIF

    ! Inviscid source terms
    ! Selection according to the geometry speciefied in input
    SELECT CASE(stag_line_geom)

      ! Cylinder case
      CASE(0)

        ! Source term(s) only
        IF (flag_Jac.EQV..FALSE.) THEN

           ALLOCATE(get_inv_source_term_SL(nb_inv_source_terms))

           DO i = 1,nb_inv_source_terms 

              SELECT CASE(source_terms(i))       

                CASE('kinetics','Kinetics')
                  IF (nb_te.EQ.1) THEN 
                    get_inv_source_term_SL(i)%source => source_term_neq_el_kinetics_SL
                  ELSE
                    get_inv_source_term_SL(i)%source => source_term_neq_kinetics_SL
                  ENDIF

                CASE('inv_stag','Inv_stag')
                  get_inv_source_term_SL(i)%source => source_term_inv_neq_1D_SL_cyl

                CASE DEFAULT
                  WRITE(*,10)'In "set_source_term_neq_SL.F90", source term not implemented yet...'
                  WRITE(*,10)source_terms(i)
                  PRINT*
                  STOP

              END SELECT 

           ENDDO

        ! Source term(s) and source term Jacobian(s)
        ELSE 

          ALLOCATE(get_inv_source_term_Jac_SL(nb_source_terms))
          ALLOCATE(get_inv_source_term_SL(nb_source_terms))

          DO i = 1,nb_inv_source_terms 

             SELECT CASE(source_terms(i))       

               CASE('kinetics','Kinetics')
                  IF (compute_source_term_Jac(i)=='numerical') THEN
                    get_inv_source_term_Jac_SL(i)%source_Jac => source_term_inv_num_Jac_SL
                    IF (nb_te.EQ.1) THEN 
                      get_inv_source_term_SL(i)%source => source_term_neq_el_kinetics_SL
                    ELSE
                      get_inv_source_term_SL(i)%source => source_term_neq_kinetics_SL
                    ENDIF 
                  ELSEIF (compute_source_term_Jac(i)=='analytical') THEN
                    IF (nb_te.EQ.1) THEN 
                      get_inv_source_term_Jac_SL(i)%source_Jac => source_term_neq_el_kinetics_Jac_SL
                    ELSE
                      get_inv_source_term_Jac_SL(i)%source_Jac => source_term_neq_kinetics_Jac_SL
                    ENDIF
                  ELSE 
                    WRITE(*,10)'In "set_source_term_neq_SL.F90", error in setting the Jacobian option for source term'
                    WRITE(*,10)source_terms(i)
                    PRINT*
                    STOP
                  ENDIF

               CASE('inv_stag','Inv_stag')
                 IF (compute_source_term_Jac(i)=='numerical') THEN
                    get_inv_source_term_Jac_SL(i)%source_Jac => source_term_inv_num_Jac_SL
                    get_inv_source_term_SL(i)%source => source_term_inv_neq_1D_SL_cyl
                 ELSEIF (compute_source_term_Jac(i)=='analytical') THEN 
                    get_inv_source_term_Jac_SL(i)%source_Jac => source_term_inv_neq_1D_SL_cyl_Jac
                 ELSE 
                    WRITE(*,10)'In "set_source_term_neq_SL.F90", error in setting the Jacobian option for source term'
                    WRITE(*,10)source_terms(i)
                    PRINT*
                    STOP
                 ENDIF

               CASE DEFAULT
                 WRITE(*,10)'In "set_source_term_neq_SL.F90", source term not implemented yet...'
                 WRITE(*,10)source_terms(i)
                 PRINT*
                 STOP

             END SELECT 

          ENDDO

        ENDIF

     ! Sphere case
     CASE(1)

       ! Source term(s) only
        IF (flag_Jac.EQV..FALSE.) THEN

           ALLOCATE(get_inv_source_term_SL(nb_inv_source_terms))

           DO i = 1,nb_inv_source_terms 

              SELECT CASE(source_terms(i))       

                CASE('kinetics','Kinetics')
                  IF (nb_te.EQ.1) THEN 
                    get_inv_source_term_SL(i)%source => source_term_neq_el_kinetics_SL
                  ELSE
                    get_inv_source_term_SL(i)%source => source_term_neq_kinetics_SL
                  ENDIF

                CASE('inv_stag','Inv_stag')
                  get_inv_source_term_SL(i)%source => source_term_inv_neq_1D_SL_sph

                CASE DEFAULT
                  WRITE(*,10)'In "set_source_term_neq_SL.F90", source term not implemented yet...'
                  WRITE(*,10)source_terms(i)
                  PRINT*
                  STOP

              END SELECT 

           ENDDO

        ! Source term(s) and source term Jacobian(s)
        ELSE 

           ALLOCATE(get_inv_source_term_Jac_SL(nb_source_terms))
           ALLOCATE(get_inv_source_term_SL(nb_source_terms))

          DO i = 1,nb_inv_source_terms 

             SELECT CASE(source_terms(i))       

               CASE('kinetics','Kinetics')
                  IF (compute_source_term_Jac(i)=='numerical') THEN
                    get_inv_source_term_Jac_SL(i)%source_Jac => source_term_inv_num_Jac_SL
                    IF (nb_te.EQ.1) THEN 
                      get_inv_source_term_SL(i)%source => source_term_neq_el_kinetics_SL
                    ELSE
                      get_inv_source_term_SL(i)%source => source_term_neq_kinetics_SL
                    ENDIF 
                  ELSEIF (compute_source_term_Jac(i)=='analytical') THEN
                    IF (nb_te.EQ.1) THEN 
                      get_inv_source_term_Jac_SL(i)%source_Jac => source_term_neq_el_kinetics_Jac_SL
                    ELSE
                      get_inv_source_term_Jac_SL(i)%source_Jac => source_term_neq_kinetics_Jac_SL
                    ENDIF
                  ELSE 
                    WRITE(*,10)'In "set_source_term_neq_SL.F90", error in setting the Jacobian option for source term'
                    WRITE(*,10)source_terms(i)
                    PRINT*
                    STOP
                  ENDIF

               CASE('inv_stag','Inv_stag')
                 IF (compute_source_term_Jac(i)=='numerical') THEN
                    get_inv_source_term_Jac_SL(i)%source_Jac => source_term_inv_num_Jac_SL
                    get_inv_source_term_SL(i)%source => source_term_inv_neq_1D_SL_sph
                 ELSEIF (compute_source_term_Jac(i)=='analytical') THEN 
                    get_inv_source_term_Jac_SL(i)%source_Jac => source_term_inv_neq_1D_SL_sph_Jac
                 ELSE 
                    WRITE(*,10)'In "set_source_term_neq_SL.F90", error in setting the Jacobian option for source term'
                    WRITE(*,10)source_terms(i)
                    PRINT*
                    STOP
                 ENDIF

               CASE DEFAULT
                 WRITE(*,10)'In "set_source_term_neq_SL.F90", source term not implemented yet...'
                 WRITE(*,10)source_terms(i)
                 PRINT*
                 STOP

             END SELECT 

          ENDDO

        ENDIF

     CASE DEFAULT
       WRITE(*,10)'In "set_source_term_neq_SL.F90", error in geometry parameter selection...'
       WRITE(*,20)stag_line_geom
       PRINT*
       STOP

   END SELECT

   ! Diffusive source term (and related Jacobians)
   IF (flag_diss.EQV..FALSE.) THEN

      get_diff_source_term_SL => null_diff_source_term_SL 
      get_diff_source_term_SL_Jac => null_diff_source_term_SL_Jac 

   ELSE

      ! Geometry selection
      SELECT CASE(stag_line_geom) 
 
        ! Cylinder
        CASE(0)

          ! 1 temperature 
          IF (nb_temp.EQ.1) THEN

             get_diff_source_term_SL => source_term_diff_neq1T_1D_SL_cyl

          ! N temperatures and eventually one for free electrons (Te)
          ELSE 

             IF (nb_te.EQ.0) THEN
                get_diff_source_term_SL => source_term_diff_neqNT_1D_SL_cyl
             ELSE
                get_diff_source_term_SL => source_term_diff_neqNT_Te_1D_SL_cyl
             ENDIF

          ENDIF

          IF (compute_source_term_Jac(nb_source_terms).EQ.'numerical') THEN
             get_diff_source_term_SL_Jac => source_term_diff_num_Jac_SL
          ELSEIF (compute_source_term_Jac(nb_source_terms).EQ.'analytical') THEN

             ! 1 temperature 
             IF (nb_temp.EQ.1) THEN

                get_diff_source_term_SL_Jac => source_term_diff_neq1T_1D_SL_cyl_Jac

             ! N temperatures and eventually one for free electrons (Te)
             ELSE 

               IF (nb_te.EQ.0) THEN
                  get_diff_source_term_SL_Jac => source_term_diff_neqNT_1D_SL_cyl_Jac
               ELSE
                  get_diff_source_term_SL_Jac => source_term_diff_neqNT_Te_1D_SL_cyl_Jac
               ENDIF

             ENDIF

          ENDIF

        ! Sphere
        CASE(1)

          ! 1 temperature 
          IF (nb_temp.EQ.1) THEN

             get_diff_source_term_SL => source_term_diff_neq1T_1D_SL_sph

          ! N temperatures and eventually one for free electrons (Te)
          ELSE 

             IF (nb_te.EQ.0) THEN
                get_diff_source_term_SL => source_term_diff_neqNT_1D_SL_sph
             ELSE
                get_diff_source_term_SL => source_term_diff_neqNT_Te_1D_SL_sph
             ENDIF

          ENDIF

          IF (compute_source_term_Jac(nb_source_terms).EQ.'numerical') THEN
             get_diff_source_term_SL_Jac => source_term_diff_num_Jac_SL
          ELSEIF (compute_source_term_Jac(nb_source_terms).EQ.'analytical') THEN

             ! 1 temperature 
             IF (nb_temp.EQ.1) THEN

                get_diff_source_term_SL_Jac => source_term_diff_neq1T_1D_SL_sph_Jac

             ! N temperatures and eventually one for free electrons (Te)
             ELSE 

               IF (nb_te.EQ.0) THEN
                  get_diff_source_term_SL_Jac => source_term_diff_neqNT_1D_SL_sph_Jac
               ELSE
                  get_diff_source_term_SL_Jac => source_term_diff_neqNT_Te_1D_SL_sph_Jac
               ENDIF

             ENDIF

          ENDIF

      END SELECT

   ENDIF


10 FORMAT(A)
20 FORMAT(I4)

  END SUBROUTINE set_source_term_neq_SL
!------------------------------------------------------------------------------!
