!------------------------------------------------------------------------------!
!> This subroutine associates the pointer for the source term (and relative Jacobian) 
!! (calorically perfect gas flows).
  SUBROUTINE set_source_term ()

    USE mod_numerics_data,      ONLY: nb_source_terms, source_terms, flag_Jac, compute_source_term_Jac
    USE mod_function_pointer 

    IMPLICIT NONE

    INTEGER :: i

    ! Source term(s) only
    IF (flag_Jac.EQV..FALSE.) THEN

       IF (nb_source_terms.GE.1) THEN

          ALLOCATE(get_source_term(nb_source_terms))

          DO i = 1,nb_source_terms 

             SELECT CASE(source_terms(i))       

               CASE('quasi1D','quasi1d')
                 get_source_term(i)%source => source_term_quasi1D

               CASE DEFAULT
                 WRITE(*,10)'In "set_source_term.F90", source term not implemented yet...'
                 WRITE(*,10)source_terms(i)
                 PRINT*
                 STOP

             END SELECT 

          ENDDO

       ENDIF

    ! Source term(s) and source term Jacobian(s)
    ELSE 

        IF (nb_source_terms.GE.1) THEN

          ! Analytical Jacobian 
          ALLOCATE(get_source_term_Jac(nb_source_terms))
          ALLOCATE(get_source_term(nb_source_terms))

             DO i = 1,nb_source_terms 

                SELECT CASE(source_terms(i))       

                 CASE('quasi1D','quasi1d')
                   IF (compute_source_term_Jac(i)=='numerical') THEN
                      get_source_term_Jac(i)%source_Jac => source_term_num_Jac
                      get_source_term(i)%source => source_term_quasi1D
                   ELSEIF (compute_source_term_Jac(i)=='analytical') THEN 
                      get_source_term_Jac(i)%source_Jac => source_term_quasi1D_Jac
                   ELSE 
                      WRITE(*,10)'In "set_source_term.F90", error in setting the Jacobian option for source term'
                      WRITE(*,10)source_terms(i)
                      PRINT*
                      ENDIF

                 CASE DEFAULT

                   WRITE(*,10)'In "set_source_term.F90", source term not implemented yet...'
                   WRITE(*,10)source_terms(i)
                   PRINT*
                   STOP

                END SELECT 

             ENDDO

        ENDIF

    ENDIF

10 FORMAT(A)

  END SUBROUTINE set_source_term
!------------------------------------------------------------------------------!
