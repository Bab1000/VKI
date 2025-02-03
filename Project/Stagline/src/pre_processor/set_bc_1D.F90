!------------------------------------------------------------------------------!
!> This subroutine associates pointers for boundary conditions (1D calorically perfect gas flows)
  SUBROUTINE set_bc_1D ()

    USE mod_numerics_data,          ONLY: poly_rec, flag_full_Impl
    USE mod_domain_boundary,        ONLY: bc
    USE mod_function_pointer

    IMPLICIT NONE

    INTEGER :: i

    ! Explicit boundary conditions 
    IF (flag_full_Impl.EQV..FALSE.) THEN 

       ! First order accuracy in space
       IF (poly_rec.EQ.'constant') THEN

          ALLOCATE(get_ghost_state_Expl_1D_1st(2))

          ! Loop over domain boundaries
          DO i = 1,2

             SELECT CASE (bc(i))

               CASE ('sup_out','trans')
                 get_ghost_state_Expl_1D_1st(i)%bc_procedure => sup_out_Expl

               CASE ('sub_in_rhoin_pin','sub_in_pin_rhoin')
                 get_ghost_state_Expl_1D_1st(i)%bc_procedure => sub_in_rhoin_pin_1DExpl

               CASE ('sup_in')
                 get_ghost_state_Expl_1D_1st(i)%bc_procedure => sup_in_1DExpl

               CASE ('sub_out_pout')
                  get_ghost_state_Expl_1D_1st(i)%bc_procedure => sub_out_pout_1DExpl 

              CASE DEFAULT 
                WRITE(*,10)'In "set_bc_1D.F90", bc not implemented yet...'
                WRITE(*,10)bc(i)
                PRINT*
                STOP 

             END SELECT 

          ENDDO

       ! Second order accuracy in space
       ELSE 

          ALLOCATE(get_ghost_state_Expl_1D_2nd(2))

          ! Loop over domain boundaries
          DO i = 1,2

             SELECT CASE (bc(i))

               CASE ('sup_out','trans')
                 get_ghost_state_Expl_1D_2nd(i)%bc_procedure => sup_out_1DExpl_2nd

               CASE ('sub_in_rhoin_pin','sub_in_pin_rhoin')
                 get_ghost_state_Expl_1D_2nd(i)%bc_procedure => sub_in_rhoin_pin_1DExpl_2nd

               CASE ('sup_in')
                 get_ghost_state_Expl_1D_2nd(i)%bc_procedure => sup_in_1DExpl_2nd

               CASE ('sub_out_pout')
                  get_ghost_state_Expl_1D_2nd(i)%bc_procedure => sub_out_pout_1DExpl_2nd 

              CASE DEFAULT 
                WRITE(*,10)'In "set_bc_1D.F90", bc not implemented yet...'
                WRITE(*,10)bc(i)
                PRINT*
                STOP 

             END SELECT 

          ENDDO

       ENDIF

    ! Implicit boundary conditions
    ELSE 

        ! First order accuracy in space
       IF (poly_rec.EQ.'constant') THEN

          ALLOCATE(get_ghost_state_Impl_1D_1st(2))

          ! Loop over domain boundaries
          DO i = 1,2

             SELECT CASE (bc(i))

               CASE ('sup_out','trans')
                 get_ghost_state_Impl_1D_1st(i)%bc_procedure => sup_out_1DImpl

               CASE ('sub_in_rhoin_pin','sub_in_pin_rhoin')
                 get_ghost_state_Impl_1D_1st(i)%bc_procedure => sub_in_rhoin_pin_1DImpl

               CASE ('sup_in')
                 get_ghost_state_Impl_1D_1st(i)%bc_procedure => sup_in_1DImpl

               CASE ('sub_out_pout')
                  get_ghost_state_Impl_1D_1st(i)%bc_procedure => sub_out_pout_1DImpl 

               CASE DEFAULT 
                 WRITE(*,10)'In "set_bc_1D.F90", bc not implemented yet...'
                 WRITE(*,10)bc(i)
                 PRINT*
                 STOP 

             END SELECT 

          ENDDO

       ! Second order accuracy in space
       ELSE 

          ALLOCATE(get_ghost_state_Impl_1D_2nd(2))

          ! Loop over domain boundaries
          DO i = 1,2

             SELECT CASE (bc(i))

               CASE ('sup_out','trans')
                 get_ghost_state_Impl_1D_2nd(i)%bc_procedure => sup_out_1DImpl_2nd

               CASE ('sub_in_rhoin_pin','sub_in_pin_rhoin')
                 get_ghost_state_Impl_1D_2nd(i)%bc_procedure => sub_in_rhoin_pin_1DImpl_2nd

               CASE ('sup_in')
                 get_ghost_state_Impl_1D_2nd(i)%bc_procedure => sup_in_1DImpl_2nd

               CASE ('sub_out_pout')
                  get_ghost_state_Impl_1D_2nd(i)%bc_procedure => sub_out_pout_1DImpl_2nd 

              CASE DEFAULT 
                WRITE(*,10)'In "set_bc_1D.F90", bc not implemented yet...'
                WRITE(*,10)bc(i)
                PRINT*
                STOP 

             END SELECT 

          ENDDO

       ENDIF

   ENDIF

10 FORMAT(A)

  END SUBROUTINE set_bc_1D
!------------------------------------------------------------------------------!
