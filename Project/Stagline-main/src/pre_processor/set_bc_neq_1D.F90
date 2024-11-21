!------------------------------------------------------------------------------!
!> This subroutine associates pointers for boundary conditions (1D nonequilibrium flows)
  SUBROUTINE set_bc_neq_1D ()

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

               CASE ('sup_in')
                 get_ghost_state_Expl_1D_1st(i)%bc_procedure => sup_in_neq_1DExpl

               CASE ('sup_out','trans')
                 get_ghost_state_Expl_1D_1st(i)%bc_procedure => sup_out_Expl

               CASE ('sub_in_rhoiin_Tin','sub_in_rhoiin_Tin_Tvin', &
                   & 'sub_in_rhoiin_Tin_Trotin',                   & 
                   & 'sub_in_rhoiin_Tin_Trotin_Tvin',              & 
                   & 'sub_in_rhoiin_Tin_Trotin_Tvin_Tein',         &
                   & 'sub_in_rhoiin_Tin_Tein',                     &
                   & 'sub_in_rhoiin_Tein')
                 get_ghost_state_Expl_1D_1st(i)%bc_procedure => sub_in_neq_rhoiin_T_1DExpl
              
               CASE ('sub_out_Tout')
                 get_ghost_state_Expl_1D_1st(i)%bc_procedure =>  sub_out_neq_Tout_1DExpl 

               CASE DEFAULT 
                 WRITE(*,10)'In "set_bc_neq_1D.F90", bc not implemented yet...'
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

               CASE ('sup_in')
                 get_ghost_state_Expl_1D_2nd(i)%bc_procedure => sup_in_neq_1DExpl_2nd

               CASE ('sup_out','trans')
                 get_ghost_state_Expl_1D_2nd(i)%bc_procedure => sup_out_neq_1DExpl_2nd 

               CASE ('sub_in_rhoiin_Tin','sub_in_rhoiin_Tin_Tvin',  & 
                   & 'sub_in_rhoiin_Tin_Trotin',                    &  
                   & 'sub_in_rhoiin_Tin_Trotin_Tvin',               & 
                   & 'sub_in_rhoiin_Tin_Trotin_Tvin_Tein',          & 
                   & 'sub_in_rhoiin_Tin_Tein',                      &
                   & 'sub_in_rhoiin_Tein' )
                 get_ghost_state_Expl_1D_2nd(i)%bc_procedure => sub_in_neq_rhoiin_T_1DExpl_2nd 

               CASE DEFAULT 
                 WRITE(*,10)'In set_bc_neq_1D.F90, bc not implemented yet'
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

               CASE('sup_in')
                 get_ghost_state_Impl_1D_1st(i)%bc_procedure => sup_in_neq_1DImpl

               CASE ('sup_out','trans')
                 get_ghost_state_Impl_1D_1st(i)%bc_procedure => sup_out_neq_1DImpl 

               CASE ('sub_in_rhoiin_Tin','sub_in_rhoiin_Tin_Tvin', &
                   & 'sub_in_rhoiin_Tin_Trotin',                   &  
                   & 'sub_in_rhoiin_Tin_Trotin_Tvin',              & 
                   & 'sub_in_rhoiin_Tin_Trotin_Tvin_Tein',         & 
                   & 'sub_in_rhoiin_Tin_Tein',                     &
                   & 'sub_in_rhoiin_Tein')
                 get_ghost_state_Impl_1D_1st(i)%bc_procedure => sub_in_neq_rhoiin_T_1DImpl 
             
               CASE('sub_out_Tout')
                 get_ghost_state_Impl_1D_1st(i)%bc_procedure => sub_out_neq_Tout_1DImpl
 
               CASE DEFAULT 
                 WRITE(*,10)'In "set_bc_neq_1D.F90", bc not implemented yet...'
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

               CASE ('sup_in')
                 get_ghost_state_Impl_1D_2nd(i)%bc_procedure => sup_in_neq_1DImpl_2nd

               CASE ('sup_out','trans')
                 get_ghost_state_Impl_1D_2nd(i)%bc_procedure => sup_out_neq_1DImpl_2nd

               CASE ('sub_in_rhoiin_Tin','sub_in_rhoiin_Tin_Tvin', & 
                   & 'sub_in_rhoiin_Tin_Trotin',                   &
                   & 'sub_in_rhoiin_Tin_Trotin_Tvin',              & 
                   & 'sub_in_rhoiin_Tin_Trotin_Tvin_Tein',         & 
                   & 'sub_in_rhoiin_Tin_Tein',                     &
                   & 'sub_in_rhoiin_Tein')
                 get_ghost_state_Impl_1D_2nd(i)%bc_procedure => sub_in_neq_rhoiin_T_1DImpl_2nd 
              
               CASE DEFAULT 
                 WRITE(*,10)'In "set_bc_neq_1D.F90", bc not implemented yet...'
                 WRITE(*,10)bc(i)
                 PRINT*
                 STOP 

             END SELECT 

         ENDDO

       ENDIF

   ENDIF

10 FORMAT(A)

  END SUBROUTINE set_bc_neq_1D 
!------------------------------------------------------------------------------!
