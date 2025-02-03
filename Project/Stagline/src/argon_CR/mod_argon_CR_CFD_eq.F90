!------------------------------------------------------------------------------!
! This modules provides subroutines for the computation of equilibrium properties of
! the em, Ar, Arp system.
  MODULE mod_argon_CR_CFD_eq 

    IMPLICIT NONE

    ! Subroutines for equilibrium properties
    CONTAINS 

      !----------------------------------------------------!
      ! This subroutine computes the equilibrium composition (in terms of molar)
      ! fraction given pressure and temperature
      SUBROUTINE eq_composition (p, T, x_em, x_Ar, x_Arp)

        USE mod_argon_CR_initialize_CFD,      ONLY: ge, ukb, una, mm_em, mm_Ar, mm_Arp, & 
                                                  & R_Ar, Eion, gi_Ar, EiJ_Ar 

        INTEGER :: i
        REAL(KIND=8), PARAMETER :: xi_tol = 1.d-25
        REAL(KIND=8) :: mm, phi, rhs, rho, no, y_em, y_Ar, y_Arp
        REAL(KIND=8) :: Qt_em, Qt_Ar, Qt_Arp, Qel_Ar, Qel_Arp

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(OUT) :: x_em, x_Ar, x_Arp

        ! Translational partition functions of em, Ar and Arp
        ! em
        CALL tr_part_function (T, mm_em, Qt_em)

        ! Electron spin contribution
        Qt_em = ge*Qt_em

        ! Ar and Arp
        CALL tr_part_function (T, mm_Ar, Qt_Ar)
        CALL tr_part_function (T, mm_Arp, Qt_Arp)

        ! Internal partition functions of Ar and Arp species
        CALL int_part_function (T, 'Ar', Qel_Ar)
        CALL int_part_function (T, 'Arp', Qel_Arp) 
        
        ! Computation of ionization degree (Saha equation)
        rhs = (ukb*T/p)*Qt_em*(Qt_Arp*Qel_Arp)/(Qt_Ar*Qel_Ar)*DEXP(-Eion/(ukb*T))
        phi = DSQRT(rhs/(rhs + 1.d0))
        
        ! Mixture density
        rho = p/((1.d0 + phi)*T*R_Ar)

        no = rho*una/mm_Ar

        ! Species mass fractions
        y_em  = phi*no*mm_em/(rho*una)
        y_Ar  = (1.d0 - phi)*no*mm_Ar/(rho*una)
        y_Arp = phi*no*mm_Arp/(rho*una)

        ! Conversion from mass fractions to mole fractions
        mm = y_em/mm_em + y_Ar/mm_Ar + y_Arp/mm_Arp
        mm = 1.d0/mm

        ! A fix is applied in order to avoid numerical problems
        x_em  = MAX(y_em*mm/mm_em,xi_tol)
        x_Ar  = MAX(y_Ar*mm/mm_Ar,xi_tol)
        x_Arp = MAX(y_Arp*mm/mm_Arp,xi_tol)
        
      END SUBROUTINE eq_composition

      !----------------------------------------------------!
      ! This subroutine computes the equilibrium composition in terms species densities
      SUBROUTINE eq_composition_CR (p, T, rhoi)

        USE mod_argon_CR_initialize_CFD,      ONLY: pos_em, pos_Ar, pos_Arp, levels_Ar, levels_Arp,    & 
                                                  & ukb, urg, R_em, R_Ar, R_Arp, mm_em, mm_Ar, mm_Arp, & 
                                                  & gi_Ar, gi_Arp, EiJ_Ar, EiJ_Arp, model 

        INTEGER :: i
        REAL(KIND=8) :: ov_kbT
        REAL(KIND=8) :: x_em, x_Ar, x_Arp
        REAL(KIND=8) :: rho_em, rho_Ar, rho_Arp
        REAL(KIND=8) :: Qel_Ar, Qel_Arp

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi

        ! Computation of macroscopic chemical composition (molar fractions)
        CALL eq_composition (p, T, x_em, x_Ar, x_Arp)

        ! Density of em, Ar and Arp
        rho_em  = p*x_em/(R_em*T)
        rho_Ar  = p*x_Ar/(R_Ar*T)
        rho_Arp = p*x_Arp/(R_Arp*T)
  
        ! Internal partition functions of Ar and Ap species
        CALL int_part_function (T, 'Ar', Qel_Ar)
        CALL int_part_function (T, 'Arp', Qel_Arp)

        ! Useful quantities
        ov_kbT  = 1.d0/(ukb*T)
        Qel_Ar  = rho_Ar/Qel_Ar
        Qel_Arp = rho_Arp/Qel_Arp

        ! em
        rhoi(pos_em) = rho_em

        SELECT CASE(model)

          ! CR model is being used
          CASE('CR')

            ! Ar - electronic levels
            DO i = 1,levels_Ar
              rhoi(pos_Ar + i - 1) = gi_Ar(i)*DEXP(-EiJ_Ar(i)*ov_kbT)*Qel_Ar
            ENDDO

            ! Arp - electronic levels
            DO i = 1,levels_Arp
               rhoi(pos_Arp + i - 1) = gi_Arp(i)*DEXP(-EiJ_Arp(i)*ov_kbT)*Qel_Arp 
            ENDDO
        
          ! MT is being used
          CASE('MT')
       
            rhoi(pos_Ar)  = rho_Ar
            rhoi(pos_Arp) = rho_Arp 

        END SELECT
        
      END SUBROUTINE eq_composition_CR

      !----------------------------------------------------!
      ! This subroutine computes the translational patition function
      SUBROUTINE tr_part_function (T, mm, Qtr)

        USE mod_argon_CR_initialize_CFD,     ONLY: fac_Q

        REAL(KIND=8), INTENT(IN)  :: T, mm
        REAL(KIND=8), INTENT(OUT) :: Qtr

        Qtr = fac_Q*(mm*T)**1.5

      END SUBROUTINE tr_part_function 

      !----------------------------------------------------!
      ! This subroutine computes the internal partition function 
      SUBROUTINE int_part_function (T, species, Q)

        USE mod_argon_CR_initialize_CFD,     ONLY: ukb, levels_Ar, levels_Arp, gi_Ar, gi_Arp, & 
                                                 & EiJ_Ar, EiJ_Arp 

        INTEGER :: i
        REAL(KIND=8) :: ov_kbT

        REAL(KIND=8), INTENT(IN) :: T
        CHARACTER*(*), INTENT(IN) :: species
        REAL(KIND=8), INTENT(OUT) :: Q

        ! Initialization
        Q = 0.d0

        ! Useful quantity
        ov_kbT = 1.d0/(ukb*T)

        SELECT CASE(species)
 
          ! Ar 
          CASE('Ar')
            DO i = 1,levels_Ar
               Q = Q + gi_Ar(i)*DEXP(-EiJ_Ar(i)*ov_kbT)
            ENDDO

          ! Arp
          CASE('Arp')
            DO i = 1,levels_Arp
               Q = Q + gi_Arp(i)*DEXP(-EiJ_Arp(i)*ov_kbT)
            ENDDO

          CASE DEFAULT
            PRINT*
            WRITE(*,10)'in "mod_argon_CR_CFD_eq"'
            WRITE(*,10)'error in species selection for internal partition function'
            STOP

        END SELECT 

10    FORMAT(A)

      END SUBROUTINE int_part_function

  END MODULE mod_argon_CR_CFD_eq
!------------------------------------------------------------------------------!
