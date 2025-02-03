!------------------------------------------------------------------------------!
! This module provides subroutine and functions for computing thermodynamic and transport 
! properties for the Polyat_gas_WCU library
  MODULE mod_polyat_gas_wcu_CFD_prop 

    IMPLICIT NONE
    
    CONTAINS 

      !----------------------------------------------------!
      ! This subroutine computes the specific energy of the gas
      SUBROUTINE get_energy (T, e)

        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: e_table

        INTEGER :: left, right
        REAL(KIND=8) :: a, b, dT_norm

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT) :: e

        ! Look-up table search
        CALL table_search(T, left, right, dT_norm)

        ! Specific energy
        a = e_table(left)
        b = e_table(right)
        e = a + (b - a)*dT_norm 

      END SUBROUTINE get_energy

      !----------------------------------------------------!
      ! This subroutine computes the constant volume specific of the gas
      SUBROUTINE get_cv (T, cv)

        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: cv_table

        INTEGER :: left, right
        REAL(KIND=8) :: a, b, dT_norm

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT) :: cv

        ! Look-up table search
        CALL table_search(T, left, right, dT_norm)

        ! Constant volume heat
        a  = cv_table(left)
        b  = cv_table(right)
        cv = a + (b - a)*dT_norm 

      END SUBROUTINE get_cv

      !----------------------------------------------------!
      ! This subroutine computes transport coefficients 
      SUBROUTINE get_transport_coeff (T, mu, kappa, lambda)
 
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: mu_table, kappa_table, lambda_table 

        INTEGER :: left, right
        REAL(KIND=8) :: a, b, dT_norm
 
        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT) :: mu, kappa, lambda

        ! Look-up table search
        CALL table_search(T, left, right, dT_norm)

        ! Shear viscosity
        a  = mu_table(left)
        b  = mu_table(right)
        mu = a + (b - a)*dT_norm 
   
        ! Bulk viscosity
        a  = kappa_table(left)
        b  = kappa_table(right)
        kappa = a + (b - a)*dT_norm 

        ! Thermal conductivity
        a  = lambda_table(left)
        b  = lambda_table(right)
        lambda = a + (b - a)*dT_norm

      END SUBROUTINE get_transport_coeff

      !----------------------------------------------------!
      ! This subroutine computes the temperature from the specific energy
      SUBROUTINE get_temperature(eint, T)

        REAL(KIND=8), PARAMETER :: tol = 1.d-8
        REAL(KIND=8) :: cv, e, res

        REAL(KIND=8), INTENT(IN) :: eint
        REAL(KIND=8), INTENT(OUT) :: T

        ! Guess value
        T = 1000.d0

        ! Start Newton's loop
        res = 1.d0
        DO WHILE (res.GT.tol)
 
           CALL get_cv(T, cv)
           CALL get_energy(T, e)

           ! Residual            
           res  = (e - eint)/cv

           ! Solution update
           T    = T - res
           res  = (ABS(res))
          
        ENDDO

      END SUBROUTINE get_temperature

      !----------------------------------------------------!
      ! This subroutine computes the post-shock equilibrium conditions 
      SUBROUTINE get_post_shock(p1, u1, T1, p2, u2, T2)

        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: solver, Rgas

        INTEGER :: i, length
        REAL(KIND=8), PARAMETER :: tol = 1.d-8
        REAL(KIND=8) :: ratio
        REAL(KIND=8) :: cp, h, mdot, res_rho, res_T
        REAL(KIND=8) :: rho1, rho2, rho2_old, ek1, ek2, h1, h2
        REAL(KIND=8), DIMENSION(3) :: flux_in, flux_out, res

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2

        ! Compute mass, momentum and energy fluxes
        CALL get_energy(T1, h1)
        h1   = h1 + Rgas*T1
        ek1  = 0.5d0*u1**2 
        rho1 = p1/(Rgas*T1)
 
        mdot = rho1*u1
        flux_in(1) = mdot
        flux_in(2) = p1 + mdot*u1
        flux_in(3) = h1 + ek1

        ! Guess value of the density ratio
        ratio   = 0.1d0
        res_rho = 1.d0
        rho2    = rho1*ratio
        DO WHILE (res_rho.GT.tol)

           u2 = u1*ratio
           p2 = p1 + mdot*u1*(1.d0 - ratio)
           h2 = h1 + ek1*(1.d0 - ratio**2)

           ! Newton's loop for the temperature
           h  = h2
           T2 = 1000.d0
           res_T = 1.d0
           DO WHILE (res_T.GT.tol)

              CALL get_cv(T2, cp)
              cp = cp + Rgas
              CALL get_energy(T2, h2)
              h2 = h2 + Rgas*T2 

              ! Residual            
              res_T = (h2 - h)/cp

              ! Solution update
              T2    = T2 - res_T
              res_T = (ABS(res_T))

           ENDDO

           ! Density ratio update
           rho2_old = rho2
           rho2  = p2/(Rgas*T2)
           ratio = rho1/rho2

           ! Density residual
           res_rho = ABS(rho2_old - rho2)/rho2
           
        ENDDO

        CALL get_energy(T2, h2)
        h2   = h2 + Rgas*T2
        ek2  = 0.5d0*u2**2 
        rho2 = p2/(Rgas*T2)
 
        mdot = rho2*u2
        flux_out(1) = mdot
        flux_out(2) = p2 + mdot*u2
        flux_out(3) = h2 + ek2

        DO i = 1,3
           res(i) = ABS(flux_in(i) - flux_out(i))/flux_in(i)*100.d0
        ENDDO

        length = LEN_TRIM(solver)
        WRITE(*,5)solver(1:length),':: Polyat_gas_WCU -> post-shock conditions'
        PRINT*
        WRITE(*,10)'Residuals on mass, momentum and energy fluxes:'
        PRINT*
        WRITE(*,15)'Mass      ',res(1),' [%]'
        PRINT*
        WRITE(*,15)'Momentum  ',res(2),' [%]'
        PRINT*
        WRITE(*,15)'Energy    ',res(3),' [%]'
        PRINT* 

5     FORMAT(A,A)
10    FORMAT(A)
15    FORMAT(A,E14.6,A)

      END SUBROUTINE get_post_shock

      !----------------------------------------------------!
      ! This subroutine searches the left and right states in the temperature table
      SUBROUTINE table_search(T, left, right, dT_norm)

        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: table_points, Tmin, Tmax, dT, T_table

        REAL(KIND=8), INTENT(IN) :: T
        INTEGER, INTENT(OUT) :: left, right 
        REAL(KIND=8), INTENT(OUT) :: dT_norm

        ! Left and right states
        left  = INT((T - Tmin)/dT) + 1 
        right = left + 1
        
        ! Normalized temperature difference
        dT_norm = (T  - T_table(left))/dT

      END SUBROUTINE table_search

  END MODULE mod_polyat_gas_wcu_CFD_prop
!------------------------------------------------------------------------------!
