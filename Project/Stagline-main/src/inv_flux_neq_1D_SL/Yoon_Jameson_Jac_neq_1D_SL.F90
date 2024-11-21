!------------------------------------------------------------------------------!
!> This subrotine provides the numerical flux Jacobian for 1D stagnation line nonequilibrium flows
!! according to the approximation proposed by Yoon and Jameson.
  SUBROUTINE Yoon_Jameson_Jac_neq_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, jfl, jfr)

    USE mod_general_data,             ONLY: nb_eq, pos_u_cell, pos_c_cell
    USE mod_function_pointer,         ONLY: inv_flux_Jac_neq_1D_SL

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8) :: cl, cr, Mnl, Mnr, spl, spr, ul, ur, unl, unr 

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfl      !< numerical flux Jacobian with respect to the left state 
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfr      !< numerical flux Jacobian with respect to the right state 

    ! Speed of sound and flow speed of left and right states 
    cl = left_data(pos_c_cell)
    ul = left_data(pos_u_cell)
  
    cr = right_data(pos_c_cell)
    ur = right_data(pos_u_cell)

    ! Normal velocities of left and right states
    unl = ul*nx
    unr = ur*nx
  
    ! Normal Mach numbers of left and right states
    Mnl = unl/cl
    Mnr = unr/cr
   
    ! Spectral radii of left and right state inviscid flux Jacobians
    spl = ABS(unl) + cl
    spr = ABS(unr) + cr 

    ! Left state 
    IF (Mnl.GE.1.d0) THEN

       ! Supersonic flow
       CALL inv_flux_Jac_neq_1D_SL(nx, u_left, left_data, jfl) 

    ELSE 

       ! Subsonic flow
       IF (Mnl.GT.-1.d0) THEN

          CALL inv_flux_Jac_neq_1D_SL(nx, u_left, left_data, jfl) 

          DO j = 1,nb_eq 

             DO i = 1,j - 1
                jfl(i,j) = 0.5d0*jfl(i,j)
             ENDDO

             jfl(j,j) = 0.5d0*(jfl(j,j) + spl)

             DO i = j + 1,nb_eq 
                jfl(i,j) = 0.5d0*jfl(i,j)
             ENDDO

          ENDDO

       ELSE 

          DO j = 1,nb_eq
             DO i = 1,nb_eq 
                jfl(i,j) = 0.d0
             ENDDO
          ENDDO 

       ENDIF

    ENDIF

    ! Right state 
    IF (Mnr.LE.-1.d0) THEN

       ! Supersonic flow
       CALL inv_flux_Jac_neq_1D_SL(nx, u_right, right_data, jfr) 

    ELSE 

       ! Subsonic flow
       IF (Mnr.LT.1.d0) THEN

          CALL inv_flux_Jac_neq_1D_SL(nx, u_right, right_data, jfr)

          DO j = 1,nb_eq 

             DO i = 1,j - 1
               jfr(i,j) = 0.5d0*jfr(i,j)
             ENDDO

             jfr(j,j) = 0.5d0*(jfr(j,j) - spr)

             DO i = j + 1,nb_eq 
                jfr(i,j) = 0.5d0*jfr(i,j)
             ENDDO

          ENDDO

       ELSE 

          DO j = 1,nb_eq 
             DO i = 1,nb_eq 
                jfr(i,j) = 0.d0
             ENDDO
          ENDDO 

       ENDIF
     
    ENDIF

  END SUBROUTINE Yoon_Jameson_Jac_neq_1D_SL
!------------------------------------------------------------------------------!

