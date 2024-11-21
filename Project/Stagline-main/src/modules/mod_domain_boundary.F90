!------------------------------------------------------------------------------!
!> This module stores data for the application of boundary conditions.
  MODULE mod_domain_boundary

    IMPLICIT NONE
 
    ! User defined type boundary
    TYPE boundary 
      PRIVATE
      INTEGER :: id 
      REAL(KIND=8) :: p0, T0, pin, Tin, uin, vin, rhoin, pout, alpha, Tout, Tein, Twall, dv_dyin 
    END TYPE boundary
 
    INTEGER, SAVE :: nb_bc
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: yi_bound, rhoi_bound, Tr_bound, Tv_bound, & 
                                                   & inlet, inlet_temp
    CHARACTER*80, ALLOCATABLE, DIMENSION(:), SAVE :: bc, name_bc
    TYPE(boundary), ALLOCATABLE, DIMENSION(:), SAVE :: boundary_data
 
    ! Functions for handling boundary data
    CONTAINS 
  
      !------------------------------------------------------!
      !> This subroutine sets the specified boundary identifier.
      SUBROUTINE set_boundary_id (id, domain_boundary)
  
        INTEGER, INTENT(IN) :: id
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%id = id 
      
      END SUBROUTINE set_boundary_id 
  
      !------------------------------------------------------!
      !> This subroutine sets the wall temperature at the specified boundary.
      SUBROUTINE set_boundary_alpha (alpha, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: alpha
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%alpha = alpha 
    
      END SUBROUTINE set_boundary_alpha 
  
      !------------------------------------------------------!
      !> This subroutine sets the wall temperature at the specified boundary.
      SUBROUTINE set_boundary_Twall (T, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: T
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%Twall = T 
    
      END SUBROUTINE set_boundary_Twall 
  
      !------------------------------------------------------!
      !> This subroutine sets the inlet total pressure at the specified boundary.
      SUBROUTINE set_boundary_p0 (p0, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: p0
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%p0 = p0 
    
      END SUBROUTINE set_boundary_p0 
  
      !------------------------------------------------------!
      !> This subroutine sets the inlet total temperature at the specified boundary.
      SUBROUTINE set_boundary_T0 (T0, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: T0
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%T0 = T0 
    
      END SUBROUTINE set_boundary_T0
  
      !------------------------------------------------------!
      !> This subroutine sets the inlet static density at the specified boundary.
      SUBROUTINE set_boundary_rhoin (rho, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: rho
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%rhoin = rho 
    
      END SUBROUTINE set_boundary_rhoin 
  
      !------------------------------------------------------!
      !> This subroutine sets the x-component of flow velocity at the specified boundary.
      SUBROUTINE set_boundary_uin (u, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: u
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%uin = u 
  
      END SUBROUTINE set_boundary_uin
  
      !------------------------------------------------------!
      !> This subroutine sets the y-component of flow velocity at the specified boundary.
      SUBROUTINE set_boundary_vin (v, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: v
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%vin = v 
  
      END SUBROUTINE set_boundary_vin
  
      !------------------------------------------------------!
      !> This subroutine sets the tangential-velocity derivative at the specified boundary.
      SUBROUTINE set_boundary_dv_dyin(dv_dy, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: dv_dy
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%dv_dyin = dv_dy 
  
      END SUBROUTINE set_boundary_dv_dyin
  

      !------------------------------------------------------!
      !> This subroutine sets the inlet static pressure at the specified boundary.
      SUBROUTINE set_boundary_pin (p, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: p
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%pin = p 
  
      END SUBROUTINE set_boundary_pin
  
      !------------------------------------------------------!
      !> This subroutine sets the outlet static pressure at the specified boundary.
      SUBROUTINE set_boundary_pout (p, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: p
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%pout = p 
  
      END SUBROUTINE set_boundary_pout
  
      !------------------------------------------------------!
      !> This subroutine sets the inlet static temperature at the specified boundary.
      SUBROUTINE set_boundary_Tin (T, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: T
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%Tin = T 
  
      END SUBROUTINE set_boundary_Tin
  
      !------------------------------------------------------!
      !> This subroutine sets the inlet free-electron electronic temperature at the specified boundary.
      SUBROUTINE set_boundary_Tein (T, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: T
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%Tein = T 
  
      END SUBROUTINE set_boundary_Tein
  
      !------------------------------------------------------!
      !> This subroutine sets the outlet static temperature at the specified boundary.
      SUBROUTINE set_boundary_Tout (T, domain_boundary)
  
        REAL(KIND=8), INTENT(IN) :: T
        TYPE(boundary), INTENT(INOUT) :: domain_boundary
    
        domain_boundary%Tout = T 
  
      END SUBROUTINE set_boundary_Tout
  
      !------------------------------------------------------!
      !> This function gives the variable imposed at a given boundary for a supersonic inlet 
      !! (the number of equation must be provided by the user). 
      FUNCTION get_boundary_inlet (neq, side) 
  
        INTEGER, INTENT(IN) :: neq
        INTEGER :: i, id
        REAL(KIND=8) , DIMENSION(neq) :: get_boundary_inlet
        TYPE(boundary), INTENT(IN) :: side
  
        ! Boundary id
        id = side%id
  
        ! Data for the supersonic inlet boundary condition
        DO i = 1,neq
           get_boundary_inlet(i) = inlet(neq*(id - 1) + i)
        ENDDO
  
      END FUNCTION get_boundary_inlet
  
      !------------------------------------------------------!
      !> This function provides species mass fractions at a given boundary 
      !! (the number of species must be provided in input by the user).
      FUNCTION get_boundary_yi (ns, side)
  
        INTEGER, INTENT(IN) :: ns
        INTEGER :: i, id
        REAL(KIND=8), DIMENSION(ns) :: get_boundary_yi
        TYPE(boundary), INTENT(IN) :: side
  
        ! Boundary id
        id = side%id 
    
        DO i = 1,ns 
           get_boundary_yi(i) = yi_bound(ns*(id - 1) + i)
        ENDDO
  
      END FUNCTION 
  
      !------------------------------------------------------!
      !> This function provides species densities at a given boundary 
      !! (the number of species must be provided in input by the user).
      FUNCTION get_boundary_rhoi (ns, side)
  
        INTEGER, INTENT(IN) :: ns
        INTEGER :: i, id
        REAL(KIND=8), DIMENSION(ns) :: get_boundary_rhoi
        TYPE(boundary), INTENT(IN) :: side
  
        ! Boundary id
        id = side%id
    
        DO i = 1,ns 
           get_boundary_rhoi(i) = rhoi_bound(ns*(id - 1) + i)
        ENDDO
  
      END FUNCTION get_boundary_rhoi
 
      !------------------------------------------------------!
      !> This function provides rotational temperatures at a given boundary 
      !! (the number of rotational temperatures must be provided in input by the user).
      FUNCTION get_boundary_Tr (nr, side)
  
        INTEGER, INTENT(IN) :: nr
        INTEGER :: i, id, pos
        REAL(KIND=8), DIMENSION(nr) :: get_boundary_Tr
        TYPE(boundary), INTENT(IN) :: side
    
        ! Boundary id
        id = side%id 
        
        pos = nr*(id - 1)
    
        DO i = 1,nr
           get_boundary_Tr(i) = Tr_bound(pos + i)
        ENDDO
  
      END FUNCTION get_boundary_Tr
 
      !------------------------------------------------------!
      !> This function provides vibrational temperatures at a given boundary
      !! (the number of vibrational temperatures must be provided in input by the user).
      FUNCTION get_boundary_Tv (nv, side)
  
        INTEGER, INTENT(IN) :: nv
        INTEGER :: i, id, pos
        REAL(KIND=8), DIMENSION(nv) :: get_boundary_Tv
        TYPE(boundary), INTENT(IN) :: side
    
        ! Boundary id
        id = side%id 
    
        pos = nv*(id - 1)
        DO i = 1,nv 
           get_boundary_Tv(i) = Tv_bound(pos + i)
        ENDDO
  
      END FUNCTION get_boundary_Tv 
  
      !------------------------------------------------------!
      !> This function provides the inlet static density at a given boundary.
      FUNCTION get_boundary_rhoin (side)
  
        REAL(KIND=8) :: get_boundary_rhoin
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_rhoin = side%rhoin
  
      END FUNCTION get_boundary_rhoin
  
 
      !------------------------------------------------------!
      !> This function provides the derivative of the velocity component normal
      !to the stagnation line (v physical) in the direction normal to the stagnation line 
      FUNCTION get_boundary_dv_dyin (side)
  
        REAL(KIND=8) :: get_boundary_dv_dyin
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_dv_dyin = side%dv_dyin
  
      END FUNCTION get_boundary_dv_dyin
  
      !------------------------------------------------------!
      !> This function provides the inlet static temperature at a given boundary.
      FUNCTION get_boundary_Tin (side)
  
        REAL(KIND=8) :: get_boundary_Tin
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_Tin = side%Tin
  
      END FUNCTION get_boundary_Tin 
  
      !------------------------------------------------------!
      !> This function provides the inlet free-electron electronic temperature at a given boundary.
      FUNCTION get_boundary_Tein (side)
  
        REAL(KIND=8) :: get_boundary_Tein
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_Tein = side%Tein
  
      END FUNCTION get_boundary_Tein   
 
      !------------------------------------------------------!
      !> This function provides the intel temperature vector 
      !! (the number of temperatures must be provides in input by the user)
      FUNCTION get_boundary_Tvec (nt, side)

        INTEGER :: i, pos, id
        REAL(KIND=8), DIMENSION(nt):: get_boundary_Tvec
        INTEGER, INTENT(IN) :: nt
        TYPE(boundary), INTENT(IN) :: side

        ! Boundary id
        id = side%id 

        pos = (id - 1)*nt 
        DO i = 1,nt
           get_boundary_Tvec(i) = inlet_temp(pos + i) 
        ENDDO

      END FUNCTION get_boundary_Tvec
 
      !------------------------------------------------------!
      !> This function provides the inlet static pressure at a given boundary.
      FUNCTION get_boundary_pin (side)
  
        REAL(KIND=8) :: get_boundary_pin
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_pin = side%pin
  
      END FUNCTION get_boundary_pin 
  
      !------------------------------------------------------!
      !> This function provides the inlet total pressure at a given boundary.
      FUNCTION get_boundary_p0 (side)
  
        REAL(KIND=8) :: get_boundary_p0
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_p0 = side%p0
  
      END FUNCTION get_boundary_p0 
  
      !------------------------------------------------------!
      !> This function provides the inlet total temperature at a given boundary.
      FUNCTION get_boundary_T0 (side)
  
        REAL(KIND=8) :: get_boundary_T0
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_T0 = side%T0
  
      END FUNCTION get_boundary_T0
  
      !------------------------------------------------------!
      !> This function provides the x-component of inlet velocity at a given boundary.
      FUNCTION get_boundary_uin (side)
  
        REAL(KIND=8) :: get_boundary_uin
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_uin = side%uin
  
      END FUNCTION get_boundary_uin
  
      !------------------------------------------------------!
      !> This function provides the y-component of inlet velocity at a given boundary.
      FUNCTION get_boundary_vin (side)
  
        REAL(KIND=8) :: get_boundary_vin
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_vin = side%vin
  
      END FUNCTION get_boundary_vin
  
      !------------------------------------------------------!
      !> This function provides the outlet pressure at a given boundary.
      FUNCTION get_boundary_pout (side)
  
        REAL(KIND=8) :: get_boundary_pout
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_pout = side%pout
  
      END FUNCTION get_boundary_pout
  
      !------------------------------------------------------!
      !> This function provides the outlet temperature at a given boundary.
      FUNCTION get_boundary_Tout (side)
  
        REAL(KIND=8) :: get_boundary_Tout
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_Tout = side%Tout
  
      END FUNCTION get_boundary_Tout
  
      !------------------------------------------------------!
      !> This function provides the wall temperature at a given boundary.
      FUNCTION get_boundary_Twall (side)
  
        REAL(KIND=8) :: get_boundary_Twall
        TYPE(boundary), INTENT(IN) :: side
    
        get_boundary_Twall = side%Twall
  
      END FUNCTION get_boundary_Twall
  
  END MODULE mod_domain_boundary
!------------------------------------------------------------------------------!
