  MODULE mod_ablation_pointers

#include "../config.h"

#ifdef CARBONABLA 

    IMPLICIT NONE

    !------------------------------------------------------!
    ! Interfaces for subroutines and related pointers
    !------------------------------------------------------!
    ABSTRACT INTERFACE

      !----------------------!
      !Procedure subroutines !
      !----------------------!

      ! Interface for subroutine initializing terms for wall ablation 
      SUBROUTINE proc_ablation_initialize (ns)
        INTEGER, INTENT(IN) :: ns
      END SUBROUTINE proc_ablation_initialize

      ! Interface for deallocating vectors and nullifying  pointers for wall ablation
      SUBROUTINE proc_ablation_finalize()
      END SUBROUTINE proc_ablation_finalize

      ! Interface for subroutine writing ablation output file 
      SUBROUTINE proc_ablation_write_sol ()
      END SUBROUTINE proc_ablation_write_sol

      ! Interface for subroutine computing the surface mass balance
      SUBROUTINE proc_ablation_SMB (phys_cons1,p_wall,rhoi_wall)
        REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: phys_cons1      !< conservative variables of the first physical cell
        REAL(KIND=8), INTENT(IN)                   :: p_wall          !< wall pressure 
        REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  :: rhoi_wall       !< wall species density vector 
      END SUBROUTINE proc_ablation_SMB

      ! Interface for subroutine computing the source terms in case of
      ! ablation
      SUBROUTINE proc_ablation_compute_wall_source_terms (rhoi_wall,Twall)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi_wall
        REAL(KIND=8), INTENT(IN) ::  Twall
      END SUBROUTINE proc_ablation_compute_wall_source_terms

      ! Interface for subroutine computing the surface energy balance
      SUBROUTINE proc_ablation_compute_wall_SEB (rhoi_wall,Twall,T_first_phys_state,seb)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi_wall
        REAL(KIND=8), INTENT(IN)  :: Twall
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: T_first_phys_state
        REAL(KIND=8), INTENT(OUT) :: seb
      END SUBROUTINE proc_ablation_compute_wall_SEB

      ! Interface for subroutine computing the surface transport fluxes at the
      ! end of the simulation
      SUBROUTINE proc_ablation_comp_surface_flux (qdiff_wall,qcond_wall, qabla_wall)
        REAL(KIND=8), INTENT(OUT) :: qdiff_wall, qabla_wall
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: qcond_wall
      END SUBROUTINE proc_ablation_comp_surface_flux

      ! Interface for subroutine computing the wall diffusive fluxes in
      ! case of
      ! ablation 
      SUBROUTINE proc_ablation_wall_diffusive_flux (r_l, r_r, voll, volr,left_data, right_data, u_left, u_right)
        REAL(KIND=8), INTENT(IN) :: r_l, r_r
        REAL(KIND=8), INTENT(IN) :: voll, volr
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left, u_right
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data
      END SUBROUTINE proc_ablation_wall_diffusive_flux

      !--------------------!
      !Library subroutines !
      !--------------------!

      ! Interface for subroutine computing the mass blowing terms in
      ! case of ablation
      SUBROUTINE lib_ablation_compute_wall_mass_blowing_rate (rhoi_wall,Twall,mdot_wall_tot)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi_wall
        REAL(KIND=8), INTENT(IN) :: Twall
        REAL(KIND=8), INTENT(OUT) :: mdot_wall_tot
      END SUBROUTINE lib_ablation_compute_wall_mass_blowing_rate 

    END INTERFACE


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

      !----------------------!
      !Procedure subroutines !
      !----------------------!

    PROCEDURE(proc_ablation_initialize), POINTER, SAVE ::procedure_ablation_initialize
    PROCEDURE(proc_ablation_finalize), POINTER, SAVE ::procedure_ablation_finalize
    PROCEDURE(proc_ablation_write_sol), POINTER, SAVE::procedure_ablation_write_sol
    PROCEDURE(proc_ablation_compute_wall_source_terms), POINTER, SAVE ::procedure_ablation_compute_wall_source_terms
    PROCEDURE(proc_ablation_compute_wall_SEB), POINTER, SAVE ::procedure_ablation_compute_wall_SEB
    PROCEDURE(proc_ablation_comp_surface_flux), POINTER, SAVE ::procedure_ablation_comp_surface_flux
    PROCEDURE(proc_ablation_wall_diffusive_flux), POINTER, SAVE :: procedure_ablation_wall_diff_flux
    PROCEDURE(proc_ablation_SMB), POINTER, SAVE :: procedure_ablation_SMB



      !--------------------!
      !Library subroutines !
      !--------------------!

    PROCEDURE(lib_ablation_compute_wall_mass_blowing_rate), POINTER, SAVE ::library_ablation_compute_wall_mass_blowing_rate


#endif

  END MODULE mod_ablation_pointers
