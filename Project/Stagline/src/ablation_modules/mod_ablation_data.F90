MODULE mod_ablation_data

#include "../config.h"

#ifdef CARBONABLA

        IMPLICIT NONE 

        CHARACTER(LEN=:), ALLOCATABLE :: abla_dir                                                 !<Directory containing useful data for ablation computations

        INTEGER, ALLOCATABLE, DIMENSION(:,:)          :: ion_species_matrix                       !< ion-species pair index matrix
        INTEGER                                       :: C_pos      = -1                          !< C position in the mixture 
        INTEGER                                       :: C3_pos     = -1                          !< C3 position in the mixture 
        INTEGER                                       :: CN_pos     = -1                          !< CN position in the mixture 
        INTEGER                                       :: CO_pos     = -1                          !< CO position in the mixture 
        INTEGER                                       :: N_pos      = -1                          !< N position in the mixture 
        INTEGER                                       :: N2_pos     = -1                          !< N position in the mixture 
        INTEGER                                       :: O_pos      = -1                          !< O position in the mixture 
        INTEGER                                       :: O2_pos     = -1                          !< O2 position in the mixture 
        INTEGER                                       :: e_minus_pos = -1                         !< e- position in the mixture 


        REAL(KIND=8), ALLOCATABLE, DIMENSION(:)       :: omega_wall                               !< wall source term vector
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:)       :: omega_wall_rec                           !< recombination source term vector 
        REAL(KIND=8), DIMENSION(4)                    :: mdot_wall                                !< char mass blowing rate vector
        REAL(KIND=8)                                  :: tau_rr_wall, tau_rt_wall                 !< stress tensor components 
        REAL(KIND=8)                                  :: cond_flux_wall                           !< wall conductive heat flux
        REAL(KIND=8)                                  :: mdot_wall_tot                            !< total char mass blowing rate 
        REAL(KIND=8)                                  :: uwall                                    !< blowing velocity 
        REAL(KIND=8)                                  :: Twall                                    !< wall temperature 
        REAL(KIND=8)                                  :: seb_stored                               !< stored Surface Energy Balance
        REAL(KIND=8)                                  :: gamma_rec = 0.d0                         !< N-N2 recombination probability (default value if not given in input)
        REAL(KIND=8)                                  :: O_prob                                   !< O+C_s=>CO  pre-exponential coefficient of the reaction probability (needed only if user defined model is used)         
        REAL(KIND=8)                                  :: O2_prob                                  !< O2+2C_s=>2CO reaction probability (needed only if user defined model is used)         
        REAL(KIND=8)                                  :: N_prob                                   !< N+C_s=>CN reaction probability (needed only if user defined model is used)         
        REAL(KIND=8)                                  :: C_prob                                   !< 3C=>C3 vaporization coefficient (needed only if user defined model is used)         
        REAL(KIND=8), SAVE                            :: wall_emiss = 0.9d0                       !< wall emissivity (default value if not given in input)
        REAL(KIND=8), SAVE                            :: phi_pyro = 0.0d0                         !< ratio between pyrolysis and ablation mass injection
        REAL(KIND=8), SAVE                            :: exp_mdot = 0.0d0                         !< experimental mass blowing rate (needed when the reaction probability rebuilding is performed)

#endif

END MODULE mod_ablation_data
