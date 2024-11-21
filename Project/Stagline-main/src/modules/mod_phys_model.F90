!------------------------------------------------------------------------------!
!> This modules stores general data for the simulation under the user def
  MODULE mod_physical_model

    IMPLICIT NONE
  
    ! User defined physical model type
    TYPE model
      PRIVATE
      INTEGER :: nb_dim, nb_ns, nb_tvib, nb_te, nb_trot, nb_temp, nb_eq
      REAL(KIND=8) :: gamma, pr, rgas
      CHARACTER*80 :: label, mixture, transf, reaction, library
      LOGICAL :: diss
    END TYPE model
 
    TYPE(model), SAVE :: physical_model
 
    ! Functions for initializing and handling the physical model under use
    CONTAINS
  
      !------------------------------------------------------!
      !> Function initializing physical model data read in the input file
      FUNCTION set_phys_model (phys, diss, ndim, ns, ntvib, nte, ntrot, mixture, transf, reaction, library, pr, rgas, gamma) 
  
        INTEGER, OPTIONAL, INTENT(IN) :: ndim, ns, ntvib, nte, ntrot 
        REAL(KIND=8), OPTIONAL, INTENT(IN) :: rgas, gamma, pr
        CHARACTER*(*), OPTIONAL, INTENT(IN) :: phys, mixture, transf, reaction, library
        LOGICAL, OPTIONAL, INTENT(IN) :: diss
        INTEGER :: neq, ntemp 
        TYPE(model) :: set_phys_model
  
        ! Set default values
        set_phys_model = model (0, 0, 0, 0, 0, 1, 0, 1.4d0, 0.72d0, 287.06d0, 'none', 'none', 'none', 'none', 'none', .FALSE.) 
       
        ! Filling the physical model data
        IF (PRESENT(ndim)) set_phys_model%nb_dim = ndim
        IF (PRESENT(ns)) set_phys_model%nb_ns = ns
        IF (PRESENT(ntvib)) set_phys_model%nb_tvib = ntvib
        IF (PRESENT(nte)) set_phys_model%nb_te = nte
        IF (PRESENT(ntrot)) set_phys_model%nb_trot = ntrot
        
        ! Computing the number of temperatures
        ntemp = 1 + nte + ntvib + ntrot
        set_phys_model%nb_temp = ntemp  
  
        ! Computing the number of conservations equations
        neq = ntemp + ndim + ns
        set_phys_model%nb_eq = neq
  
        ! Filling the physical model data (calorically perfect gas only)
        IF (PRESENT(gamma)) set_phys_model%gamma = gamma
        IF (PRESENT(pr)) set_phys_model%pr = pr 
        IF (PRESENT(rgas)) set_phys_model%rgas = rgas
  
        ! Label of the present physical model
        IF (PRESENT(phys)) set_phys_model%label = phys
  
        ! Mixture, transfer mechanism and reaction files
        IF (PRESENT(mixture)) set_phys_model%mixture = mixture
        IF (PRESENT(transf)) set_phys_model%transf = transf
        IF (PRESENT(reaction)) set_phys_model%reaction = reaction
  
        ! Flag for dissipation effects (Navier-Stokes model)
        IF (PRESENT(diss)) set_phys_model%diss = diss
  
        ! Name of the thermodynamic library
        IF (PRESENT(library)) set_phys_model%library = library
  
      END FUNCTION set_phys_model
  
      !------------------------------------------------------!
      !> This function gives the specific heat ratio (it should be used only for calorically perfect gas simulation)
      FUNCTION get_gamma (phys_model)
   
         TYPE(model), INTENT(IN) :: phys_model   
         REAL(KIND=8) :: get_gamma
   
         get_gamma = phys_model%gamma
  
      END FUNCTION get_gamma
  
      !------------------------------------------------------!
      !> This function gives the gas constant r (ru/mm) (it should be used only for calorically perfect gas simulation)
      FUNCTION get_rgas (phys_model)
   
         TYPE(model), INTENT(IN) :: phys_model   
         REAL(KIND=8) :: get_rgas
   
         get_rgas = phys_model%rgas
  
      END FUNCTION get_rgas
  
      !------------------------------------------------------!
      !> This function gives the prandtl number (it should be used only for calorically perfect gas simulation)
      FUNCTION get_prandtl (phys_model)
   
         TYPE(model), INTENT(IN) :: phys_model   
         REAL(KIND=8) :: get_prandtl
   
         get_prandtl = phys_model%pr
  
      END FUNCTION get_prandtl 
  
      !------------------------------------------------------!
      !> This function gives the number of conservation equations
      FUNCTION get_nb_eq (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        INTEGER :: get_nb_eq
  
        get_nb_eq = phys_model%nb_eq
  
      END FUNCTION get_nb_eq
  
      !------------------------------------------------------!
      !> This function gives the number of dimensions
      FUNCTION get_nb_dim (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        INTEGER :: get_nb_dim
  
        get_nb_dim = phys_model%nb_dim
  
      END FUNCTION get_nb_dim
  
      !------------------------------------------------------!
      !> This function gives the number of temperatures
      FUNCTION get_nb_temp (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        INTEGER :: get_nb_temp
  
        get_nb_temp = phys_model%nb_temp
  
      END FUNCTION get_nb_temp
  
      !------------------------------------------------------!
      !> This function gives the number of species
      FUNCTION get_nb_species (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        INTEGER :: get_nb_species
  
        get_nb_species = phys_model%nb_ns
  
      END FUNCTION get_nb_species
  
      !------------------------------------------------------!
      !> This function gives the number of vibrational temperatures
      FUNCTION get_nb_tvib (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        INTEGER :: get_nb_tvib
  
        get_nb_tvib = phys_model%nb_tvib
  
      END FUNCTION get_nb_tvib
  
      !------------------------------------------------------!
      !> This function gives the number of rotational temperatures
      FUNCTION get_nb_trot (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        INTEGER :: get_nb_trot
  
        get_nb_trot = phys_model%nb_trot
  
      END FUNCTION get_nb_trot
  
      !------------------------------------------------------!
      !> This function gives the number of free-electron electronic temperatures
      FUNCTION get_nb_te (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        INTEGER :: get_nb_te
  
        get_nb_te = phys_model%nb_te
  
      END FUNCTION get_nb_te
  
      !------------------------------------------------------!
      !> This function tells if the physical model makes of a free-electron electronic temperature
      !! (-0 no Te in use -1 Te in use)
      FUNCTION presence_te (phys_model)
  
        TYPE(model), INTENT(IN) :: phys_model   
        INTEGER :: presence_te 
  
        presence_te = 0
        IF (phys_model%nb_te.NE.0) presence_te = 1
  
      END FUNCTION presence_te 
    
      !------------------------------------------------------!
      !> This function tells if the physical model makes of a rotational temperatures
      !! (-0 no Tvib in use -1 Tvib in use)
      FUNCTION presence_trot (phys_model)
  
        TYPE(model), INTENT(IN) :: phys_model   
        INTEGER :: presence_trot 
  
        presence_trot = 0
        IF (phys_model%nb_trot.NE.0) presence_trot = 1
  
      END FUNCTION presence_trot 
  
      !------------------------------------------------------!
      !> This function tells if the physical model makes of a vibrational temperatures
      !! (-0 no Tvib in use -1 Tvib in use)
      FUNCTION presence_tvib (phys_model)
  
        TYPE(model), INTENT(IN) :: phys_model   
        INTEGER :: presence_tvib 
  
        presence_tvib = 0
        IF (phys_model%nb_tvib.NE.0) presence_tvib = 1
  
      END FUNCTION presence_tvib 
  
      !------------------------------------------------------!
      !> This function gives the name thermodynamic library name
      FUNCTION get_library (phys_model)
  
        TYPE(model), INTENT(IN) :: phys_model   
        CHARACTER*80 :: get_library 
  
        get_library = phys_model%library
  
      END FUNCTION get_library
  
      !------------------------------------------------------!
      !> This function gives the mixture name
      FUNCTION get_mixture (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        CHARACTER*80 :: get_mixture
  
        get_mixture = phys_model%mixture
  
      END FUNCTION get_mixture
  
      !------------------------------------------------------!
      !> This function gives the reaction mechanism name
      FUNCTION get_reaction (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        CHARACTER*80 :: get_reaction
  
        get_reaction = phys_model%reaction
  
      END FUNCTION get_reaction
  
      !------------------------------------------------------!
      !> This function gives the transfer mechanism name
      FUNCTION get_transfer (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        CHARACTER*80 :: get_transfer
  
        get_transfer =  phys_model%transf
  
      END FUNCTION get_transfer 
  
      !------------------------------------------------------!
      !> This function gives the physical name
      FUNCTION get_name (phys_model)    
  
        TYPE(model), INTENT(IN) :: phys_model   
        CHARACTER*80 :: get_name
  
        get_name = phys_model%label
  
      END FUNCTION get_name
  
      !-------------------------------------------------------!
      !> This function gives a logical flag related to the inclusion of dissipation phenomena in the flowfield 
      FUNCTION get_diss_flag (phys_model)
  
        TYPE(model), INTENT(IN) :: phys_model 
        LOGICAL :: get_diss_flag
   
        get_diss_flag = phys_model%diss
  
      END FUNCTION
  
  END MODULE mod_physical_model
!------------------------------------------------------------------------------!

