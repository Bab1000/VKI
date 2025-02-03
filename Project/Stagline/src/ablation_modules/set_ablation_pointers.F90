 SUBROUTINE set_ablation_pointers ()

#include "../config.h"

#ifdef CARBONABLA

 USE mod_general_data,        ONLY: ablation_library_name
 USE mod_ablation_pointers
 USE mod_ablation_procedure_carbon_library
 USE mod_ablation_procedures
                                    
 IMPLICIT NONE

   ! Ablation inizialization
   procedure_ablation_initialize => ablation_initialize   

   ! Ablation vector deallocation and pointers nullification 
   procedure_ablation_finalize => ablation_finalize   

   ! Ablation write solution file 
   procedure_ablation_write_sol => ablation_write_sol  

   ! Subroutine for computing ablation source terms of each species
   procedure_ablation_compute_wall_source_terms => ablation_compute_source_terms

   ! Subroutine for computing the surface energy balance
   procedure_ablation_compute_wall_SEB => ablation_compute_SEB

   ! Subroutine for computing the surface transport fluxes at the end of the
   ! simualtion
   procedure_ablation_comp_surface_flux => ablation_comp_surface_flux

   ! Subroutine for computing the wall conductive heat flux  and the stress tensor 
   procedure_ablation_wall_diff_flux => ablation_wall_stress_tensor_and_conductive_heat_flux

   ! Subroutine for computing the surface mass balance
   procedure_ablation_SMB => ablation_SMB

   ! Subroutine for computing ablation mass flow rates from the considered heterogeneous reactions (library chosen in input file)
   SELECT CASE (ablation_library_name)

           CASE('Park_original')
                   library_ablation_compute_wall_mass_blowing_rate => carbon_ablation_park_original_mdot_wall

           CASE('Park_Suzuki_nitri')
                   library_ablation_compute_wall_mass_blowing_rate => carbon_ablation_park_suzuki_nitri_mdot_wall 

           CASE('Park_no_nitri')
                   library_ablation_compute_wall_mass_blowing_rate => carbon_ablation_park_no_nitri_mdot_wall 

           CASE('Custom')
                   library_ablation_compute_wall_mass_blowing_rate => carbon_ablation_user_defined_mdot_wall 

           CASE('Custom_nitri')
                   library_ablation_compute_wall_mass_blowing_rate => carbon_nitridation_user_defined_mdot_wall

           CASE('Nitri_probability_rebuild')
                   library_ablation_compute_wall_mass_blowing_rate => carbon_nitridation_mdot_wall_imposed

           CASE('Park_oxi_rebuilt_nitri_and_cata')
                   library_ablation_compute_wall_mass_blowing_rate => carbon_oxidation_park_rebuilt_nitri_mdot_wall
   END SELECT

#endif

 END SUBROUTINE set_ablation_pointers

