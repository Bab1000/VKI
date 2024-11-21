  MODULE mod_ablation_procedure_carbon_library


#include "../config.h"

#ifdef CARBONABLA

    IMPLICIT NONE

    CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE carbon_ablation_park_original_mdot_wall(rhoi_wall,T_wall,mdot_tot)

        USE mod_general_data,          ONLY: mi, Ri
        USE mod_ablation_data,         ONLY: mdot_wall, mdot_wall_tot,           &
                                           & O_pos, C_pos, O2_pos, N_pos, C3_pos 

        REAL(KIND=8) :: p_eq, rho_eq

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_wall       !< species density vector at the wall
        REAL(KIND=8), INTENT(IN)                :: T_wall          !< wall temperature
        REAL(KIND=8), INTENT(OUT)               :: mdot_tot        !< total char mass blowing rate  at the wall

        !O + C_s => CO
        mdot_wall(1)=rhoi_wall(O_pos)*dsqrt(Ri(O_pos)*T_wall/ &
        &(2.d0*acos(-1.d0)))*0.63d0*dexp(-1160.d0/T_wall)*mi(C_pos)/mi(O_pos)

        !O2 + 2C_s => 2CO
        mdot_wall(2)=rhoi_wall(O2_pos)*dsqrt(2.d0*Ri(O2_pos)*T_wall/(acos(-1.d0)))*0.5d0*mi(C_pos)/mi(O2_pos)

        !N + C_s => CN
        mdot_wall(3)=rhoi_wall(N_pos)*dsqrt(Ri(N_pos)*T_wall/(2.d0*acos(-1.d0)))*0.3d0*mi(C_pos)/mi(N_pos) 

        !3C_s => C3 
        p_eq=5.19d14*dexp(-90845.d0/T_wall)
        rho_eq=p_eq/(Ri(C3_pos)*T_wall)
        mdot_wall(4)=(rho_eq-rhoi_wall(C3_pos))*dsqrt(Ri(C3_pos)*T_wall/(2.d0*acos(-1.d0)))

        mdot_tot = SUM(mdot_wall)
        mdot_wall_tot = mdot_tot 

        END SUBROUTINE carbon_ablation_park_original_mdot_wall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE carbon_ablation_park_suzuki_nitri_mdot_wall(rhoi_wall,T_wall,mdot_tot)

        USE mod_general_data,          ONLY: mi, Ri
        USE mod_ablation_data,         ONLY: mdot_wall, mdot_wall_tot,           &
                                           & O_pos, C_pos, O2_pos, N_pos, C3_pos 

        REAL(KIND=8) :: p_eq, rho_eq

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_wall       !< species density vector at the wall
        REAL(KIND=8), INTENT(IN)                :: T_wall          !< wall temperature
        REAL(KIND=8), INTENT(OUT)               :: mdot_tot        !< total char mass blowing rate  at the wall

        !O + C_s => CO
        mdot_wall(1)=rhoi_wall(O_pos)*dsqrt(Ri(O_pos)*T_wall/ &
        &(2.d0*acos(-1.d0)))*0.63d0*dexp(-1160.d0/T_wall)*mi(C_pos)/mi(O_pos)

        !O2 + 2C_s => 2CO
        mdot_wall(2)=rhoi_wall(O2_pos)*dsqrt(2.d0*Ri(O2_pos)*T_wall/(acos(-1.d0)))*0.5d0*mi(C_pos)/mi(O2_pos)

        !N + C_s => CN
        mdot_wall(3)=rhoi_wall(N_pos)*dsqrt(Ri(N_pos)*T_wall/(2.d0*acos(-1.d0)))*0.3d-2*mi(C_pos)/mi(N_pos)      

        !3C_s => C3 
        p_eq=5.19d14*dexp(-90845.d0/T_wall)
        rho_eq=p_eq/(Ri(C3_pos)*T_wall)
        mdot_wall(4)=(rho_eq-rhoi_wall(C3_pos))*dsqrt(Ri(C3_pos)*T_wall/(2.d0*acos(-1.d0)))

        mdot_tot = SUM(mdot_wall)
        mdot_wall_tot = mdot_tot 

        END SUBROUTINE carbon_ablation_park_suzuki_nitri_mdot_wall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE carbon_ablation_park_no_nitri_mdot_wall(rhoi_wall,T_wall,mdot_tot)

        USE mod_general_data,          ONLY: mi, Ri
        USE mod_ablation_data,         ONLY: mdot_wall, mdot_wall_tot,           &
                                           & O_pos, C_pos, O2_pos, N_pos, C3_pos 

        REAL(KIND=8) :: p_eq, rho_eq

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_wall       !< species density vector at the wall
        REAL(KIND=8), INTENT(IN)                :: T_wall          !< wall temperature
        REAL(KIND=8), INTENT(OUT)               :: mdot_tot        !< total char mass blowing rate  at the wall

        !O + C_s => CO
        mdot_wall(1)=rhoi_wall(O_pos)*dsqrt(Ri(O_pos)*T_wall/ &
        &(2.d0*acos(-1.d0)))*0.63d0*dexp(-1160.d0/T_wall)*mi(C_pos)/mi(O_pos)

        !O2 + 2C_s => 2CO
        mdot_wall(2)=rhoi_wall(O2_pos)*dsqrt(2.d0*Ri(O2_pos)*T_wall/(acos(-1.d0)))*0.5d0*mi(C_pos)/mi(O2_pos)

        !N + C_s => CN
        mdot_wall(3)= 0.d0                                                                         ! No nitridation

        !3C_s => C3 
        p_eq=5.19d14*dexp(-90845.d0/T_wall)
        rho_eq=p_eq/(Ri(C3_pos)*T_wall)
        mdot_wall(4)=(rho_eq-rhoi_wall(C3_pos))*dsqrt(Ri(C3_pos)*T_wall/(2.d0*acos(-1.d0)))

        mdot_tot = SUM(mdot_wall)
        mdot_wall_tot = mdot_tot 

        END SUBROUTINE carbon_ablation_park_no_nitri_mdot_wall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE carbon_ablation_user_defined_mdot_wall(rhoi_wall,T_wall,mdot_tot)

        USE mod_general_data,          ONLY: mi, Ri
        USE mod_ablation_data,         ONLY: mdot_wall, mdot_wall_tot,           &
                                           & O_pos, C_pos, O2_pos, N_pos, C3_pos, &
                                           & O_prob, O2_prob, N_prob, C_prob

        REAL(KIND=8) :: p_eq, rho_eq

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_wall       !< species density vector at the wall
        REAL(KIND=8), INTENT(IN)                :: T_wall          !< wall temperature
        REAL(KIND=8), INTENT(OUT)               :: mdot_tot        !< total char mass blowing rate  at the wall

        !O + C_s => CO
        mdot_wall(1)=rhoi_wall(O_pos)*dsqrt(Ri(O_pos)*T_wall/ &
        &(2.d0*acos(-1.d0)))*O_prob*dexp(-1160.d0/T_wall)*mi(C_pos)/mi(O_pos)

        !O2 + 2C_s => 2CO
        mdot_wall(2)=rhoi_wall(O2_pos)*dsqrt(2.d0*Ri(O2_pos)*T_wall/(acos(-1.d0)))*O2_prob*mi(C_pos)/mi(O2_pos)

        !N + C_s => CN
        mdot_wall(3)=rhoi_wall(N_pos)*dsqrt(Ri(N_pos)*T_wall/(2.d0*acos(-1.d0)))*N_prob*mi(C_pos)/mi(N_pos) 

        !3C_s => C3 
        p_eq=5.19d14*dexp(-90845.d0/T_wall)
        rho_eq=p_eq/(Ri(C3_pos)*T_wall)
        mdot_wall(4)=(rho_eq-rhoi_wall(C3_pos))*dsqrt(Ri(C3_pos)*T_wall/(2.d0*acos(-1.d0)))*C_prob

        mdot_tot = SUM(mdot_wall)
        mdot_wall_tot = mdot_tot 

        END SUBROUTINE carbon_ablation_user_defined_mdot_wall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE carbon_nitridation_user_defined_mdot_wall(rhoi_wall,T_wall,mdot_tot)

        USE mod_general_data,          ONLY: mi, Ri
        USE mod_ablation_data,         ONLY: mdot_wall, mdot_wall_tot,           &
                                           & C_pos, N_pos, C3_pos, N_prob, C_prob

        REAL(KIND=8) :: p_eq, rho_eq

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_wall       !< species density vector at the wall
        REAL(KIND=8), INTENT(IN)                :: T_wall          !< wall temperature
        REAL(KIND=8), INTENT(OUT)               :: mdot_tot        !< total char mass blowing rate  at the wall

        !The first two element of the mdot_wall vector a null because there is not oxidation (note that they are not used at all in the code in this case)
        mdot_wall(1)=0.d0
        mdot_wall(2)=0.d0

        !N + C_s => CN
        mdot_wall(3)=rhoi_wall(N_pos)*dsqrt(Ri(N_pos)*T_wall/(2.d0*acos(-1.d0)))*N_prob*mi(C_pos)/mi(N_pos) 

        !3C_s => C3 
        p_eq=5.19d14*dexp(-90845.d0/T_wall)
        rho_eq=p_eq/(Ri(C3_pos)*T_wall)
        mdot_wall(4)=(rho_eq-rhoi_wall(C3_pos))*dsqrt(Ri(C3_pos)*T_wall/(2.d0*acos(-1.d0)))*C_prob

        mdot_tot = SUM(mdot_wall)
        mdot_wall_tot = mdot_tot 

        END SUBROUTINE carbon_nitridation_user_defined_mdot_wall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE carbon_nitridation_mdot_wall_imposed(rhoi_wall,T_wall,mdot_tot)

        USE mod_general_data,          ONLY: mi, Ri
        USE mod_ablation_data,         ONLY: mdot_wall, mdot_wall_tot,           &
                                           & C_pos, N_pos, C3_pos, N_prob, C_prob, exp_mdot

        REAL(KIND=8) :: inv 

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_wall       !< species density vector at the wall
        REAL(KIND=8), INTENT(IN)                :: T_wall          !< wall temperature
        REAL(KIND=8), INTENT(OUT)               :: mdot_tot        !< total char mass blowing rate  at the wall

        !The first two element of the mdot_wall vector a null because there is not oxidation (note that they are not used at all in the code in this case)
        mdot_wall(1)=0.d0
        mdot_wall(2)=0.0

        !N + C_s => CN
           ! The experimental value of the mass blowing rate is used.
        mdot_wall(3)=exp_mdot
           ! The reaction probability that generates this mass blowing is evaluated and stored 
        inv=1.d0/(rhoi_wall(N_pos)*dsqrt(Ri(N_pos)*T_wall/(2.d0*acos(-1.d0)))*mi(C_pos)/mi(N_pos))
        N_prob=mdot_wall(3)*inv

        !The sublimation is null because only surface nitridation has a significant role  (this has to be verified with a previous simulation using
        ! the subroutine "carbon_nitridation_user_defined_mdot_wall" to show that the sublimation is negligible in the analyzed conditions)
        mdot_wall(4)=0.d0

        mdot_tot = SUM(mdot_wall)
        mdot_wall_tot = mdot_tot 

        END SUBROUTINE carbon_nitridation_mdot_wall_imposed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE carbon_oxidation_park_rebuilt_nitri_mdot_wall(rhoi_wall,T_wall,mdot_tot)

        USE mod_general_data,          ONLY: mi, Ri
        USE mod_ablation_data,         ONLY: mdot_wall, mdot_wall_tot,                   &
                                           & O_pos, C_pos, O2_pos, N_pos, C3_pos,        &
                                           & N_prob, gamma_rec

        REAL(KIND=8) :: p_eq, rho_eq

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_wall       !< species density vector at the wall
        REAL(KIND=8), INTENT(IN)                :: T_wall          !< wall temperature
        REAL(KIND=8), INTENT(OUT)               :: mdot_tot        !< total char mass blowing rate  at the wall

        !O + C_s => CO
        mdot_wall(1)=rhoi_wall(O_pos)*dsqrt(Ri(O_pos)*T_wall/ &
        &(2.d0*acos(-1.d0)))*0.63d0*dexp(-1160.d0/T_wall)*mi(C_pos)/mi(O_pos)

        !O2 + 2C_s => 2CO
        mdot_wall(2)=rhoi_wall(O2_pos)*dsqrt(2.d0*Ri(O2_pos)*T_wall/(acos(-1.d0)))*0.5d0*mi(C_pos)/mi(O2_pos)

        !N + C_s => CN 
         ! The nitridation reaction probability is calculated using the linear
         ! interpolation of the value obtained in the pure nitrogen test rebuilding

        ! Olynick reactions w/ air NDPs  for P4 (FIRST_NDPs) and measured test static pressure
        ! in the rebuilding for P4. VKI reference value for copper catalycity
        ! 0.1 @ 1500Pa for P3 AND P4 
        !N_prob=(T_wall-2030.d0)*(0.0039d0-0.0014d0)/(2644.d0-2030.d0) + 0.0014d0


        ! Olynick reactions w/ air NDPs  for P4 (FIRST_NDPs) and updated static pressure
        ! in the rebuilding for P4. VKI reference value for copper catalycity
        ! 0.1 @ 1500Pa for P3 AND P4 
        !N_prob=(T_wall-2030.d0)*(0.0042d0-0.0015d0)/(2644.d0-2030.d0) + 0.0015d0


        ! Olynick reactions w/ nitrogen NDPs  for P4 (SECOND_NDPs) and updated static pressure
        ! in the rebuilding for P4. VKI reference value for copper catalycity
        ! 0.1 @ 1500Pa for P3 AND P4 
        ! If linear interpolation is used
        !N_prob=(T_wall-2029.3d0)*(0.00451d0-0.00152d0)/(2643.2d0-2029.3d0) + 0.00152d0
        ! If Arrhenius law interpolation is used
        !N_prob=0.00152d0/dexp(2643.2d0*dlog(0.00451d0/0.00152d0)/(2029.3d0-2643.2d0))*    &
        !       & dexp(2643.2d0*2029.3d0*dlog(0.00451d0/0.00152d0)/(2029.3d0-2643.2d0)/T_wall)


        ! Dunn-Kang reactions w/ nitrogen NDPs  for P4 (SECOND_NDPs) and updated static pressure
        ! in the rebuilding for P4. VKI reference value for copper catalycity,
        ! 0.1 @ 1500Pa for P3 AND P4 
        !N_prob=(T_wall-2029.3d0)*(0.00455d0-0.00152d0)/(2643.6d0-2029.3d0) + 0.00152d0


        ! Olynick reactions w/ nitrogen NDPs for P4 (SECOND_NDPs), updated static pressure
        ! and copper catalycicty 0.3 (Driver 2013) in the rebuilding for P3 and P4
        ! If Arrhenius law interpolation is used
        !N_prob=0.813d-3/dexp(2644.d0*dlog(0.339d-2/0.813d-3)/(2030.d0-2644.d0))*    &
        !       & dexp(2644.d0*2030.d0*dlog(0.339d-2/0.813d-3)/(2030.d0-2644.d0)/T_wall)

        ! Olynick reactions w/ nitrogen NDPs  for P4 (SECOND_NDPs), updated static pressure
        ! and copper catalycicty 0.19 (giving the right surface temperature of
        ! P4 w/o any preform catalycity for N->N2) for P3 and P4.
        ! If Arrhenius law interpolation is used
        N_prob=0.6534440709d-3/dexp(2636.d0*dlog(0.3200876880d-2/0.6534440709d-3)/(2048.d0-2636.d0))*    &
               & dexp(2636.d0*2048.d0*dlog(0.3200876880d-2/0.6534440709d-3)/(2048.d0-2636.d0)/T_wall)

        mdot_wall(3)=rhoi_wall(N_pos)*dsqrt(Ri(N_pos)*T_wall/(2.d0*acos(-1.d0)))*N_prob*mi(C_pos)/mi(N_pos) 

        !3C_s => C3 
        p_eq=5.19d14*dexp(-90845.d0/T_wall)
        rho_eq=p_eq/(Ri(C3_pos)*T_wall)
        mdot_wall(4)=(rho_eq-rhoi_wall(C3_pos))*dsqrt(Ri(C3_pos)*T_wall/(2.d0*acos(-1.d0)))

        mdot_tot = SUM(mdot_wall)
        mdot_wall_tot = mdot_tot 

        ! Also the recombination probability is computed (and stored) from the linear interpolation
        ! of the pure nitrogen test

        ! Olynick reactions w/ air NDPs  for P4 (FIRST_NDPs) and measured test static pressure
        ! in the rebuilding for P4. VKI reference value for copper catalycity
        ! 0.1 @ 1500Pa for P3 AND P4 
        !gamma_rec=(T_wall-2030.d0)*(0.05414d0-0.08015d0)/(2644.d0-2030.d0) + 0.08015d0

        ! Olynick reactions w/ air NDPs  for P4 (FIRST_NDPs) and updated static pressure
        ! in the rebuilding for P4. VKI reference value for copper catalycity
        ! 0.1 @ 1500Pa for P3 AND P4 
        !gamma_rec=(T_wall-2030.d0)*(0.06974d0-0.08515d0)/(2644.d0-2030.d0) + 0.08515d0

        ! Olynick reactions w/ nitrogen NDPs  for P4 (SECOND_NDPs) and updated static pressure
        ! in the rebuilding for P4. VKI reference value for copper catalycity
        ! 0.1 @ 1500Pa for P3 AND P4 
        ! If linear interpolation is used
        !gamma_rec=(T_wall-2029.3d0)*(0.08086d0-0.08515d0)/(2643.2d0-2029.3d0) + 0.08515d0
        ! If Arrhenius law interpolation is used
        !gamma_rec=0.08515d0/dexp(2643.2d0*dlog(0.08086d0/0.08515d0)/(2029.3d0-2643.2d0))*    &
        !       & dexp(2643.2d0*2029.3d0*dlog(0.08086d0/0.08515d0)/(2029.3d0-2643.2d0)/T_wall)

        ! Dunn-Kang reactions w/ nitrogen NDPs  for P4 (SECOND_NDPs) and updated static pressure
        ! in the rebuilding for P4. VKI reference value for copper catalycity,
        ! 0.1 @ 1500Pa for P3 AND P4 
        !gamma_rec=(T_wall-2029.3d0)*(0.09186d0-0.08515d0)/(2643.6d0-2029.3d0) + 0.08515d0

        ! Olynick reactions w/ nitrogen NDPs for P4 (SECOND_NDPs), updated static pressure
        ! and copper catalycicty 0.3 (Driver 2013) in the rebuilding for P3 and P4
        ! If T_wall~2030
        ! gamma_rec=0.04d0
        ! If T_wall~2644
        ! gamma_rec=0.01610d0


        ! Olynick reactions w/ nitrogen NDPs  for P4 (SECOND_NDPs) and updated static pressure
        ! in the rebuilding for P4. VKI reference value for copper catalycity
        ! 0.1 @ 1500Pa for P3 AND P4 
        ! If T_wall~2030
        !gamma_rec=0.08515d0 
        ! If T_wall~2644
        !gamma_rec=0.08086d0

        ! Olynick reactions w/ nitrogen NDPs  for P4 (SECOND_NDPs), updated static pressure
        ! and copper catalycicty 0.19 (giving the right surface temperature of
        ! P4 w/o any preform catalycity for N->N2) for P3 and P4.
        ! If T_wall~2030
        gamma_rec=0.03d0
        ! If T_wall~2644
        !gamma_rec=0.d0

        END SUBROUTINE carbon_oxidation_park_rebuilt_nitri_mdot_wall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif  
  END MODULE mod_ablation_procedure_carbon_library
