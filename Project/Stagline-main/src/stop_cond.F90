!--------------------------------------------------------------------------!
!> This function evaluates the stop condition according to what specified in the input file by the user. 
  FUNCTION ev_stop_cond (it, level, res)
  
    USE mod_general_data,       ONLY: nb_eq, nb_ns
    USE mod_numerics_data,      ONLY: nb_itmax, time_limit, stop_id,lim_res

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: it                         !< iteration number 
    REAL(KIND=8), INTENT(IN) :: level                 !< time-level
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: res     !< rhs residual on conservative variables
    LOGICAL :: ev_stop_cond
  
    ev_stop_cond = .FALSE.
    
    ! Selection of test on terminate condition
    SELECT CASE (stop_id)
    
      ! Condition on the residual value
      CASE (1) 
        IF (MAXVAL(res).LE.lim_res) ev_stop_cond = .TRUE.        
    
      ! Condition on the time level
      CASE (2) 
        IF (level.GE.time_limit) ev_stop_cond = .TRUE.  
    
      ! Condition on the number of interation
      CASE (3) 
        IF (it.GE.nb_itmax) ev_stop_cond = .TRUE.  
    
       CASE (4) 
        IF ((it.GE.nb_itmax) .OR. (MAXVAL(res).LE.lim_res)) Then 
          ev_stop_cond = .TRUE. 
       ENDIF  
    
    END SELECT 
   
    WRITE(*,5)it,dlog10(sum(10.0d0**res(1:nb_ns))),res(nb_ns + 1:nb_eq)

5 FORMAT(I8,1X,100F10.5)
   
  END FUNCTION ev_stop_cond  
!------------------------------------------------------------------------------!

