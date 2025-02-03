!------------------------------------------------------------------------------!
!> This subroutine reads the CFL number when changed interactively during the simulation or
!> it takes the CFL from the History log FIle
   SUBROUTINE update_cfl()

     USE mod_general_data,        ONLY: cfl_file, log_cfl_int, log_cfl_number, inter_cfl, &
                                        & log_cfl_file_lines
     USE mod_numerics_data,       ONLY: cfl

     IMPLICIT NONE

     INTEGER, PARAMETER :: in1 = 1
     REAL(KIND=8), SAVE :: cfl_old

     ! Initialization
     cfl_old = cfl
    
     OPEN(UNIT=in1,FILE=cfl_file(1:LEN_TRIM(cfl_file)),STATUS='unknown')

     IF (inter_cfl.EQV..TRUE.) THEN 
       ! Read interactive CFL number
     READ(in1,*)cfl
     ELSE
       cfl = log_cfl_number(log_cfl_int)
       log_cfl_int = MIN(log_cfl_int + 1,log_cfl_file_lines)
     ENDIF

     CLOSE(in1)

     IF (cfl.NE.cfl_old) THEN
       WRITE(*,5)'solver_fvmcc_F90:: CFL updated'
       WRITE(*,10)'CFL = ',cfl
     ENDIF

5  FORMAT(A)
10 FORMAT(A,E14.6)

   END SUBROUTINE update_cfl
!------------------------------------------------------------------------------!
