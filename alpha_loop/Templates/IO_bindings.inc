      PROGRAM IO_binding

      implicit none
        
      integer i,j

      double precision P(0:7)
      double complex ANS

      %(truncated_mom_list)s

      INTEGER MAX_PROC_NUMBER
      PARAMETER( MAX_PROC_NUMBER=%(max_proc_number)d)
      INTEGER PROC_ID_TO_POSITION(0:MAX_PROC_NUMBER)
      DATA (PROC_ID_TO_POSITION(I),I=0,MAX_PROC_NUMBER)/%(proc_number_to_position_map)s/

      integer proc_ID
      integer SDL, SDR

      CHARACTER(512) PATH

      READ(*,'(A)') PATH

      IF (PATH(1:1).eq.'') THEN
          CALL SETPARA('param_card.dat')
      ELSE
          CALL SETPARA2(PATH)
      ENDIF 

9997  continue
      READ(*,*) PROC_ID, SDL, SDR
      IF (PROC_ID .lt. 0) goto 9998

      GOTO (%(proc_positions_list)s) PROC_ID_TO_POSITION(PROC_ID)

      %(matrix_element_call_dispatch)s

9999  CONTINUE
      WRITE(*,*) DBLE(ANS), DIMAG(ANS)
      GOTO 9997
9998  CONTINUE

      END PROGRAM IO_binding

