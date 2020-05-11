C It sets up all C bindings for interfacing to the fortran Matrix Element

C Additional interface subroutines to add:

      SUBROUTINE %(binding_prefix)sGET_NUMERATOR(P,PROC_ID,SELECTED_DIAGRAM_LEFT,SELECTED_DIAGRAM_RIGHT,ANS_RE,ANS_IM) bind(c)
        USE iso_c_binding, ONLY: c_int, c_double
        INTEGER I,J

        INTEGER MAX_PROC_NUMBER
        PARAMETER( MAX_PROC_NUMBER=%(max_proc_number)d)
        INTEGER PROC_ID_TO_POSITION(0:MAX_PROC_NUMBER)
        DATA (PROC_ID_TO_POSITION(I),I=0,MAX_PROC_NUMBER)/%(proc_number_to_position_map)s/

        %(truncated_mom_list)s

C       Format ia P1E_re, P1E_im, P1x_re, P1x_im, ..., P2E_re, P2E_im, ...
        real(c_double), intent(in)  :: P(*)
        integer(c_int), intent(in)  :: PROC_ID
        integer(c_int), intent(in)  :: SELECTED_DIAGRAM_LEFT
        integer(c_int), intent(in)  :: SELECTED_DIAGRAM_RIGHT
        real(c_double), intent(out) :: ANS_RE
        real(c_double), intent(out) :: ANS_IM
        double complex ANS
        INTEGER SDL, SDG

        SDL = SELECTED_DIAGRAM_LEFT
        SDG = SELECTED_DIAGRAM_RIGHT

        GOTO (%(proc_positions_list)s) PROC_ID_TO_POSITION(PROC_ID)

        %(matrix_element_call_dispatch)s

        ANS_RE = DBLE(ANS)
        ANS_IM = DIMAG(ANS)

9999  continue

      END SUBROUTINE %(binding_prefix)sGET_NUMERATOR

      SUBROUTINE %(binding_prefix)sINITIALISE(PATH) bind(c)
        USE iso_c_binding
        CHARACTER(c_char) :: PATH(512)
        CHARACTER(512) :: PATH_IN
        DO i=1,512
          PATH_IN(i:i) = PATH(i)
        ENDDO
        CALL %(first_proc_prefix)sINITIALISEMODEL(PATH_IN)
      END SUBROUTINE %(binding_prefix)sINITIALISE
