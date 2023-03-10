      SUBROUTINE %(proc_prefix)sSMATRIXHEL(P,HEL, SELECTED_DIAG_LEFT,SELECTED_DIAG_RIGHT, ANS)
      IMPLICIT NONE
C
C CONSTANT
C
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=%(ncomb)d)
CF2PY INTENT(OUT) :: ANS                                                                                                                                                                                     
CF2PY INTENT(IN) :: HEL                                                                                                                                                                                    
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL)  

C  
C ARGUMENTS 
C  
      DOUBLE COMPLEX P(0:3,NEXTERNAL),ANS
	INTEGER HEL, SELECTED_DIAG_LEFT,SELECTED_DIAG_RIGHT
C
C GLOBAL VARIABLES
C
      INTEGER USERHEL
	  COMMON/%(proc_prefix)sHELUSERCHOICE/USERHEL
C ----------
C BEGIN CODE
C ----------
      USERHEL=HEL
	  CALL %(proc_prefix)sSMATRIX(P, SELECTED_DIAG_LEFT,SELECTED_DIAG_RIGHT,ANS)
	  USERHEL=-1

	  END

      SUBROUTINE %(proc_prefix)sSMATRIX(P, SELECTED_DIAG_LEFT,SELECTED_DIAG_RIGHT,ANS)
C  
%(info_lines)s
C 
C MadGraph5_aMC@NLO StandAlone Version
C 
C Returns amplitude squared summed/avg over colors
c and helicities
c for the point in phase space P(0:3,NEXTERNAL)
C  
%(process_lines)s
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
	  INTEGER    NINITIAL 
      PARAMETER (NINITIAL=%(nincoming)d)

      INTEGER NPOLENTRIES
      PARAMETER (NPOLENTRIES=(NEXTERNAL+1)*6)
      INTEGER                 NCOMB         
      PARAMETER (             NCOMB=%(ncomb)d)
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=%(ngraphs)d) 
      INTEGER NGOODTRYENTIRES
      PARAMETER (NGOODTRYENTIRES=NGRAPHS*NGRAPHS*NCOMB) 
      INTEGER NGRAPHSSQ
      PARAMETER (NGRAPHSSQ=NGRAPHS*NGRAPHS) 
	INTEGER HELAVGFACTOR
	PARAMETER (HELAVGFACTOR=%(hel_avg_factor)d)
      
      INTEGER MIN_TRY_FOR_HEL_FILTER
C     Setting below to -1 deactivate the helicity filtering
C     TODO this filter should be precomputed and read from an external file at runtime.
      PARAMETER(MIN_TRY_FOR_HEL_FILTER=-1)
C  
C ARGUMENTS 
C  
      double complex P(0:3,NEXTERNAL),ANS
CF2PY INTENT(OUT) :: ANS
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL)
      INTEGER SELECTED_DIAG_LEFT,SELECTED_DIAG_RIGHT
CF2PY INTENT(IN) :: SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT

C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY(NGRAPHS,NGRAPHS)
C  put in common block to expose this variable to python interface
      COMMON/%(proc_prefix)sPROCESS_NHEL/NHEL 
      double complex T
      double complex %(proc_prefix)sMATRIX
      INTEGER IHEL,IDEN, I, J
C     For a 1>N process, them BEAMTWO_HELAVGFACTOR would be set to 1.
      INTEGER BEAMS_HELAVGFACTOR(2)
	  DATA (BEAMS_HELAVGFACTOR(I),I=1,2)/%(beamone_helavgfactor)d,%(beamtwo_helavgfactor)d/
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NGRAPHS,NGRAPHS,NCOMB)
      DATA NTRY/NGRAPHSSQ*0/
      DATA GOODHEL/NGOODTRYENTIRES*.FALSE./

C
C GLOBAL VARIABLES
C
      INTEGER USERHEL
	  COMMON/%(proc_prefix)sHELUSERCHOICE/USERHEL
	  DATA USERHEL/-1/	  

%(helicity_lines)s
%(den_factor_line)s

	  INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/%(proc_prefix)sBORN_BEAM_POL/POLARIZATIONS
      DATA ((POLARIZATIONS(I,J),I=0,NEXTERNAL),J=0,5)/NPOLENTRIES*-1/

      INTEGER NHEL_COMPUTED

C
C     FUNCTIONS
C
      LOGICAL %(proc_prefix)sIS_BORN_HEL_SELECTED

C ----------
C BEGIN CODE
C ----------
      IF(USERHEL.EQ.-1) NTRY(SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT)=NTRY(SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT)+1
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
C     When spin-2 particles are involved, the Helicity filtering is dangerous for the 2->1 topology.
C     This is because depending on the MC setup the initial PS points have back-to-back initial states
C     for which some of the spin-2 helicity configurations are zero. But they are no longer zero
C     if the point is boosted on the z-axis. Remember that HELAS helicity amplitudes are no longer
C     lorentz invariant with expternal spin-2 particles (only the helicity sum is).
C     For this reason, we simply remove the filterin when there is only three external particles.
C     ALSO DO THIS IF WE WANTED TO DEACTIVATE THE FILTER
      IF ((NTRY(SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT).EQ.1).AND.((NEXTERNAL.LE.3).OR.(MIN_TRY_FOR_HEL_FILTER.lt.0))) THEN
        DO IHEL=1,NCOMB
		GOODHEL(SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT, IHEL)=.TRUE.
	  ENDDO
	ENDIF
      ANS = 0D0
      NHEL_COMPUTED = 0
          DO IHEL=1,NCOMB
	         IF (USERHEL.EQ.-1.OR.USERHEL.EQ.IHEL) THEN		  
                IF (GOODHEL(SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT,IHEL) .OR. (NTRY(SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT) .LT. MIN_TRY_FOR_HEL_FILTER).OR.(USERHEL.NE.-1)) THEN
				    IF(NTRY(SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT).GE.2.AND.POLARIZATIONS(0,0).ne.-1.and.(.not.%(proc_prefix)sIS_BORN_HEL_SELECTED(IHEL))) THEN
                        CYCLE
                    ENDIF
                    NHEL_COMPUTED = NHEL_COMPUTED+1 
                    T=%(proc_prefix)sMATRIX(P ,NHEL(1,IHEL),JC(1), SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT)
                    IF(POLARIZATIONS(0,0).eq.-1.or.%(proc_prefix)sIS_BORN_HEL_SELECTED(IHEL)) THEN
                        ANS=ANS+T
                    ENDIF
				    IF (((DBLE(T) .NE. 0D0).OR.(DIMAG(T) .NE. 0D0)) .AND. .NOT.    GOODHEL(SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT, IHEL)) THEN
                        GOODHEL(SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT, IHEL)=.TRUE.
                    ENDIF
                ENDIF
		     ENDIF
          ENDDO

c      WRITE(*,*) 'N helicity config computed:',NHEL_COMPUTED,'/',NCOMB

C     For LTD squared we must *not* divide by the final state symmetry factor
C      ANS=ANS/DBLE(IDEN)
      ANS=ANS/DBLE(HELAVGFACTOR)

C	  IF(USERHEL.NE.-1) THEN
C          ANS=ANS*HELAVGFACTOR	  
C	  ELSE
C		DO J=1,NINITIAL
C		  IF (POLARIZATIONS(J,0).ne.-1) THEN
C            ANS=ANS*BEAMS_HELAVGFACTOR(J)
C            ANS=ANS/POLARIZATIONS(J,0)
C		  ENDIF
C		ENDDO
C	  ENDIF
	  END
       
      SUBROUTINE %(proc_prefix)sGET_N_DIAGRAMS( NDIAG )

      INTEGER NDIAG
      INTEGER    NDIAGS_IN_THIS_PROC
      PARAMETER (NDIAGS_IN_THIS_PROC=%(ndiags_in_this_proc)d) 
      
      NDIAG = NDIAGS_IN_THIS_PROC

      END SUBROUTINE %(proc_prefix)sGET_N_DIAGRAMS

      DOUBLE COMPLEX FUNCTION %(proc_prefix)sMATRIX(P,NHEL,IC, SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT)
C  
%(info_lines)s
C
C Returns amplitude squared summed/avg over colors
c for the point with external lines W(0:6,NEXTERNAL)
C  
%(process_lines)s
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=%(ngraphs)d) 
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=%(nwavefuncs)d, NCOLOR=%(ncolor)d) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1(2)

      COMPLEX*16 CMPLX_ONE
      PARAMETER (CMPLX_ONE=DCMPLX(1.0d0,0.0d0))

C  
C ARGUMENTS 
C  
      DOUBLE COMPLEX P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL), SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT
C  
C LOCAL VARIABLES 
C  
      INTEGER SELECTED_DIAGS(2)
      INTEGER I,J,K
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(2,NCOLOR)
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1

      DATA DUM0, DUM1/(0d0, 0d0), (1d0, 0d0)/
C  
C GLOBAL VARIABLES
C  
      include 'coupl.inc'
%(global_variable)s
C  
C COLOR DATA
C  
%(color_data_lines)s
C ----------
C BEGIN CODE
C ----------

IMAG1(1)=DCMPLX(0D0,1D0)
IMAG1(2)=DCMPLX(0D0,-1D0)

SELECTED_DIAGS(1) = SELECTED_DIAG_LEFT
SELECTED_DIAGS(2) = SELECTED_DIAG_RIGHT

DO K=1,2
DO I=1,NCOLOR
JAMP(K,I) = DCMPLX(0.0d0, 0.0d0)
ENDDO
DO I=1,NGRAPHS
   AMP(I) = DCMPLX(0.0d0, 0.0d0)
ENDDO
C Jump to the code for the selected diagram
GOTO (%(diagram_indices)s) SELECTED_DIAGS(K) 
%(helas_calls)s
C End of the list of code defining the computation for each diagram
9999 CONTINUE
%(jamp_lines)s
ENDDO

      %(proc_prefix)sMATRIX = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(1,J)
          ENDDO
C Complex conjugation must not be performed here.
C It is already automatically applied to the external incoming polarisation vectors
          %(proc_prefix)sMATRIX = %(proc_prefix)sMATRIX+ZTEMP*JAMP(2,I)/DENOM(I)   
      ENDDO
      END

      SUBROUTINE %(proc_prefix)sGET_value(P, ALPHAS, NHEL, SELECTED_DIAG_LEFT,SELECTED_DIAG_RIGHT ,ANS)
      IMPLICIT NONE   
C
C CONSTANT
C
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
C  
C ARGUMENTS 
C  
      DOUBLE COMPLEX P(0:3,NEXTERNAL),ANS
      INTEGER NHEL, SELECTED_DIAG_LEFT,SELECTED_DIAG_RIGHT
      DOUBLE PRECISION ALPHAS 
      real*8 pi     
CF2PY INTENT(OUT) :: ANS  
CF2PY INTENT(IN) :: NHEL   
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL) 
CF2PY INTENT(IN) :: ALPHAS
C     ROUTINE FOR F2PY to read the benchmark point.    
C     the include file with the values of the parameters and masses 
      include "coupl.inc"
      
      pi = 3.141592653589793d0
      G = 2* DSQRT(ALPHAS*pi)
      call UPDATE_AS_PARAM()
      if (NHEL.ne.0) then
        CALL %(proc_prefix)sSMATRIXHEL(P, NHEL, SELECTED_DIAG_LEFT,SELECTED_DIAG_RIGHT, ANS)
      else
        CALL %(proc_prefix)sSMATRIX(P, SELECTED_DIAG_LEFT, SELECTED_DIAG_RIGHT, ANS)
      endif
      return 
      end

      SUBROUTINE %(proc_prefix)sINITIALISEMODEL(PATH)
C     ROUTINE FOR F2PY to read the benchmark point.    
      IMPLICIT NONE   
      CHARACTER*512 PATH
CF2PY INTENT(IN) :: PATH 
C     USE SETPARA2 and not SETPARA so that ident_card.dat can be
C      automatically found.
      call setpara2(PATH)  !first call to setup the paramaters    
      return 
      end      

      LOGICAL FUNCTION %(proc_prefix)sIS_BORN_HEL_SELECTED(HELID)
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=%(nexternal)d)
      INTEGER    NCOMB
      PARAMETER (NCOMB=%(ncomb)d)
C
C     ARGUMENTS
C
      INTEGER HELID
C
C     LOCALS
C
      INTEGER I,J
      LOGICAL FOUNDIT
C
C     GLOBALS
C
      INTEGER HELC(NEXTERNAL,NCOMB)
      COMMON/%(proc_prefix)sPROCESS_NHEL/HELC

      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/%(proc_prefix)sBORN_BEAM_POL/POLARIZATIONS
C     ----------
C     BEGIN CODE
C     ----------
      
      %(proc_prefix)sIS_BORN_HEL_SELECTED = .True.
      if (POLARIZATIONS(0,0).eq.-1) THEN
        RETURN
      ENDIF

      DO I=1,NEXTERNAL
        IF (POLARIZATIONS(I,0).eq.-1) THEN
          CYCLE
        ENDIF
        FOUNDIT = .FALSE.
        DO J=1,POLARIZATIONS(I,0)
          IF (HELC(I,HELID).eq.POLARIZATIONS(I,J)) THEN
            FOUNDIT = .True.
            EXIT
          ENDIF
        ENDDO
        IF(.not.FOUNDIT) THEN
          %(proc_prefix)sIS_BORN_HEL_SELECTED = .False.
          RETURN
        ENDIF
      ENDDO

      RETURN
      END
