      program PolSumTest
C**************************************************************************
C     Program testing polarisation sums
C**************************************************************************
      implicit none
C     
C     constants  
C     
      double precision rzero
      parameter (rzero=0.0d0)

      double complex czero
      parameter (czero=DCMPLX(0.0d0,0.0d0))

C     
C     LOCAL
C     
      integer i,j, k

      double complex pTest(0:3)      
      double complex pTestRealOnshell(0:3)
      double complex pTestRealOffshell(0:3)
      double complex pTestCmplxOnshell(0:3)
      double complex pTestCmplxOffshell(0:3)

      double precision test_mass
      logical did_test_pass, verbosity

C-----
C     BEGIN CODE
C-----

      pTestRealOnshell(1) = DCMPLX(1.0d0,0.0d0)
      pTestRealOnshell(2) = DCMPLX(2.0d0,0.0d0)
      pTestRealOnshell(3) = DCMPLX(3.0d0,0.0d0)
      pTestRealOnshell(0) = sqrt(
     &      pTestRealOnshell(1)**2+
     &      pTestRealOnshell(2)**2+
     &      pTestRealOnshell(3)**2
     & )

      pTestRealOffshell(1) = DCMPLX(1.0d0,0.0d0)
      pTestRealOffshell(2) = DCMPLX(2.0d0,0.0d0)
      pTestRealOffshell(3) = DCMPLX(3.0d0,0.0d0)
      pTestRealOffshell(0) = sqrt(
     &      pTestRealOffshell(1)**2+
     &      pTestRealOffshell(2)**2+
     &      pTestRealOffshell(3)**2
     & ) + DCMPLX(7.0,0.0d0)

      pTestCmplxOnshell(1) = DCMPLX(1.0d0,4.0d0)
      pTestCmplxOnshell(2) = DCMPLX(3.0d0,5.0d0)
      pTestCmplxOnshell(3) = DCMPLX(7.0d0,2.0d0)
      pTestCmplxOnshell(0) = sqrt(
     &      pTestCmplxOnshell(1)**2+
     &      pTestCmplxOnshell(2)**2+
     &      pTestCmplxOnshell(3)**2
     & )

      pTestCmplxOffshell(1) = DCMPLX(1.0d0,4.0d0)
      pTestCmplxOffshell(2) = DCMPLX(3.0d0,5.0d0)
      pTestCmplxOffshell(3) = DCMPLX(7.0d0,2.0d0)
      pTestCmplxOffshell(0) = sqrt(
     &      pTestCmplxOffshell(1)**2+
     &      pTestCmplxOffshell(2)**2+
     &      pTestCmplxOffshell(3)**2
     & ) + DCMPLX(3.0d0, 11.0d0)

      verbosity = .False.

      Write(*,*) ''
      Write(*,*) '>>>> Now testing a massless on-shell real spinor.'
      WRITE(*,*) '================================================='
      test_mass = rzero
      do i=0,3
        pTest(i) = pTestRealOnshell(i)
      enddo
      pTest(0) = sqrt(pTest(0)**2+DCMPLX(test_mass**2,0.0d0))
      WRITE(*,*) 'P=',(pTest(i),i=0,3)
      WRITE(*,*) 'sqrt(P^2)=',sqrt(pTest(0)**2-pTest(1)**2-pTest(2)**2-pTest(3)**2)
      WRITE(*,*) 'mass=',test_mass      
      call testCase(pTest,test_mass, verbosity,did_test_pass)
      if (did_test_pass) then
          WRITE(*,*) '-> PASSED!'
      else
          WRITE(*,*) '-> FAILED!'          
      endif

      Write(*,*) ''
      Write(*,*) '>>>> Now testing a massive on-shell real spinor.'
      WRITE(*,*) '================================================='
      test_mass = 3.5d0
      do i=0,3
        pTest(i) = pTestRealOnshell(i)
      enddo
      pTest(0) = sqrt(pTest(0)**2+DCMPLX(test_mass**2,0.0d0))
      WRITE(*,*) 'P=',(pTest(i),i=0,3)
      WRITE(*,*) 'sqrt(P^2)=',sqrt(pTest(0)**2-pTest(1)**2-pTest(2)**2-pTest(3)**2)
      WRITE(*,*) 'mass=',test_mass
      call testCase(pTest,test_mass, verbosity,did_test_pass)
      if (did_test_pass) then
          WRITE(*,*) '-> PASSED!'
      else
          WRITE(*,*) '-> FAILED!'          
      endif

      Write(*,*) ''      
      Write(*,*) '>>>> Now testing a massless off-shell real spinor.'
      WRITE(*,*) '===================================================='
      test_mass = rzero
      do i=0,3
        pTest(i) = pTestRealOffshell(i)
      enddo
c      pTest(0) = sqrt(pTest(0)**2+DCMPLX(test_mass**2,0.0d0))
      WRITE(*,*) 'P=',(pTest(i),i=0,3)
      WRITE(*,*) 'sqrt(P^2)=',sqrt(pTest(0)**2-pTest(1)**2-pTest(2)**2-pTest(3)**2)
      WRITE(*,*) 'mass=',test_mass
      call testCase(pTest,test_mass, verbosity,did_test_pass)
      if (did_test_pass) then
          WRITE(*,*) '-> PASSED!'
      else
          WRITE(*,*) '-> FAILED!'          
      endif

      Write(*,*) ''      
      Write(*,*) '>>>> Now testing a massive off-shell real spinor.'
      WRITE(*,*) '===================================================='
      test_mass = 3.5d0
      do i=0,3
        pTest(i) = pTestRealOffshell(i)
      enddo
c      pTest(0) = sqrt(pTest(0)**2+DCMPLX(test_mass**2,0.0d0))
      WRITE(*,*) 'P=',(pTest(i),i=0,3)
      WRITE(*,*) 'sqrt(P^2)=',sqrt(pTest(0)**2-pTest(1)**2-pTest(2)**2-pTest(3)**2)
      WRITE(*,*) 'mass=',test_mass
      call testCase(pTest,test_mass, verbosity,did_test_pass)
      if (did_test_pass) then
          WRITE(*,*) '-> PASSED!'
      else
          WRITE(*,*) '-> FAILED!'          
      endif


      Write(*,*) ''
      Write(*,*) '>>>> Now testing a massless on-shell complex spinor.'
      WRITE(*,*) '================================================='
      test_mass = rzero
      do i=0,3
        pTest(i) = pTestCmplxOnshell(i)
      enddo
      pTest(0) = sqrt(pTest(0)**2+DCMPLX(test_mass**2,0.0d0))
      WRITE(*,*) 'P=',(pTest(i),i=0,3)
      WRITE(*,*) 'sqrt(P^2)=',sqrt(pTest(0)**2-pTest(1)**2-pTest(2)**2-pTest(3)**2)
      WRITE(*,*) 'mass=',test_mass      
      call testCase(pTest,test_mass, verbosity,did_test_pass)
      if (did_test_pass) then
          WRITE(*,*) '-> PASSED!'
      else
          WRITE(*,*) '-> FAILED!'          
      endif

      Write(*,*) ''
      Write(*,*) '>>>> Now testing a massive on-shell complex spinor.'
      WRITE(*,*) '================================================='
      test_mass = 3.5d0
      do i=0,3
        pTest(i) = pTestCmplxOnshell(i)
      enddo
      pTest(0) = sqrt(pTest(0)**2+DCMPLX(test_mass**2,0.0d0))
      WRITE(*,*) 'P=',(pTest(i),i=0,3)
      WRITE(*,*) 'sqrt(P^2)=',sqrt(pTest(0)**2-pTest(1)**2-pTest(2)**2-pTest(3)**2)
      WRITE(*,*) 'mass=',test_mass
      call testCase(pTest,test_mass, verbosity,did_test_pass)
      if (did_test_pass) then
          WRITE(*,*) '-> PASSED!'
      else
          WRITE(*,*) '-> FAILED!'          
      endif

      Write(*,*) ''      
      Write(*,*) '>>>> Now testing a massless off-shell complex spinor.'
      WRITE(*,*) '===================================================='
      test_mass = rzero
      do i=0,3
        pTest(i) = pTestCmplxOffshell(i)
      enddo
c      pTest(0) = sqrt(pTest(0)**2+DCMPLX(test_mass**2,0.0d0))
      WRITE(*,*) 'P=',(pTest(i),i=0,3)
      WRITE(*,*) 'sqrt(P^2)=',sqrt(pTest(0)**2-pTest(1)**2-pTest(2)**2-pTest(3)**2)
      WRITE(*,*) 'mass=',test_mass
      call testCase(pTest,test_mass, verbosity,did_test_pass)
      if (did_test_pass) then
          WRITE(*,*) '-> PASSED!'
      else
          WRITE(*,*) '-> FAILED!'          
      endif

      Write(*,*) ''      
      Write(*,*) '>>>> Now testing a massive off-shell complex spinor.'
      WRITE(*,*) '===================================================='
      test_mass = 3.5d0
      do i=0,3
        pTest(i) = pTestCmplxOffshell(i)
      enddo
c      pTest(0) = sqrt(pTest(0)**2+DCMPLX(test_mass**2,0.0d0))
      WRITE(*,*) 'P=',(pTest(i),i=0,3)
      WRITE(*,*) 'sqrt(P^2)=',sqrt(pTest(0)**2-pTest(1)**2-pTest(2)**2-pTest(3)**2)
      WRITE(*,*) 'mass=',test_mass
      call testCase(pTest,test_mass, verbosity,did_test_pass)
      if (did_test_pass) then
          WRITE(*,*) '-> PASSED!'
      else
          WRITE(*,*) '-> FAILED!'          
      endif

      end program PolSumTest

      subroutine testCase(p,mass, verbose, test_passed)

          implicit none

          double complex p(0:3), pE, px, py, pz
          double precision preal(0:3)
          double precision mass

          integer i,j,k,ihel
          integer nhel(2)

          integer left, right

c         Particle
          double complex i_uspinor(8)
          double complex i_ubarspinor(8),i_ubarspinorg0(8)
          double complex o_uspinor(8)
          double complex o_ubarspinor(8),o_ubarspinorg0(8)

c         AntiParticle
          double complex o_vspinor(8)
          double complex o_vbarspinor(8),o_vbarspinorg0(8)
          double complex i_vspinor(8)
          double complex i_vbarspinor(8),i_vbarspinorg0(8)

          double precision tolerance

          double complex i_polsumU(0:3,0:3)
          double complex i_polsumV(0:3,0:3)       
          double complex o_polsumU(0:3,0:3)
          double complex o_polsumV(0:3,0:3)
          double complex i_target_polsumU(0:3,0:3)
          double complex i_target_polsumV(0:3,0:3)
          double complex o_target_polsumU(0:3,0:3)
          double complex o_target_polsumV(0:3,0:3)

          logical passTest(0:3,0:3)
          logical verbose, test_passed, did_printout

          logical use_original_routines
          logical use_mod_routines
          logical use_trick_routines

          double precision rzero
          parameter (rzero=0.0d0)

          double complex czero
          parameter (czero=DCMPLX(0.0d0,0.0d0))

          double complex ci
          parameter (ci=DCMPLX(0.0d0,1.0d0))

c         Left / Right of Cutkosky cut.
          left = 1
          right = 2


          do i=0,3
            do j=0,3
              i_polsumU(i,j) = czero
              i_polsumV(i,j) = czero            
              o_polsumU(i,j) = czero
              o_polsumV(i,j) = czero
            enddo
          enddo

          nhel(1) = -1
          nhel(2) = +1

          use_original_routines = .False.
          use_mod_routines = .False.
          use_trick_routines = .True.
       
          do i=0,3
            preal(i) = DBLE(p(i))
          enddo

          if (.not. use_trick_routines) then
          do ihel=1,2
            IF (use_original_routines) then
c             Particle
              call ixxxxx(preal, mass, nhel(ihel), +1, i_uspinor)
              call ixxxxx(preal, mass, nhel(ihel), +1, i_ubarspinor)
              call oxxxxx(preal, mass, nhel(ihel), +1, o_uspinor)
              call oxxxxx(preal, mass, nhel(ihel), +1, o_ubarspinor)
c             AntiParticle
              call ixxxxx(preal, mass, nhel(ihel), -1, i_vspinor)
              call ixxxxx(preal, mass, nhel(ihel), -1, i_vbarspinor)
              call oxxxxx(preal, mass, nhel(ihel), -1, o_vspinor)
              call oxxxxx(preal, mass, nhel(ihel), -1, o_vbarspinor)
            ELSE IF (use_mod_routines) then
c             Particle
              call PROPPOL_ixxxxx(p, mass, nhel(ihel), +1, left, i_uspinor)
              call PROPPOL_ixxxxx(p, mass, nhel(ihel), +1, right, i_ubarspinor)
              call PROPPOL_oxxxxx(p, mass, nhel(ihel), +1, left, o_uspinor)
              call PROPPOL_oxxxxx(p, mass, nhel(ihel), +1, right, o_ubarspinor)
c             AntiParticle
              call PROPPOL_ixxxxx(p, mass, nhel(ihel), -1, left, i_vspinor)
              call PROPPOL_ixxxxx(p, mass, nhel(ihel), -1, right, i_vbarspinor)
              call PROPPOL_oxxxxx(p, mass, nhel(ihel), -1, left, o_vspinor)
              call PROPPOL_oxxxxx(p, mass, nhel(ihel), -1, right, o_vbarspinor)
            ENDIF
C           We need to multipy by Gamma^0 the bar spinors now

            i_ubarspinorg0(5) = i_ubarspinor(7)
            i_ubarspinorg0(6) = i_ubarspinor(8)
            i_ubarspinorg0(7) = i_ubarspinor(5)
            i_ubarspinorg0(8) = i_ubarspinor(6)
            
            i_vbarspinorg0(5) = i_vbarspinor(7)
            i_vbarspinorg0(6) = i_vbarspinor(8)
            i_vbarspinorg0(7) = i_vbarspinor(5)
            i_vbarspinorg0(8) = i_vbarspinor(6)


            o_ubarspinorg0(5) = o_ubarspinor(7)
            o_ubarspinorg0(6) = o_ubarspinor(8)
            o_ubarspinorg0(7) = o_ubarspinor(5)
            o_ubarspinorg0(8) = o_ubarspinor(6)
            
            o_vbarspinorg0(5) = o_vbarspinor(7)
            o_vbarspinorg0(6) = o_vbarspinor(8)
            o_vbarspinorg0(7) = o_vbarspinor(5)
            o_vbarspinorg0(8) = o_vbarspinor(6)

c           Compute contribution to pol sum 
            do i=0,3
              do j=0,3
                IF (use_original_routines) then
                  i_polsumU(i,j) = i_polsumU(i,j) + i_uspinor(5+i)*DCONJG(i_ubarspinor(5+j))
                  i_polsumV(i,j) = i_polsumV(i,j) + i_vspinor(5+i)*DCONJG(i_vbarspinor(5+j))
                  o_polsumU(i,j) = o_polsumU(i,j) + o_uspinor(5+i)*DCONJG(o_ubarspinor(5+j))
                  o_polsumV(i,j) = o_polsumV(i,j) + o_vspinor(5+i)*DCONJG(o_vbarspinor(5+j))
                ELSE
                  i_polsumU(i,j) = i_polsumU(i,j) + i_uspinor(5+i)*i_ubarspinor(5+j)
                  i_polsumV(i,j) = i_polsumV(i,j) + i_vspinor(5+i)*i_vbarspinor(5+j) 
                  o_polsumU(i,j) = o_polsumU(i,j) + o_uspinor(5+i)*o_ubarspinor(5+j)
                  o_polsumV(i,j) = o_polsumV(i,j) + o_vspinor(5+i)*o_vbarspinor(5+j)  
                ENDIF
              enddo
            enddo
          enddo

          else if (use_trick_routines) then

          do ihel=0,3
          
c             Particle
              call PROP_ixxxxx(p, mass, ihel, +1, left, i_uspinor)
              call PROP_ixxxxx(p, mass, ihel, +1, right, i_ubarspinor)
              call PROP_oxxxxx(p, mass, ihel, +1, left, o_uspinor)
              call PROP_oxxxxx(p, mass, ihel, +1, right, o_ubarspinor)
c             AntiParticle                                                       
              call PROP_ixxxxx(p, mass, ihel, -1, left, i_vspinor)
              call PROP_ixxxxx(p, mass, ihel, -1, right, i_vbarspinor)
              call PROP_oxxxxx(p, mass, ihel, -1, left, o_vspinor)
              call PROP_oxxxxx(p, mass, ihel, -1, right, o_vbarspinor)

c           Compute contribution to pol sum 
            do i=0,3
              do j=0,3
                i_polsumU(i,j) = i_polsumU(i,j) + i_uspinor(5+i)*i_ubarspinor(5+j)
                i_polsumV(i,j) = i_polsumV(i,j) + i_vspinor(5+i)*i_vbarspinor(5+j) 
                o_polsumU(i,j) = o_polsumU(i,j) + o_uspinor(5+i)*o_ubarspinor(5+j)
                o_polsumV(i,j) = o_polsumV(i,j) + o_vspinor(5+i)*o_vbarspinor(5+j)  
              enddo
            enddo

          enddo

          endif
          
c         Compute target polarisation sum tensor now (pslash +/-m).gamma^0
          
          pE = p(0)
          px = p(1)
          py = p(2)
          pz = p(3)

          i_target_polsumU(0,0) = pE-pz 
          i_target_polsumU(0,1) = -px + ci*py
          i_target_polsumU(0,2) = czero
          i_target_polsumU(0,3) = czero 

          i_target_polsumU(1,0) = -px -ci*py
          i_target_polsumU(1,1) = pE + pz
          i_target_polsumU(1,2) = czero
          i_target_polsumU(1,3) = czero

          i_target_polsumU(2,0) = czero
          i_target_polsumU(2,1) = czero
          i_target_polsumU(2,2) = pE+pz
          i_target_polsumU(2,3) = px -ci*py

          i_target_polsumU(3,0) = czero
          i_target_polsumU(3,1) = czero
          i_target_polsumU(3,2) = px+ci*py
          i_target_polsumU(3,3) = pE-pz

          do i=0,3
            do j=0,3
              i_target_polsumV(i,j) = i_target_polsumU(i,j)
            enddo
          enddo

          i_target_polsumU(0,2) = i_target_polsumU(0,2) + DCMPLX(mass)
          i_target_polsumU(1,3) = i_target_polsumU(1,3) + DCMPLX(mass)
          i_target_polsumU(2,0) = i_target_polsumU(2,0) + DCMPLX(mass)
          i_target_polsumU(3,1) = i_target_polsumU(3,1) + DCMPLX(mass)

          i_target_polsumV(0,2) = i_target_polsumV(0,2) - DCMPLX(mass)
          i_target_polsumV(1,3) = i_target_polsumV(1,3) - DCMPLX(mass)
          i_target_polsumV(2,0) = i_target_polsumV(2,0) - DCMPLX(mass)
          i_target_polsumV(3,1) = i_target_polsumV(3,1) - DCMPLX(mass)


c         Compute target polarisation sum tensor now DCONJG(gamma^0.(pslash +/-m))
          
          o_target_polsumU(0,0) = pE+pz 
          o_target_polsumU(0,1) = px + ci*py
          o_target_polsumU(0,2) = czero
          o_target_polsumU(0,3) = czero 

          o_target_polsumU(1,0) = px - ci*py
          o_target_polsumU(1,1) = pE - pz
          o_target_polsumU(1,2) = czero
          o_target_polsumU(1,3) = czero

          o_target_polsumU(2,0) = czero
          o_target_polsumU(2,1) = czero
          o_target_polsumU(2,2) = pE-pz
          o_target_polsumU(2,3) = -px - ci*py

          o_target_polsumU(3,0) = czero
          o_target_polsumU(3,1) = czero
          o_target_polsumU(3,2) = -px + ci*py
          o_target_polsumU(3,3) = pE+pz

          do i=0,3
            do j=0,3
              o_target_polsumV(i,j) = o_target_polsumU(i,j)
            enddo
          enddo

          o_target_polsumU(0,2) = o_target_polsumU(0,2) + DCMPLX(mass)
          o_target_polsumU(1,3) = o_target_polsumU(1,3) + DCMPLX(mass)
          o_target_polsumU(2,0) = o_target_polsumU(2,0) + DCMPLX(mass)
          o_target_polsumU(3,1) = o_target_polsumU(3,1) + DCMPLX(mass)

          o_target_polsumV(0,2) = o_target_polsumV(0,2) - DCMPLX(mass)
          o_target_polsumV(1,3) = o_target_polsumV(1,3) - DCMPLX(mass)
          o_target_polsumV(2,0) = o_target_polsumV(2,0) - DCMPLX(mass)
          o_target_polsumV(3,1) = o_target_polsumV(3,1) - DCMPLX(mass)

          tolerance = 1.0d-10

          test_passed = .True.

C         Finally compare
          IF (verbose) WRITE(*,*) ''
          IF (verbose) WRITE(*,*) '---'
          IF (verbose) WRITE(*,*) 'I:Testing \sum_h u_h(p) \bar{u}_h(p)'
          IF (verbose) WRITE(*,*) 'Tensor computed:'
          DO I=0,3
            IF (verbose) Write(*,*) 'Row ',I,' = ',(i_polsumU(I,J),J=0,3)
          ENDDO
          IF (verbose) WRITE(*,*) 'Target tensor:'
          DO I=0,3
            IF (verbose) Write(*,*) 'Row ',I,' = ',(i_target_polsumU(I,J),J=0,3)
          ENDDO
          call compare(i_polsumU, i_target_polsumU, verbose, tolerance, test_passed)
          IF (verbose) WRITE(*,*) 'O:Testing \sum_h u_h(p) \bar{u}_h(p)'
          IF (verbose) WRITE(*,*) 'Tensor computed:'
          DO I=0,3
            IF (verbose) Write(*,*) 'Row ',I,' = ',(o_polsumU(I,J),J=0,3)
          ENDDO
          IF (verbose) WRITE(*,*) 'Target tensor:'
          DO I=0,3
            IF (verbose) Write(*,*) 'Row ',I,' = ',(o_target_polsumU(I,J),J=0,3)
          ENDDO
          call compare(o_polsumU, o_target_polsumU, verbose, tolerance, test_passed)

          IF (verbose) WRITE(*,*) '---'

          IF (verbose) WRITE(*,*) 'I:Testing \sum_h v_h(p) \bar{v}_h(p)'
          IF (verbose) WRITE(*,*) 'Tensor computed:'
          DO I=0,3
            IF (verbose) Write(*,*) 'Row ',I,' = ',(i_polsumV(I,J),J=0,3)
          ENDDO
          IF (verbose) WRITE(*,*) 'Target tensor:'
          DO I=0,3
            IF (verbose) Write(*,*) 'Row ',I,' = ',(i_target_polsumV(I,J),J=0,3)
          ENDDO
          call compare(i_polsumV, i_target_polsumV, verbose, tolerance, test_passed)
          IF (verbose) WRITE(*,*) 'O:Testing \sum_h v_h(p) \bar{v}_h(p)'
          IF (verbose) WRITE(*,*) 'Tensor computed:'
          DO I=0,3
            IF (verbose) Write(*,*) 'Row ',I,' = ',(o_polsumV(I,J),J=0,3)
          ENDDO
          IF (verbose) WRITE(*,*) 'Target tensor:'
          DO I=0,3
            IF (verbose) Write(*,*) 'Row ',I,' = ',(o_target_polsumV(I,J),J=0,3)
          ENDDO
          call compare(o_polsumV, o_target_polsumV, verbose, tolerance, test_passed)


          IF (verbose) WRITE(*,*) ''

      end subroutine testCase

      subroutine compare(tensA, tensB, verbose, tolerance, test_passed)

          double complex tensA(0:3,0:3)
          double complex tensB(0:3,0:3)

          double precision tolerance
          logical did_printout, test_passed, verbose
          integer i,j

          do i=0,3
            do j=0,3
              did_printout = .False.
              if (dble(tensB(I,J)).eq.0.0d0) then
                  if (abs(dble(tensA(I,J))).gt.tolerance) then
                      test_passed = .False.
                      did_printout = .True.
                      IF (verbose) 
     & WRITE(*,*) I,J,' : ',tensA(I,J),' =/= ',tensB(I,J)
                  endif
              else
                  if (abs(dble(tensA(I,J))
     &  /dble(tensB(I,J))-1.0d0).gt.tolerance) then
                      test_passed = .False.
                      did_printout = .True.
                      IF (verbose) 
     & WRITE(*,*) I,J,' : ',tensA(I,J),' =/= ',tensB(I,J)
                  endif
              endif
              if (dimag(tensB(I,J)).eq.0.0d0) then
                  if (abs(dimag(tensA(I,J))).gt.tolerance) then
                      test_passed = .False.
                      if (.not.did_printout) then
                         IF (verbose) 
     & WRITE(*,*) I,J,' : ',tensA(I,J),' =/= ',tensB(I,J)
                      ENDIF
                  endif
              else
                  if (abs(dimag(tensA(I,J))
     &  /dimag(tensB(I,J))-1.0d0).gt.tolerance) then
                      test_passed = .False.
                      if (.not.did_printout) then                      
                         IF (verbose) 
     & WRITE(*,*) I,J,' : ',tensA(I,J),' =/= ',tensB(I,J)
                      ENDIF                         
                  endif
              endif       
            enddo
          enddo

      end subroutine compare

