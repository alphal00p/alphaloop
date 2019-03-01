      include 'max_length.inc'

      double complex ONE, CI, ZERO 
      parameter (ZERO=DCMPLX(0.0d0,0.0d0), &
		         ONE=DCMPLX(1.0d0,0.0d0), &
                 CI=DCMPLX(0.0d0,1.0d0))

      integer n_vectors
      double complex vbar(0:3)
      double complex u(0:3)

!     Index zero will hold the length of the gamma chain
      integer indices(0:max_length)
      double complex vectors(max_length,0:3)

      integer repeated_indices(0:max_length)          
      integer repeated_indices_values(-max_length:-1)

      common/chainglobals/n_vectors,vbar,u, &
       indices,vectors,repeated_indices, &
       repeated_indices_values
