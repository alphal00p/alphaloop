
      subroutine gamma_multiply(flow, gamma_index, factor, res_flow)
          implicit None

!         -------
!         Variables declaration
!         -------

          double complex ONE, CI 
          parameter (ONE=DCMPLX(1.0d0,0.0d0), &
                     CI=DCMPLX(0.0d0,1.0d0))

!         The flow is modified in place
          double complex flow(0:3), res_flow(0:3)
          integer gamma_index
          double complex factor

!         -------
!         Implementation
!         ------

          if (gamma_index.eq.0) then
              res_flow(0) = res_flow(0) + flow(0) * factor
              res_flow(1) = res_flow(1) + flow(1) * factor
              res_flow(2) = res_flow(2) - flow(2) * factor
              res_flow(3) = res_flow(3) - flow(3) * factor
          elseif (gamma_index.eq.1) then
              res_flow(0) = res_flow(0) - flow(3) * factor
              res_flow(1) = res_flow(1) - flow(2) * factor
              res_flow(2) = res_flow(2) + flow(1) * factor
              res_flow(3) = res_flow(3) + flow(0) * factor
          elseif (gamma_index.eq.2) then
              res_flow(0) = res_flow(0) - CI * flow(3) * factor
              res_flow(1) = res_flow(1) + CI * flow(2) * factor
              res_flow(2) = res_flow(2) + CI * flow(1) * factor
              res_flow(3) = res_flow(3) - CI * flow(0) * factor
          elseif (gamma_index.eq.3) then
              res_flow(0) = res_flow(0) - flow(2) * factor
              res_flow(1) = res_flow(1) + flow(3) * factor
              res_flow(2) = res_flow(2) + flow(0) * factor
              res_flow(3) = res_flow(3) - flow(1) * factor
          elseif (gamma_index.eq.5) then
              res_flow(0) = res_flow(0) + flow(2) * factor
              res_flow(1) = res_flow(1) + flow(3) * factor
              res_flow(2) = res_flow(2) + flow(0) * factor
              res_flow(3) = res_flow(3) + flow(1) * factor
          endif

      end subroutine gamma_multiply

      subroutine compute_chain_with_dummy_indices_fixed(res)
          implicit None

!         -------
!         Variables declaration
!         -------

          integer i,j,k

          double complex res
          double complex flow(0:3), next_flow(0:3), res_flow(0:3)
          include 'globals.inc'

!         -------
!         Implementation
!         ------

!         Follow the spinor flow
          do k=0,3
            flow(k) = vbar(k)
          enddo

          do i=1,indices(0)
            do k=0,3
               next_flow(k) = ZERO 
            enddo
            if (indices(i).gt.0) then
!              Contraction with an external vector
               call gamma_multiply(&
                   flow, 0, vectors(indices(i),0), next_flow)
               do j=1,3
                 call gamma_multiply(&
                   flow, j, -vectors(indices(i),j), next_flow)
               enddo
            else
!              Internal summation
               if(repeated_indices_values(indices(i)).eq.0) then
                 call gamma_multiply(flow, &
                   repeated_indices_values(indices(i)), ONE, next_flow) 
               else
                 call gamma_multiply(flow, &
                   repeated_indices_values(indices(i)), CI, next_flow) 
               endif
            endif
!           Update the running flow
            do k=0,3
               flow(k) = next_flow(k)
            enddo
          enddo

!         Build and return the result
          res = DCMPLX(0.0d0,0.0d0)
          do i=0,3
            res = res + flow(i)*u(i)
          enddo

      end subroutine compute_chain_with_dummy_indices_fixed

      recursive subroutine compute_chain_element( &
       in_repeated_index_pointer, res)
          implicit None

!         -------
!         Variables declaration
!         -------

          integer i,j

          include 'globals.inc'
          integer repeated_index_pointer, in_repeated_index_pointer

          double complex res, partial_res

!         -------
!         Implementation
!         -------

!         Instantiate the result
          res = DCMPLX(0.0d0,0.0d0)

!         If the repeated_index_pointer is beyond the number of
!         repeated indices then return
          repeated_index_pointer = in_repeated_index_pointer
          if (repeated_index_pointer .gt. repeated_indices(0)) then
              call compute_chain_with_dummy_indices_fixed(res)
          else
!             move the dummy index pointer
              repeated_index_pointer = repeated_index_pointer + 1
!             Then perform the summation over that dummy index
              do i=0,3
                repeated_indices_values( &
                     repeated_indices(repeated_index_pointer-1)) = i
                call compute_chain_element( &
                     repeated_index_pointer, partial_res)
                res = res +  partial_res
              enddo
          endif

      end subroutine compute_chain_element

      subroutine compute_chain( &
            in_vbar, in_u, in_indices, n_indices, in_vectors, &
            in_n_vectors, res)
          implicit None

!f2py intent(out) :: res
!f2py integer intent(hide),depend(in_vectors) :: in_n_vectors=shape(in_vectors,0)
!f2py integer intent(hide),depend(in_indices) :: n_indices=shape(in_indices,0)

!         -------
!         Variables declaration
!         -------
          integer i,j

          include 'globals.inc'

          integer in_n_vectors, n_indices
          double complex in_vbar(4)
          double complex in_u(4)

!         Index zero will hold the length of the gamma chain
          integer in_indices(n_indices)
          double complex in_vectors(in_n_vectors,4)
          
          integer repeated_index_pointer

          double complex res

!         -------
!         Implementation
!         -------
          
!         Set the globals
          do i=0,3
            vbar(i) = in_vbar(i+1)
            u(i) = in_u(i+1)
            n_vectors = in_n_vectors
            do j = 1,n_vectors
              vectors(j,i) = in_vectors(j,i+1)
            enddo
          enddo
          
          indices(0) = n_indices
          do i=1,n_indices
            indices(i) = in_indices(i)
          enddo


!         First initialise the value of all repeated indices values to
!         -1, i.e. not used
          do i=-max_length,-1
            repeated_indices_values(i) = -1
          enddo

!         Then list all reapeated indices used
          repeated_indices(0) = 0
          do i=1,indices(0)
            if ((indices(i).lt.0).and. &
               (repeated_indices_values(indices(i)).lt.0)) then
!             Add it to the list
              repeated_indices(0) = repeated_indices(0)+1
              repeated_indices(repeated_indices(0)) = indices(i)
!             Flag this repeated as found setting its initial value to 0
              repeated_indices_values(indices(i)) = 0
            endif
          enddo

!         Instantiate the result
          res = DCMPLX(0.0d0,0.0d0)

!         Place the repeated index pointer at 1 to initialise the recursion
          repeated_index_pointer = 1
          
          call compute_chain_element(repeated_index_pointer, res)

      end subroutine compute_chain

      program ChainExample
          implicit None

!         -------
!         Variables declaration
!         -------

          integer n_vectors, n_indices
          parameter(n_vectors=2, n_indices=14)

          include 'max_length.inc'
          double complex final_result

          double complex vectors(2,0:3)
          double complex vbar(0:3)
          double complex u(0:3)

          integer indices(n_indices)
 
!         -------
!         Implementation
!         -------
          
!         Specify the number of vectors and their values
          vectors(1,0) = DCMPLX(1.0d0,0.d0)
          vectors(1,1) = DCMPLX(1.1d0,0.d0)
          vectors(1,2) = DCMPLX(1.2d0,0.d0)
          vectors(1,3) = DCMPLX(1.3d0,0.d0)
          vectors(2,0) = DCMPLX(2.0d0,0.d0)
          vectors(2,1) = DCMPLX(2.1d0,0.d0)
          vectors(2,2) = DCMPLX(2.2d0,0.d0)
          vectors(2,3) = DCMPLX(2.3d0,0.d0)

!         Specify the number of indices and their values (negative
!         repeated are summed)

!         Index of vector to contract with first gamma 
          indices(1) = 1
!         Dummy index to contract gamma matrices within a chain
          indices(2) = -1
!         Second dummy index to contract gamma matrices within a chain
          indices(3) = -2
!         Index of vector to contract with third gamma
          indices(4) = 2
!         Contract this gamma matrix with the second one
          indices(5) = -1
!         Contract this gamma matrix with the third one
          indices(6) = -2
!         And for fun add to this product 8 interleaved gamma matrices
!         :)
          indices(7) = -3
          indices(8) = -4
          indices(9) = -5
          indices(10) = -6
          indices(11) = -5
          indices(12) = -4
          indices(13) = -6
          indices(14) = -3


!         Definition of the endpoint of the gamma chain
          vbar(0) = DCMPLX(3.0d0,0.d0)
          vbar(1) = DCMPLX(3.1d0,0.d0)
          vbar(2) = DCMPLX(3.2d0,0.d0)
          vbar(3) = DCMPLX(3.3d0,0.d0)
          u(0) = DCMPLX(4.0d0,0.d0)
          u(1) = DCMPLX(4.1d0,0.d0)
          u(2) = DCMPLX(4.2d0,0.d0)
          u(3) = DCMPLX(4.3d0,0.d0)


!         Main routine call
          call compute_chain( &
              vbar, u, indices, n_indices, vectors, n_vectors,& 
              final_result)

          write(*,*) 'Result from the numerical evaluation of the'//&
                     ' gamma matrix chain evaluation:'
          write(*,*) final_result
          write(*,*) 'Target result from test chain should be:'
          write(*,*) DCMPLX(-78354.8416d0,1312.2560000000d0)

      end program ChainExample
