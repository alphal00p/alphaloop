       subroutine compute_one_loop_tadpole()
            use avh_olo_dp
            write(*,*) 'Tadpoles not implemented yet in AVHOneLoopHook'
       end subroutine compute_one_loop_tadpole

       subroutine compute_one_loop_bubble()
            use avh_olo_dp
            write(*,*) 'Bubbles not implemented yet in AVHOneLoopHook'
       end subroutine compute_one_loop_bubble

       subroutine compute_one_loop_triangle(p1_in,p2_in,p3_in,
     &                     m1_in,m2_in,m3_in,res)
            use avh_olo_dp
CF2PY INTENT(OUT) :: res
CF2PY INTENT(IN)  :: p1_in,p2_in,p3_in,m1_in,m2_in,m3_in

            double complex res(0:2)
            double precision p1_in,p2_in,p3_in,m1_in,m2_in,m3_in
            double complex p1,p2,p3,m1,m2,m3

            p1 = DCMPLX(p1_in, 0.0d0)
            p2 = DCMPLX(p2_in, 0.0d0)
            p3 = DCMPLX(p3_in, 0.0d0)
            m1 = DCMPLX(m1_in, 0.0d0)
            m2 = DCMPLX(m2_in, 0.0d0)
            m3 = DCMPLX(m3_in, 0.0d0)

!           avh_olo_d0c
            call olo( res ,p1,p2,p3,m1,m2,m3 )
!            write(*,*) 'res=',res(0)
       end subroutine compute_one_loop_triangle

       subroutine compute_one_loop_box(p1_in,p2_in,p3_in,p4_in,
     &                     p12_in,p23_in,m1_in,m2_in,m3_in,m4_in,res)
            use avh_olo_dp

CF2PY INTENT(OUT) :: res
CF2PY INTENT(IN)  :: p1_in,p2_in,p3_in,p4_in,p12_in,p23_in,m1_in,m2_in,m3_in,m4_in

            double complex res(0:2)
            double precision p1_in,p2_in,p3_in,p4_in,p12_in,p23_in,
     &                                        m1_in,m2_in,m3_in,m4_in
            double complex p1,p2,p3,p4,p12,p23,m1,m2,m3,m4

            p1 = DCMPLX(p1_in, 0.0d0)
            p2 = DCMPLX(p2_in, 0.0d0)
            p3 = DCMPLX(p3_in, 0.0d0)
            p4 = DCMPLX(p4_in, 0.0d0)
            p12 = DCMPLX(p12_in, 0.0d0)
            p23 = DCMPLX(p23_in, 0.0d0)
            m1 = DCMPLX(m1_in, 0.0d0)
            m2 = DCMPLX(m2_in, 0.0d0)
            m3 = DCMPLX(m3_in, 0.0d0)
            m4 = DCMPLX(m4_in, 0.0d0)

!           avh_olo_d0c
            call olo( res ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
!            write(*,*) 'res=',res(0)
       end subroutine compute_one_loop_box
