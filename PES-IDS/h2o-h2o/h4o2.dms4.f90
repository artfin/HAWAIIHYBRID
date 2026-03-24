        subroutine predip()

        implicit double precision (a-h,o-z)
        implicit integer (i-n)

        double precision V, cart_in(6,3), v0
        double precision dc0(0:5,0:5),dw0(0:5,0:5), coef(0:2879)

!        integer,parameter :: m=1976, mr=15
!        double precision coef(0:m+3*mr-1)
!        double precision dc0,dc1,dc2,dw0,dw1,dw2

        integer i,j,k,i1,j1,k1

        common/NCOEdip/ms,mr
        common/h4o2coefdip/dc0,dw0,coef

        ms=1827 ; mr=1053

        open(20,file='h4o2.dms4.coeff.dat',status='old')
        read(20,*)
        read(20,*)
        read(20,*)(coef(i1),i1=0,ms+mr-1)
        close(20)

        return
        end

!***************************************************************
        subroutine calcdip(dip,cart_in)

        implicit none

        integer ms,mr
        double precision dc0(0:5,0:5),dw0(0:5,0:5),coef(0:2879)

        common/NCOEdip/ms,mr
        common/h4o2coefdip/dc0,dw0,coef
        double precision V, cart_in(6,3), dip(3)

        double precision cart0(6,3),cart1(6,3)
        integer i,j,k,l,i1,j1,k1,l1,i2,j2,k2,l2


        double precision rvec(0:3),d0(0:5,0:5),r0(0:5,0:5)
        double precision xnuc(0:2,0:5),vec(0:3,0:2879)

         do j=1,3
          xnuc(j-1,0:1)=cart_in(2:3,j)
          xnuc(j-1,2:3)=cart_in(5:6,j)
          xnuc(j-1,4)=cart_in(1,j)
          xnuc(j-1,5)=cart_in(4,j)
        end do

        call getdvec (ms, mr, xnuc(0:2,0:5), vec)
!        dip(1:3) = matmul(vec(0:2,0:2879),coef)

        do i=1,3
          dip(i)=dot_product(vec(i-1,:),coef)
        end do

        return
        end subroutine calcdip


