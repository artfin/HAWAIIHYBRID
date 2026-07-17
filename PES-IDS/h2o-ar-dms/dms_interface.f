      subroutine RadauToCartesian(Rd,cart,mA,mB,mC)
      implicit none

      real*8,parameter:: PI=dacos(-1.d0)
      real*8,intent(in):: Rd(6),mA,mB,mC        ! in a.u. and degree !
      real*8,intent(out):: cart(4,3)    ! cart(i,j) is ith atom jth axis !
      real*8 RR,r1,r2,th1,th2,phi
      real*8 mAC,mABC
      real*8 absOB,sth,cth,th
      real*8,dimension(3):: OO,OY,COM,AC,AE,OE
      real*8,dimension(3):: OApp,OAp,OA,OCpp,OCp,OC,OBpp,OB
      data OO /0.d0, 0.d0, 0.d0/
      data OY /0.d0, 1.d0, 0.d0/

      RR=Rd(1); r1=Rd(2); r2=Rd(3)
      th1=Rd(4)*PI/180.d0
      th2=Rd(5)*PI/180.d0
      phi=Rd(6)*PI/180.d0
      mAC=mA+mC; mABC=mA+mB+mC

      OApp(1)=0.d0; OApp(2)=0.d0; OApp(3)=r1
      OCpp(1)=r2*dsin(th1); OCpp(2)=0.d0; OCpp(3)=r2*dcos(th1)
      AC(:)=OCpp(:)-OApp(:)
      AE(:)=AC(:)*mC/mAC
      OE(:)=OApp(:)+AE(:)
      OBpp(:)=OE(:)*(1.d0-dsqrt(mABC/mB))
      ! Translate the mass center of ABC to the origin !
      COM(:)=(mAC*OE(:)+mB*OBpp(:))/mABC
      OApp(:)=OApp(:)-COM(:)
      OBpp(:)=OBpp(:)-COM(:)
      OCpp(:)=OCpp(:)-COM(:)

      ! Rotate ABC in xOz plane !
      call vector_length(OBpp,absOB)
      cth=OBpp(3)/absOB
      sth=OBpp(1)/absOB
      th=dacos(cth)
      if(sth.lt.0.d0) th=2.d0*PI-th
      th=2.d0*PI-th+th2
      call rotate_point(OApp,OY,OO,th,OAp)
      call rotate_point(OBpp,OY,OO,th,OB)
      call rotate_point(OCpp,OY,OO,th,OCp)

      ! Rotate AC by phi along OB !
      call rotate_point(OAp,OB,OO,phi,OA)
      call rotate_point(OCp,OB,OO,phi,OC)

      cart(1,:)=OA(:)
      cart(2,:)=OB(:)
      cart(3,:)=OC(:)
      cart(4,1)=0.d0; cart(4,2)=0.d0; cart(4,3)=RR

      return
      endsubroutine



      subroutine rotate_point(p,vec,center,theta,np)
      implicit none

      real*8,intent(in):: p(3),vec(3),center(3),theta   ! in rad !
      real*8,intent(out):: np(3)
      real*8 s,c,c1,nx,ny,nz
      real*8 xx,xy,xz,yy,yz,zz,sx,sy,sz
      real*8 r,trans(3,3),pm(1,3),npm(1,3)

      s=dsin(theta)
      c=dcos(theta)
      c1=1.d0-c
      pm(1,:)=p(:)-center(:)

      call vector_length(vec,r)
      nx=vec(1)/r
      ny=vec(2)/r
      nz=vec(3)/r
      xx=nx*nx; xy=nx*ny; xz=nx*nz
      yy=ny*ny; yz=ny*nz; zz=nz*nz
      sx=nx*s; sy=ny*s; sz=nz*s

      trans(1,1)=xx*c1+c
      trans(2,1)=xy*c1-sz
      trans(3,1)=xz*c1+sy

      trans(1,2)=xy*c1+sz
      trans(2,2)=yy*c1+c
      trans(3,2)=yz*c1-sx

      trans(1,3)=xz*c1-sy
      trans(2,3)=yz*c1+sx
      trans(3,3)=zz*c1+c

      npm=matmul(pm,trans)
      np(:)=npm(1,:)+center(:)

      return
      endsubroutine


      subroutine normal_vector(a,b,ab)
      implicit none

      real*8,intent(in),dimension(3):: a,b
      real*8,intent(out):: ab(3)
      real*8 l,m,n,o,p,q

      l=a(1); m=a(2); n=a(3)
      o=b(1); p=b(2); q=b(3)
      ab(1)=m*q-n*p
      ab(2)=n*o-l*q
      ab(3)=l*p-m*o

      return
      endsubroutine



      subroutine cosine_theom(a,b,cth)
      implicit none

      real*8,intent(in),dimension(3):: a,b
      real*8,intent(out):: cth
      real*8 ax,ay,az,bx,by,bz,rhoa,rhob,c

      ax=a(1); ay=a(2); az=a(3)
      bx=b(1); by=b(2); bz=b(3)
      c=ax*bx+ay*by+az*bz
      rhoa=dsqrt(ax*ax+ay*ay+az*az)
      rhob=dsqrt(bx*bx+by*by+bz*bz)
      cth=c/(rhoa*rhob)
      if(isNaN(cth)) cth=0.d0
      if(cth.lt.-1.d0) cth=-1.d0
      if(cth.gt. 1.d0) cth= 1.d0

      return
      endsubroutine



      subroutine vector_length(vec,rho)
      implicit none

      real*8,intent(in):: vec(3)
      real*8,intent(out):: rho
      real*8 x,y,z

      x=vec(1); y=vec(2); z=vec(3)
      rho=dsqrt(x*x+y*y+z*z)

      return
      endsubroutine



