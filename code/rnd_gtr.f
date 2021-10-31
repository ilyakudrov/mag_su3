C**********************************************************************
C**********************************************************************
      subroutine rnd_gtr1
c not vectorized !!!
C this subrotine generates and applies random gauge transf.
C for 1st SU(2) subgroup 
c input:                                       
c u(nsite,3,3,nd) - usual lattice nonabel. config. 
c output: u(nsite,3,3,nd) after gauge transf.
C**********************************************************************
C
      include 'paravp3'

      common /var3/ u(nsite,3,3,nd)
      common /sdb2/ im(nsite,nd,n0),nm(nsite,nd)
      common /wrk2/ wg(nsite,3,3),z(nsite,3,3,nd)
c     dimension wg(nsite,3,3)
c     dimension z(nsite,3,3,nd)

      pi=4.*atan(1.d0)
      omegamax=pi/2.
      omegarange = 2. * omegamax

      do m=1,nsite
       omega1=RND(0.D0)
       omega2=RND(0.D0)
       omega3=RND(0.D0)

       omega1 = omegarange*omega1 - omegamax
       omega2 = omegarange*omega2 - omegamax
       omega3 = omegarange*omega3 - omegamax
       tmp1   = sqrt (omega1*omega1+omega2*omega2+omega3*omega3)
       if (tmp1.ne.0) then
        tmp2 = sin(tmp1)/tmp1
       else
        tmp2 = 0.
       endif
       g0= cos(tmp1)
       g1= omega1*tmp2
       g2= omega2*tmp2
       g3= omega3*tmp2

       wg(m,1,1)=dcmplx(g0,g3) 
       wg(m,1,2)=dcmplx(g2,g1) 
       wg(m,2,1)=dcmplx(-g2,g1) 
       wg(m,2,2)=dcmplx(g0,-g3) 
      enddo
       do mu=1,nd
       do i=1,3
       do j=1,3
       do m=1,nsite
         z(m,i,j,mu)=u(m,i,j,mu)
       enddo
       enddo
       enddo
       enddo
c multipl. from the left: 
       do mu=1,nd
       do i=1,2
       do j=1,3
       do m=1,nsite
        u(m,i,j,mu)=wg(m,i,1)*z(m,1,j,mu)+wg(m,i,2)*z(m,2,j,mu)
       enddo
       enddo
       enddo
       enddo
       do mu=1,nd
       do i=1,3
       do j=1,3
       do m=1,nsite
         z(m,i,j,mu)=u(m,i,j,mu)
       enddo
       enddo
       enddo
       enddo
c multipl. from the right:
       do mu=1,nd
       do i=1,2
       do j=1,3
       do m=1,nsite
        m0=nm(m,mu)
        u(m0,j,i,mu)=z(m0,j,1,mu)*dconjg(wg(m,i,1))+
     *               z(m0,j,2,mu)*dconjg(wg(m,i,2))
       enddo
       enddo
       enddo
       enddo

      RETURN
      END
C**********************************************************************
C**********************************************************************
      subroutine rnd_gtr2
c not vectorized !!!
C this subrotine generates and applies random gauge transf.
C for 2nd SU(2) subgroup 
c input:                                       
c u(nsite,3,3,nd) - usual lattice nonabel. config. 
c output: u(nsite,3,3,nd) after gauge transf.
C**********************************************************************
C
      include 'paravp3'

      common /var3/ u(nsite,3,3,nd)
      common /sdb2/ im(nsite,nd,n0),nm(nsite,nd)
      dimension wg(nsite,3,3)
      dimension z(nsite,3,3,nd)

      pi=4.*atan(1.d0)
      omegamax=pi/2.
      omegarange = 2. * omegamax

      do m=1,nsite
       omega1=RND(0.D0)
       omega2=RND(0.D0)
       omega3=RND(0.D0)

       omega1 = omegarange*omega1 - omegamax
       omega2 = omegarange*omega2 - omegamax
       omega3 = omegarange*omega3 - omegamax
       tmp1   = sqrt (omega1*omega1+omega2*omega2+omega3*omega3)
       if (tmp1.ne.0) then
        tmp2 = sin(tmp1)/tmp1
       else
        tmp2 = 0.
       endif
       g0= cos(tmp1)
       g1= omega1*tmp2
       g2= omega2*tmp2
       g3= omega3*tmp2

       wg(m,1,1)=dcmplx(g0,g3) 
       wg(m,1,3)=dcmplx(g2,g1) 
       wg(m,3,1)=dcmplx(-g2,g1) 
       wg(m,3,3)=dcmplx(g0,-g3) 
      enddo
       do mu=1,nd
       do i=1,3
       do j=1,3
       do m=1,nsite
         z(m,i,j,mu)=u(m,i,j,mu)
       enddo
       enddo
       enddo
       enddo
c multipl. from the left: 
       do mu=1,nd
       do i=1,3,2
       do j=1,3
       do m=1,nsite
        u(m,i,j,mu)=wg(m,i,1)*z(m,1,j,mu)+wg(m,i,3)*z(m,3,j,mu)
       enddo
       enddo
       enddo
       enddo
       do mu=1,nd
       do i=1,3
       do j=1,3
       do m=1,nsite
         z(m,i,j,mu)=u(m,i,j,mu)
       enddo
       enddo
       enddo
       enddo
c multipl. from the right:
       do mu=1,nd
       do i=1,3,2
       do j=1,3
       do m=1,nsite
        m0=nm(m,mu)
        u(m0,j,i,mu)=z(m0,j,1,mu)*dconjg(wg(m,i,1))+
     *               z(m0,j,3,mu)*dconjg(wg(m,i,3))
       enddo
       enddo
       enddo
       enddo

      RETURN
      END
C**********************************************************************
C**********************************************************************
      subroutine rnd_gtr3
c not vectorized !!!
C this subrotine generates and applies random gauge transf.
C for 3rd SU(2) subgroup 
c input:                                       
c u(nsite,3,3,nd) - usual lattice nonabel. config. 
c output: u(nsite,3,3,nd) after gauge transf.
C**********************************************************************
C
      include 'paravp3'

      common /var3/ u(nsite,3,3,nd)
      common /sdb2/ im(nsite,nd,n0),nm(nsite,nd)
      dimension wg(nsite,3,3)
      dimension z(nsite,3,3,nd)

      pi=4.*atan(1.d0)
      omegamax=pi/2.
      omegarange = 2. * omegamax

      do m=1,nsite
       omega1=RND(0.D0)
       omega2=RND(0.D0)
       omega3=RND(0.D0)

       omega1 = omegarange*omega1 - omegamax
       omega2 = omegarange*omega2 - omegamax
       omega3 = omegarange*omega3 - omegamax
       tmp1   = sqrt (omega1*omega1+omega2*omega2+omega3*omega3)
       if (tmp1.ne.0) then
        tmp2 = sin(tmp1)/tmp1
       else
        tmp2 = 0.
       endif
       g0= cos(tmp1)
       g1= omega1*tmp2
       g2= omega2*tmp2
       g3= omega3*tmp2

       wg(m,2,2)=dcmplx(g0,g3) 
       wg(m,2,3)=dcmplx(g2,g1) 
       wg(m,3,2)=dcmplx(-g2,g1) 
       wg(m,3,3)=dcmplx(g0,-g3) 
      enddo
       do mu=1,nd
       do i=1,3
       do j=1,3
       do m=1,nsite
         z(m,i,j,mu)=u(m,i,j,mu)
       enddo
       enddo
       enddo
       enddo
c multipl. from the left: 
       do mu=1,nd
       do i=2,3
       do j=1,3
       do m=1,nsite
        u(m,i,j,mu)=wg(m,i,2)*z(m,2,j,mu)+wg(m,i,3)*z(m,3,j,mu)
       enddo
       enddo
       enddo
       enddo
       do mu=1,nd
       do i=1,3
       do j=1,3
       do m=1,nsite
         z(m,i,j,mu)=u(m,i,j,mu)
       enddo
       enddo
       enddo
       enddo
c multipl. from the right:
       do mu=1,nd
       do i=2,3
       do j=1,3
       do m=1,nsite
        m0=nm(m,mu)
        u(m0,j,i,mu)=z(m0,j,2,mu)*dconjg(wg(m,i,2))+
     *               z(m0,j,3,mu)*dconjg(wg(m,i,3))
       enddo
       enddo
       enddo
       enddo

      RETURN
      END
C**********************************************************************

