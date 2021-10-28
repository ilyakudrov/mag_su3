C**********************************************************************
C**********************************************************************
      subroutine sweep_or(am,fi)
c not vectorized !!!!
C updating for adjoint spin model with local couplings
C this subrotine uses overrelaxationprocedure 
c input:                                       +
c am(nsite,nd,3,3) - matrix (1/2)*tr(U(x,id)*sigma_i*U(x,id)*sigma_j)
c                                       
c fi(nsite,3) - related tp gauge transformation: 
c                + 
c               g(x)*sigma_3*g(x)=i*fi(x,i)*sigma_i
c output:                                       +
c new values of fi(nsite,3)
C**********************************************************************
C
C routine deals with adjoint spin model defined by
C S=SUM(1-0.5*(TR F(x)*U(x,mue)*F(x+mue)*U(x+mue,mue))
C
C
      include 'paravp3'
c     include 'su2mon.cmn'
      common /sdb2 / im(nsite,nd,n0),nm(nsite,nd)

      dimension am(nsite,3,3,nd),fi(nsite,3)
      dimension si(3)


c     write(6,*) 'sweep_or started'
c     if(is.eq.0) then
c fi starting values:
c     do j=1,nsite
c       fi(j,1)=0.0
c       fi(j,2)=0.0
c       fi(j,3)=1.0
c     enddo
c     endif

c initial value (the lines to be removed):
c     R=0.d0
c     do i=1,3
c     do j=1,3
c     do mu=1,nd
c     do m=1,nsite
c      ns=im(m,mu,1) ! positive direct. step
c      R=R+am(m,i,j,mu)*fi(m,i)*fi(ns,j)
c     enddo
c     enddo
c     enddo
c     enddo
c     write(2,*) 'Ri=',R

       do ns1=1,nsite

c   compute the 'open plaquette'
       do i=1,3
         si(i)=0.d0
       enddo
       do id=1,nd  ! over directions
        ns2=nm(ns1,id)   ! negative direct. step
        ns3=im(ns1,id,1) ! positive direct. step
C positive and negative directions :
        do i=1,3
        do j=1,3
        si(i)=si(i)+am(ns1,i,j,id)*fi(ns3,j)+am(ns2,j,i,id)*fi(ns2,j)
        enddo
        enddo
       enddo

C*****************************************************************************
c the transformation : overrel. 
C*****************************************************************************
       dm=sqrt(si(1)*si(1)+si(2)*si(2)+si(3)*si(3))
c****************
       if(dm.eq.0.) then
         write(6,*) 'dm=',dm
         stop
       endif
c****************
       dm=1./dm
c transformation matrix:
c normalization:
      dvn=0.d0
      do i=1,3
        dvn=dvn+fi(ns1,i)*si(i)*dm
      enddo
      vn2=2.*dvn
c        it's the output of overrelax.:
      do i=1,3
        fi(ns1,i)=-fi(ns1,i)+vn2*si(i)*dm
      enddo
C*****************************************************************************
      enddo

c final value:
c     R=0.d0
c     do i=1,3
c     do j=1,3
c     do mu=1,nd
c     do m=1,nsite
c      ns=im(m,mu,1) ! positive direct. step
c      R=R+am(m,i,j,mu)*fi(m,i)*fi(ns,j)
c     enddo
c     enddo
c     enddo
c     enddo
c     write(2,*) 'Rf=',R/nlink

c     write(6,*) 'sweep_or finished'
      RETURN
      END
C**********************************************************************
