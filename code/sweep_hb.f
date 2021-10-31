C**********************************************************************
C**********************************************************************
      subroutine sweep_hb(am,fi,temp)
c not vectorized !!!
C updating for adjoint spin model with local couplings
C this subrotine uses standard heatbath procedure
C without overrelaxation.
c input:a                                       +
c am(nsite,nd,3,3) - matrix (1/2)*tr(U(x,id)*sigma_i*U(x,id)*sigma_j)
c
c fi(nsite,3) - related tp gauge transformation:
c                +
c               g(x)*sigma_3*g(x)=i*fi(x,i)*sigma_i
C**********************************************************************
C
C routine deals with adjoint spin model defined by
C S=SUM(1-0.5*(TR F(x)*U(x,mue)*F(x+mue)*U(x+mue,mue))
C
C
      include 'paravp3'
c     include 'su2mon.cmn'

      common /var3/ u(nsite,3,3,nd)
      common /sdb2 / im(nsite,nd,n0),nm(nsite,nd)
      dimension am(nsite,3,3,nd),fi(nsite,3)
      dimension si(3)

      tem=temp
c     write(6,*) 'tem=',tem
       bet=1./tem
c      write(6,*) 'beta=',bet

c     tem=1./betsp

c     if(is.eq.0) then
c fi starting values:
c     do j=1,nsite
c      fi(j,1)=0.0
c      fi(j,2)=0.0
c      fi(j,3)=1.0
c     enddo
c     endif

c     write(6,*)'sweep_hb started'
c initial value:
c     r=0.d0
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
c     write(2,*) 'Ri=',r/nlink

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
c to choose the transformation : overrel. or heatbath
C*****************************************************************************
c      if(mod(NS1,7).ne.0) then
C*****************************************************************************
C   OVERRELAXATION :
C*****************************************************************************
c      DM=sqrt(SI(1)*SI(1)+SI(2)*SI(2)+SI(3)*SI(3))
c****************
c      if(DM.eq.0.) then
c        write(6,*) 'DM=',DM
c        stop
c      endif
c****************
c      DM=1./DM
Ctransformation matrix:
c normalization:
c     DVN=0.d0
c     do i=1,3
c       DVN=DVN+FI(NS1,i)*SI(i)*DM
c     enddo
c     VN2=2.*DVN
c     do i=1,3
c       FI(NS1,i)=-FI(NS1,i)+VN2*SI(i)*DM
c     enddo
C*****************************************************************************
c     ELSE
C*****************************************************************************
C   updating by the heat bath method:
C*****************************************************************************

       dm=sqrt(si(1)*si(1)+si(2)*si(2)+si(3)*si(3))
c****************
       if(dm.eq.0.) then
         write(6,*) 'dm=',dm
         stop
       endif
c****************
       d=dm*bet
       dm=1./dm
       ed=exp(d)
       emd=1./ed
       de=ed-emd
 6     r3=rnd(0.d0)*de+emd
       r3=dlog(r3)/d
       r=1.-r3*r3
       if(r.lt.0) go to 6
 7     r1=rnd(0.d0)-0.5
       r2=rnd(0.d0)-0.5
       dv=r1*r1+r2*r2
       if(dv.gt.0.25) go to 7

c       X3=X3
       r=sqrt(r/dv)
       r1=r1*r
       r2=r2*r
c  T
c o *si = (0,0,1)
       cna=si(3)*dm
c****************
       if(abs(cna).gt.0.9999999) then
         write(6,*) 'cna=',cna
         fi(ns1,1)= r1
         fi(ns1,2)= r2
         fi(ns1,3)= r3*sign(1.d0,cna)
c****************
       else
         sna=sqrt(1.-cna*cna)
         if(sna.eq.0.) write(6,*) 'sna=',sna
         sna1=1./sna
         snf=si(1)*sna1*dm
         cnf=si(2)*sna1*dm
         o11=cnf
         o12=snf*cna
         o13=snf*sna
         o21=-snf
         o22=cnf*cna
         o23=cnf*sna
c         o31=0.
         o32=-sna
         o33=cnA
C******************************************************************************
C        IT'S THE OUTPUT OF heatbath:
C******************************************************************************
          fi(ns1,1)= o11*r1+o12*r2+o13*r3
          fi(ns1,2)= o21*r1+o22*r2+o23*r3
          fi(ns1,3)=        o32*r2+o33*r3
        endif

c     tem=tem-ddtemp
c     ENDIF

      enddo

c final value:
c      R=0.d0
c      do i=1,3
c      do j=1,3
c      do mu=1,nd
c      do m=1,nsite
c       ns=im(m,mu,1) ! positive direct. step
c       R=R+am(m,i,j,mu)*fi(m,i)*fi(ns,j)
c      enddo
c      enddo
c      enddo
c      enddo
c      write(2,*) 'Rf=',R/nlink

c      write(6,*) 'temp=',tem,'  beta=',bet
c      write(6,*)'sweep_hb finished'

      RETURN
      END
C**********************************************************************
