c*********************************************************************c
      subroutine monop_def
c  04.12.00
c     * for SU(3) Lattice on VP.
c     --------------------------
c references:
c
c     * A common block, /abplaq/ is abel. config. gauge invar. plaq.
c                                to be used in mag. currents definition
c     * Before this routine is called, you must once call
c         Dirct.
c
c     * Include file
c         paravp3
c
c     * Input arrays and parameters
c         /abplaq/ : the gauge-fixed configuration.
c         /sdb2/ : the list vectors which are returned by Dirct.
c
c     * Output arrays and variables
c
c     * by V.Bornyakov 04.12.00
c----------------------------------------------------------------------

      include'paravp3'
      parameter(npl=nd*(nd-1)/2)
      common /sdb2 /im(nsite,nd,n0), nm(nsite,nd)
c     common /abplaq/ abpl(nsite,3,npl)    ! abel. gauge invar. plaq.
      common /wrk2/ abpl(nsite,3,npl)
      common /mcurr/ cmag(nsite,nd,3)      ! mag. currents

      dimension     flx(nsite,2,3)         ! ab. link 
      dimension isgn(12),nr(2)
      dimension itab(6,2)
      dimension nmon(4,3)                  ! number of monop.
      data isgn/+1,-1,-1,+1,+1,-1,+1,-1,-1,+1,+1,-1/
      data itab/3,2,2,1,1,1,4,4,3,4,3,2/


      twopi=8.*datan(1.d0)

      do i=1,3
      do mu=1,nd
      do m=1,nsite
         cmag(m,mu,i)=0.d0
      enddo
      enddo
      enddo

      do i=1,3
        np=0
	icp=0
      do mu=1,nd-1
      do nu=mu+1,nd
        np=np+1
	l=0
       do kp=1,2
	 k=itab(np,kp)
         l=l+1
	 icp=icp+1
       do m=1,nsite
         flx(m,l,i)=isgn(icp)*abpl(m,i,np)
         nr(l)=k 
         cmag(m,k,i)=cmag(m,k,i)+flx(m,l,i)
       enddo
       enddo
	 nr1=nr(1)
	 nr2=nr(2)
       do m=1,nsite
	 m1=nm(m,nr1)
         cmag(m1,nr2,i)=cmag(m1,nr2,i)-flx(m,2,i)
	 m2=nm(m,nr2)
         cmag(m2,nr1,i)=cmag(m2,nr1,i)-flx(m,1,i)
       enddo
      enddo
      enddo
      enddo

      do i=1,3
      do m=1,nsite
      do mu=1,nd
        cmag(m,mu,i)=nint(cmag(m,mu,i)/twopi)
c       cmag(m,mu,i)=cmag(m,mu,i)/twopi
      enddo
      enddo
      enddo
      do i=1,3
      do mu=1,nd
	nmon(mu,i)=0
      do m=1,nsite
       if(cmag(m,mu,i).ne.0.) then
c       write(8,*) m,mu,i,cmag(m,mu,i)
	nmon(mu,i)=nmon(mu,i)+abs(cmag(m,mu,i))
       endif
      enddo
      enddo
      enddo
      do i=1,3
       write(8,*) (nmon(mu,i),mu=1,nd)
      enddo
      do mu=1,nd
      do m=1,nsite
         nf=cmag(m,mu,1)+cmag(m,mu,2)+cmag(m,mu,3)
        if(nf.ne.0) write(8,*) 'error', m,mu,(cmag(m,mu,i),i=1,3)
      enddo
      enddo
      return
      end
