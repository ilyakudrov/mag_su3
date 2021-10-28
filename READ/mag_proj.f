c*********************************************************************c
      subroutine mag_proj
c  22.11.00
c     * for SU(3) Lattice on VP.
c     --------------------------
c references:
c M.Polikarpov, K.Yee, hep-lat/9306014
c
c     * A common block, /wrk1/ is the gauge-fixed configuration.
c
c     * Before this routine is called, you must once call
c         Dirct.
c
c     * Include file
c         paravp3
c
c     * Input arrays and parameters
c         /var3/ : the gauge-fixed configuration.
c         /tabl/ : the list vectors which are returned by Dirct.
c
c     * Output arrays and variables
c
c     * by V.Bornyakov 22.11.00
c----------------------------------------------------------------------

      include'paravp3'
      parameter(npl=nd*(nd-1)/2)
      common /var3/ z(nsite,3,3,nd)
c     common /wrk1/  x(nvect,3,3,nd),y(nvect,3,3,nd)
c     common /tb1/  mioe(nsite), mioo(nsite)
c     common /tabl /mo1(nvect,nd), mo2(nvect,nd),
c    &              me1(nvect,nd), me2(nvect,nd),
c    &              ldir2(nd-1,nd)
      common /sdb2 /im(nsite,nd,n0), nm(nsite,nd)
c     common /abplaq/ abpl(nsite,3,npl)  ! ab. plaq
      common/abli/  ab(nsite,3,nd)         ! ab. link

      common /wrk2/ abpl(nsite,3,npl), ph(nsite,nd),fabp(nsite,npl)
c     dimension     zx(nvect,3,3,nd),zy(nvect,3,3,nd),wg(nvect,3,3)
c    c             ,qdi(nvect),qod(nvect),qre(nvect),qim(nvect)
c     dimension     ab(nsite,3,nd)         ! ab. link 
c     dimension     abp(nsite,3,npl),araux(3)         ! auxil. 
c     dimension     ph(nsite,nd),fabp(nsite,npl)  ! auxil.
      dimension     av(3)


      pi=4.*datan(1.d0)
      twopi=8.*datan(1.d0)
      sum=0.
      do i=1,3
      do l=1,nd
      do m=1,nsite 
       sum=sum+z(m,i,i,l)*conjg(z(m,i,i,l))
      enddo
      enddo
      enddo
      sum=sum/(8.*nsite) - 0.5
      write(2,*)'R=',sum
c abel. link angle calculation:
      do l=1,nd
      do m=1,nsite 
       ph(m,l) = (dimag(log(z(m,1,1,l)))+dimag(log(z(m,2,2,l)))+
     &             dimag(log(z(m,3,3,l))))
c    &             dimag(log(z(m,3,3,l))))/3. 
      enddo
      enddo
      do l=1,nd
      do m=1,nsite
       ph(m,l) = ph(m,l) - twopi*nint( ph(m,l)/twopi )
       if(abs(ph(m,l)).gt.pi) write(2,*) 'error, phi'
      enddo
      enddo

      do l=1,nd
      do i=1,3
      do m=1,nsite
        ab(m,i,l) = dimag(log(z(m,i,i,l))) - ph(m,l)/3.
      enddo
      enddo
      enddo
c abel. plaq. calculation:
      np=0
      do mu=1,nd-1
      do nu=mu+1,nd
	np=np+1
      do i=1,3
      do m1=1,nsite
        m2=im(m1,mu,1)
        m4=im(m1,nu,1)
        abpl(m1,i,np) = ab(m1,i,mu)+ab(m2,i,nu)-ab(m4,i,mu)-ab(m1,i,nu)
        abpl(m1,i,np) = abpl(m1,i,np) - twopi*nint( abpl(m1,i,np)/twopi)
      enddo
      enddo
      enddo
      enddo
      do np=1,6
      do m=1,nsite
       fabp(m,np)=abpl(m,1,np) + abpl(m,2,np) + abpl(m,3,np) 
      enddo
      enddo
      do np=1,6
      do m=1,nsite
        if(abs(fabp(m,np)).gt.0.0001) then
c         do i=1,3
c          araux(i)=abpl(m,i,np)
c         enddo
         if(fabp(m,np).gt.0.0001) then
	  abmax=0.
	  do i=1,3 
	   if(abpl(m,i,np).gt.abmax) then
	    abmax=abpl(m,i,np)
	    naux=i
           endif
	  enddo   
c          naux=maxloc(araux)              ! location of max. value(F90)
           abpl(m,naux,np)=abpl(m,naux,np) - twopi
         else
	  abmin=0.
	  do i=1,3 
	   if(abpl(m,i,np).lt.abmin) then
	    abmin=abpl(m,i,np)
	    naux=i
           endif
	  enddo   
c          naux=minloc(araux)              ! location of min. value(F90)
           abpl(m,naux,np)=abpl(m,naux,np) + twopi
         endif
        endif
      enddo
      enddo
c check:
      do np=1,6
      do m=1,nsite
       fabp(m,np)=abpl(m,1,np) + abpl(m,2,np) + abpl(m,3,np) 
       if(abs(fabp(m,np)).gt.0.0001) write(7,*) m,np,fabp(m,np)
      enddo
      enddo

c average:
      do i=1,3 
       av(i) = 0.
      do np=1,6
      do m=1,nsite
       av(i) = av(i)+abpl(m,i,np)
      enddo
      enddo
       av(i) =av(i)/(6.*nsite)
      enddo
      write(6,*) 'aver. abel. plaq. =',av

 2000 return
      end
c*********************************************************************c
c*********************************************************************c
      subroutine sp_action(temp)
c
c to compute spin action: |u(1,1)|**2 + |u(2,2)|**2 +|u(3,3)|**2 
c
c     * Input arrays and parameters
c         /var3/ : the gauge-fixed configuration.
c         /tabl/ : the list vectors which are returned by Dirct.
c
c     * Output arrays and variables
c
c     * by V.Bornyakov 11.12.00
c----------------------------------------------------------------------
      include'paravp3'
      parameter(npl=nd*(nd-1)/2)
      common /var3/ z(nsite,3,3,nd)


      sum=0.
      do i=1,3
      do l=1,nd
      do m=1,nsite 
       sum=sum+z(m,i,i,l)*conjg(z(m,i,i,l))
      enddo
      enddo
      enddo
      sum=sum/(8.*nsite) - 0.5
      write(2,*)temp,sum

      return
      end
