C**********************************************************************
C**********************************************************************
      subroutine gtr1(fi)
c not vectorized !!!
C this subrotine computes gauge transformation
c for given fi
c input:                                       +
c fi(nsite,3) - related tp gauge transformation: 
c                + 
c               wg(x)*sigma_3*wg(x)=fi(x,i)*sigma_i
c output: gauge transformation  wg(x)
C**********************************************************************
C
      include 'paravp3'

      common /var3/ u(nsite,3,3,nd)
      common /sdb2/ im(nsite,nd,n0),nm(nsite,nd)
      common /wrk2/ wg(nsite,3,3),z(nsite,3,3,nd)
      dimension fi(nsite,3)
c     dimension fi(nsite,3),wg(nsite,3,3)
c     dimension z(nsite,3,3,nd)


c initial value:
      r=0.d0
      do mu=1,nd
      do m=1,nsite
       r = r + u(m,1,1,mu)*dconjg(u(m,1,1,mu)) +
     *         u(m,2,2,mu)*dconjg(u(m,2,2,mu))          
     *       - u(m,1,2,mu)*dconjg(u(m,1,2,mu))          
     *       - u(m,2,1,mu)*dconjg(u(m,2,1,mu))          
      enddo
      enddo
c     write(2,*)'Ri=',r/(2.*nlink)

       do m=1,nsite
        if(fi(m,3).le.-1.0) then
          g0=0.
          g1=0.
          g2=1.
        else
          gn=sqrt(2.0*(1.+fi(m,3)))
          g0=gn*0.5
          gn=1.0/gn
          g1=-fi(m,2)*gn
          g2= fi(m,1)*gn
        endif
        wg(m,1,1)=dcmplx(g0) 
        wg(m,1,2)=dcmplx(g2,g1) 
        wg(m,2,1)=dcmplx(-g2,g1) 
        wg(m,2,2)=dcmplx(g0) 
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

c final value:
      r=0.d0
      do mu=1,nd
      do m=1,nsite
       r = r + u(m,1,1,mu)*dconjg(u(m,1,1,mu)) +
     *         u(m,2,2,mu)*dconjg(u(m,2,2,mu))          
     *       - u(m,1,2,mu)*dconjg(u(m,1,2,mu))          
     *       - u(m,2,1,mu)*dconjg(u(m,2,1,mu))          
      enddo
      enddo
c     write(2,*)'after gtr. Rf=',r/(2.*nlink)

      RETURN
      END
C**********************************************************************
C**********************************************************************
      subroutine gtr2(fi)
c not vectorized !!!
C this subrotine computes gauge transformation
c for given fi
c input:                                       +
c fi(nsite,3) - related tp gauge transformation: 
c                + 
c               wg(x)*sigma_3*wg(x)=fi(x,i)*sigma_i
c output: gauge transformation  wg(x)
C**********************************************************************
C
      include 'paravp3'

      common /var3/ u(nsite,3,3,nd)
      common /sdb2/ im(nsite,nd,n0),nm(nsite,nd)
      dimension fi(nsite,3),wg(nsite,3,3)
      dimension z(nsite,3,3,nd)


c initial value:
      r=0.d0
      do mu=1,nd
      do m=1,nsite
       r = r + u(m,1,1,mu)*dconjg(u(m,1,1,mu)) +
     *         u(m,3,3,mu)*dconjg(u(m,3,3,mu))          
     *       - u(m,1,3,mu)*dconjg(u(m,1,3,mu))          
     *       - u(m,3,1,mu)*dconjg(u(m,3,1,mu))          
      enddo
      enddo
c     write(2,*)'Ri=',r/(2.*nlink)

       do m=1,nsite
        if(fi(m,3).le.-1.0) then
          g0=0.
          g1=0.
          g2=1.
        else
          gn=sqrt(2.0*(1.+fi(m,3)))
          g0=gn*0.5
          gn=1.0/gn
          g1=-fi(m,2)*gn
          g2= fi(m,1)*gn
        endif
        wg(m,1,1)=dcmplx(g0) 
        wg(m,1,3)=dcmplx(g2,g1) 
        wg(m,3,1)=dcmplx(-g2,g1) 
        wg(m,3,3)=dcmplx(g0) 
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

c final value:
      r=0.d0
      do mu=1,nd
      do m=1,nsite
       r = r + u(m,1,1,mu)*dconjg(u(m,1,1,mu)) +
     *         u(m,3,3,mu)*dconjg(u(m,3,3,mu))          
     *       - u(m,1,3,mu)*dconjg(u(m,1,3,mu))          
     *       - u(m,3,1,mu)*dconjg(u(m,3,1,mu))          
      enddo
      enddo
c     write(2,*)'after gtr. Rf=',r/(2.*nlink)

      RETURN
      END
C**********************************************************************
C**********************************************************************
      subroutine gtr3(fi)
c not vectorized !!!
C this subrotine computes gauge transformation
c for given fi
c input:                                       +
c fi(nsite,3) - related tp gauge transformation: 
c                + 
c               wg(x)*sigma_3*wg(x)=fi(x,i)*sigma_i
c output: gauge transformation  wg(x)
C**********************************************************************
C
      include 'paravp3'

      common /var3/ u(nsite,3,3,nd)
      common /sdb2/ im(nsite,nd,n0),nm(nsite,nd)
      dimension fi(nsite,3),wg(nsite,3,3)
      dimension z(nsite,3,3,nd)


c initial value:
      r=0.d0
      do mu=1,nd
      do m=1,nsite
       r = r + u(m,2,2,mu)*dconjg(u(m,2,2,mu)) +
     *         u(m,3,3,mu)*dconjg(u(m,3,3,mu))          
     *       - u(m,2,3,mu)*dconjg(u(m,2,3,mu))          
     *       - u(m,3,2,mu)*dconjg(u(m,3,2,mu))          
      enddo
      enddo
c     write(2,*)'Ri=',r/(2.*nlink)

       do m=1,nsite
        if(fi(m,3).le.-1.0) then
          g0=0.
          g1=0.
          g2=1.
        else
          gn=sqrt(2.0*(1.+fi(m,3)))
          g0=gn*0.5
          gn=1.0/gn
          g1=-fi(m,2)*gn
          g2= fi(m,1)*gn
        endif
        wg(m,2,2)=dcmplx(g0) 
        wg(m,2,3)=dcmplx(g2,g1) 
        wg(m,3,2)=dcmplx(-g2,g1) 
        wg(m,3,3)=dcmplx(g0) 
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

c final value:
      r=0.d0
      do mu=1,nd
      do m=1,nsite
       r = r + u(m,2,2,mu)*dconjg(u(m,2,2,mu)) +
     *         u(m,3,3,mu)*dconjg(u(m,3,3,mu))          
     *       - u(m,2,3,mu)*dconjg(u(m,2,3,mu))          
     *       - u(m,3,2,mu)*dconjg(u(m,3,2,mu))          
      enddo
      enddo
c     write(2,*)'after gtr. Rf=',r/(2.*nlink)

      return
      end
C**********************************************************************

