c*********************************************************************c
      program vpsu3
c---------------------------------------------------------------------c
c     * This is a main routine to perform                             c
c       monte-carlo simulations of SU(3) lattice gauge theory         c
c       (in the maximally abelian gauge).                             c
c---------------------------------------------------------------------c
c
c     * Include file
c         paravp3
c
c     * The following subroutine is only available on VP and VPP.
c         clockv
c         gettod
c
c     * This routine calls
c         Seed, Dirct, Dir, inits, Monte, (Mafix1, Dcuxy,
c         Measure).
c
c     * Input parameters
c         b      : 6/(g*g). g is a coupling constant.
c         init   : 1 = start. 0 = continue.
c         ints   : 0 = hot start. 1 = cold start.
c         iter   : the number of MC sweeps for every measurement.
c         nav    : Iter * Nav is the number of thermalization sweeps.
c         niter  : the number of vacuum configuration.
c         timemax: the time limit (micro sec.).
c         nhit   : see subroutine Monte.
c         ngfmax : see subroutine Mafix1.
c         epsrzz : see subroutine Mafix1.
c         omega  : see subroutine Mafix1.
c       The file, inputvp3, describes a example of these parameters.
c
c       S.Kitahara (98.10.14)
c----------------------------------------------------------------------

      include'paravp3'


      parameter(nbet=30)
      character*2 num(100)
      character*5 numc

      data num
     *    /'00','01','02','03','04','05','06','07','08','09',
     *     '10','11','12','13','14','15','16','17','18','19',
     *     '20','21','22','23','24','25','26','27','28','29',
     *     '30','31','32','33','34','35','36','37','38','39','40',
     *     '41','42','43','44','45','46','47','48','49','50',
     *     '51','52','53','54','55','56','57','58','59','60',
     *     '61','62','63','64','65','66','67','68','69','70',
     *     '71','72','73','74','75','76','77','78','79','80',
     *     '81','82','83','84','85','86','87','88','89','90',
     *     '91','92','93','94','95','96','97','98','99'/

      common /rndk /k1,k2
      common /var3/ z(nsite,3,3,nd)
      common /var2/ u(nvect,3,3,nd),v(nvect,3,3,nd)
      common /functional/ rmax
      dimension am(nsite,3,3,nd),fi(nsite,3)
c     dimension fi(nsite,3)
      dimension nspsw1(nbet)
c     data nspsw1/10,7*4,4*6,8,15,25,50,30,15,12*10/ !325
c     data nspsw1/10,7*4,4*6,8,12,20,40,25,12,12*8/ !275
      data nspsw1/2*10,11,12,13,14,16,17,20,23,25,31,37,47,70,230,132,
     *            98,84,74,68,63,60,57,55,53,52,51,50,50/
      real*8 time0,time1,time2
      character(len=200) :: conf_path

*      call gettod(time0)
      open(1,status='old',file='inputvp3')
      read(1,*) b
      read(1,*) init
      read(1,*) ints
      read(1,*) iter
      read(1,*) nav
      read(1,*) niter
      read(1,*) timemax
      read(1,*) nhit
      read(1,*) ngfmax
      read(1,*) epsrzz
      read(1,*) omega
      read(1,*) nrnd
      read(1,*) nstart
      read(1,*) ninit
      read(1,*) nstep
      read(1,*) coeff
      read(1,*) nover
      close(1)

      open(2,status='unknown',file='fmax.res')
c     open(3,status='unknown',file='wloop.dat')
      open(4,status='unknown',file='fmax.dat')
      open(6,status='unknown',file='mess.res')
c     open(7,status='unknown',file='ploop.dat')
      open(8,status='unknown',file='monop.dat')
      open(10,status='unknown',file='plaq.dat')
c     open(12,status='unknown',file='plc.dat')
c     open(14,status='unknown',file='plc_off.dat')
      open(15,status='unknown',file='functional-temperature.dat')
      open(17,status='unknown',file='wrap.dat')
      open(18,status='unknown',file='cluster.dat')
      write(2,300) n1,n2,n3,n4,b,init,ints,nav,iter,niter,nhit
  300 format(1x,'*******************************************'/
     c       1x,' Monte carlo simulation of SU3 LGT '/
     c       1x,'   '/
     c       1x,'*******************************************'/
     c       1x,'lattice size                :',i2,' x',i2,' x',i2
     c                                            ,' x',i2/
     c       1x,'beta                        :',f12.6/
     c       1x,'init (1=start, 0=continue)  :',i8/
     c       1x,'ints (0=hot, 1=cold)        :',i8/
     c       1x,'nav                         :',i8/
     c       1x,'iter                        :',i8/
     c       1x,'niter (# of samples)        :',i8/
     c       1x,'nhit                        :',i8/
     c       1x,'*******************************************')

c     preparation of random numbers, list vectors, constants
      k1=0
      k2=147
      call seed
      call dirct
      call dir
      pi2=8.e0*datan(1.d0)
      pi =4.e0*datan(1.d0)

c**********************************************
c     generate a thermalized configuration
      if(init.eq.1)then
        call inits(ints,pi2)
        do iav=1,nav
          call monte(iter,nhit,trace,b,pi2)
          write(2,*)'nsweep=',iav
        enddo
      else
c       input file is written here.
      endif
c***********************************************
c     to read thermalized configuration
c      if(init.eq.0)then
c        conf_path='/home/ilya/soft/lattice'//
c     *   '/general_code/tests/confs/SU3_conf/'//
c     *   'nt6/conf.0501_test'
c        conf_path='conf.0501_test'
c       OPEN(1,file='../../CONFIGS/CON000.LAT',
c        OPEN(1,file=conf_path,access='stream',
c     *form='unformatted',status='unknown')
c         read(1) u,v
c         read(1) u,v,k1,k2
c          read(1) z
c        CLOSE(1)
c        call plaq
c      endif
c      stop
c----------------------------------------------------------------------
c     main loop start ------------------------------------------

      print *, 'label 1'
      do nnn = ninit,nstep*(niter-1)+ninit,nstep
c     do nnn = ninit,niter-1+ninit
c       numc=num(nnn)
        write(numc,'(i5.5)') nnn
        conf_path='../confs/'//
     *   'nt6/conf.0501_test'
c        conf_path='conf.0501_test'
c       OPEN(1,file='../../CONFIGS/CON000.LAT',
        OPEN(1,file=conf_path,access='stream',
     *form='unformatted',status='unknown')
c       OPEN(1,file='../../CONFIGS/CON0'//numc//'.LAT',
c       OPEN(1,file='CON0'//numc//'.LAT',
c    *      form='unformatted',status='unknown')
c         read(1) z,k1,k2,rand
         read(1) z

c       OPEN(1,file=
c     *'/scr/tokyo/public/configurations/BQCD/16x16x16x4/b_5.25/k_0.13605
c    */bqcd.001.1.1.'//numc//'.lat',
c    *'bqcd.001.1.1.'//numc//'.lat',
c    *      form='unformatted',status='unknown')
c         read(1) z
c       CLOSE(1)
c     print *, 'label 2'
c----------------------------------------------------------------------
c for temporal use only:
c     do ir=1,111111
c       call rndpr3(4)
c     enddo
c random gauge transformation:
c     do jr=1,5
c        call rnd_gtr1
c        call rnd_gtr2
c        call rnd_gtr3
c     enddo
c----------------------------------------------------------------------
        call plaq
      print *, 'label 3'
c----------------------------------------------------------------------
c iterative gauge fixing:
c----------------------------------------------------------------------
c       call dcuxy3          ! to rewrite usual into even-odd
c       call mafix1(ngfmax,igf_fin,omega,epsrzz)
c       go to 10
c       call dcuxy
c       call measure(b,pi,pi2)
c----------------------------------------------------------------------
c SA gauge fixing:
c----------------------------------------------------------------------
c     go to 1
c     call dcuxy2          ! to rewrite even-odd into usual
c     nswnst=nspsw1(nstart)
c     nswnst=nspsw1(nstart)/2
      nswnst=nspsw1(nstart)*coeff
      rmax0=0.
      do irnd = 1,nrnd
c to start SA gauge fixing
        di=1.5d0/nbet
        do item=1,nbet
          temp0=1.6-item*di
          write(2,*) 'temp0=',temp0
          if(item.lt.nstart) nsw=0
c         if(item.eq.nstart) nsw=50
          if(item.eq.nstart) nsw=50*coeff
c         if(item.gt.nstart) nsw=nspsw1(item)
c         if(item.gt.nstart) nsw=nspsw1(item)/2
          if(item.gt.nstart) nsw=nspsw1(item)*coeff
          write(2,*) 'temp0=',temp0,' nsw=',nsw
          if(item.eq.nstart) then
           dtemp=di/nswnst
          else
           dtemp=di/nsw
          endif
        do ifsw=1,nsw
          if(item.eq.nstart) then
           idif=nsw-ifsw
           if(idif.ge.nswnst) then
            temp=temp0
           else
            temp=temp0-dtemp*(nswnst-(nsw-ifsw))
           endif
          else
           temp=temp0-dtemp*(ifsw-1)
          endif
c 1st SU(2) subgroup:
         call order1(am)
          do m=1,nsite
           fi(m,1)=0.d0
           fi(m,2)=0.d0
           fi(m,3)=1.d0
          enddo
           call sweep_hb(am,fi,temp)
          do isw=1,nover
           call sweep_or(am,fi)
          enddo
         call gtr1(fi)


c 2nd SU(2) subgroup:
         call order2(am)
          do m=1,nsite
           fi(m,1)=0.d0
           fi(m,2)=0.d0
           fi(m,3)=1.d0
          enddo
           call sweep_hb(am,fi,temp)
          do isw=1,nover
           call sweep_or(am,fi)
          enddo
         call gtr2(fi)


c 3rd SU(2) subgroup:
         call order3(am)
          do m=1,nsite
           fi(m,1)=0.d0
           fi(m,2)=0.d0
           fi(m,3)=1.d0
          enddo
           call sweep_hb(am,fi,temp)
          do isw=1,nover
           call sweep_or(am,fi)
          enddo
         call gtr3(fi)


         call sp_action(temp)
        enddo
        enddo
c       call measure(b,pi,pi2)   ! to check gauge invar.
        call dcuxy3              ! to rewrite usual lat. into odd-even
        igf_fin=0
        call mafix1(ngfmax,igf_fin,omega,epsrzz)
        call dcuxy              ! to rewrite odd-even into usual
c       call measure(b,pi,pi2)   ! to check gauge invar.

        write(2,20)nnn,igf_fin,trace
c       call mag_proj
c       numc=num(nnn)
c       call save_ab(numc,irnd)
c       call monop_def
c       do icolor=1,3
c        write(2,*) 'mlnew5 starts'
c        call mlnew5(icolor)
c        write(2,*) 'mlnew5 finished'
c       enddo
      call mag_proj
      call monop_def
      do icolor=1,3
       write(2,*) 'mlnew5 starts'
       call MLNEW5(icolor)
       write(2,*) 'mlnew5 finished'
      enddo
c to save best copy
       if(rmax.gt.rmax0) then
        rmax0=rmax
        call save_ab(numc,irnd)
        OPEN(1,file=
     *'bqcd_mag.'//numc//'.lat',
     *      form='unformatted',status='unknown')
          write(1) z
        CLOSE(1)
       endif
      enddo                   ! loop over gauge copies

c random gauge transformation:
       if(irnd.lt.nrnd) then
         call rnd_gtr1
         call rnd_gtr2
         call rnd_gtr3
c       call measure(b,pi,pi2)   ! to check gauge invar.
       endif
c to restore MC configuration:
c       OPEN(1,file='../CON.LAT',
c    *      form='unformatted',status='unknown')
c         read(1) u,v
c       CLOSE(1)


 1    continue
      enddo
c     main loop end --------------------------------------------
 10   continue

  500 continue

*      call clockv(time1,time2,2,2)
*      write(*,*)' cpu,vpu,ratio ',time1/1.e6,time2/1.e6,time1/time2
   20   format('nnn=',i4,' igf=',i4,' trace=',f12.5)

      close(2)
      stop
      end
c*********************************************************************c
      subroutine monte(iter,nhit,trace,b,pi2)
c
c     * for SU(3) Lattice on VP.
c     --------------------------
c
c     * This routine generates a equilibrated configurations(/var2/).
c       using heat bath method.
c
c         action = 1-1/2 Tr\sum_{n,mu,nu} UUUU
c
c     * Before this routine is called, you must once call
c         Dirct.
c
c     * This routine calls
c         Rndpr3.
c
c     * Include file
c         paravp3
c
c     * Input arrays, parameters and constants.
c         /var2/ : input configuration.
c         /tabl/ : the list vectors which are returned by Dirct.
c         iter   : the number of sweeps.
c         nhit   : the number of trial in order to select a component
c                  of SU(2) matrix (usually =1).
c         b      : 6/(g*g). g is the coupling constant.
c         pi2    : 2 times pi (6.28318...).
c
c     * Output arrays and variables
c         /var2/ : updated configuration.
c         trace  : the value of the plaquette action.
c
c     * Work arrays
c         /ran1/ : random numbers generated by Rndpr3.
c
c     * programmed by KEK group.
c     * modified by S.Kitahara. (98.10.14)
c----------------------------------------------------------------------

      include'paravp3'

      common /var2/ u(nvect,3,3,nd),v(nvect,3,3,nd)
      common /ran1/ rm1(nvect),rm2(nvect),rm3(nvect),rm4(nvect)
      common /tabl /mo1(nvect,nd), mo2(nvect,nd),
     &              me1(nvect,nd), me2(nvect,nd),
     &              ldir2(nd-1,nd)

      common /wrk2/ w1(nvect,3,3,nd,nd),  w2(nvect,3,3,nd,nd),
     c              w3(nvect,3,3),w(nvect,3,3),
     c              a0(nvect),dxi2(nvect),delta(nvect)
c     dimension     w1(nvect,3,3,nd,nd),  w2(nvect,3,3,nd,nd),
c    c              w3(nvect,3,3),w(nvect,3,3),
c    c              a0(nvect),dxi2(nvect),delta(nvect)
      dimension     jhit(nvect)

      b3=b/3.e0
      hitp=0.e0

      do 1000 nit=1,iter
        ztr=(0.e0,0.e0)

c   ------- u  -------------------------------------------
      do nu = 2,nd
      do mu = 1,nu-1
      do j  = 1,3
      do i  = 1,3
        do m = 1,nvect
          m2 = mo1(m,mu)
          m4 = mo1(m,nu)
          m6 = mo2(m,nu)
          w1(m,i,j,mu,nu) = v(m2,i,1,nu)*dconjg(v(m4,j,1,mu))
     c                    + v(m2,i,2,nu)*dconjg(v(m4,j,2,mu))
     c                    + v(m2,i,3,nu)*dconjg(v(m4,j,3,mu))
          w2(m,i,j,mu,nu) = dconjg(v(m6,1,i,mu))*v(m6,1,j,nu)
     c                    + dconjg(v(m6,2,i,mu))*v(m6,2,j,nu)
     c                    + dconjg(v(m6,3,i,mu))*v(m6,3,j,nu)
        enddo
      enddo
      enddo
      enddo
      enddo
c
      do nu = 1,nd-1
      do mu = nu+1,nd
      do i  = 1,3
      do j  = 1,3
        do m = 1,nvect
          m5=me2(mo1(m,mu),nu)
          w1(m,i,j,mu,nu) = dconjg(w1(m,j,i,nu,mu))
          w3(m,i,j)       = dconjg(w2(m5,j,i,nu,mu))
        enddo
        do m = 1,nvect
          w2(m,i,j,mu,nu) = w3(m,i,j)
        enddo
      enddo
      enddo
      enddo
      enddo

      do 20 mu = 1,nd

      do j  = 1,3
      do i  = 1,3
        do m = 1,nvect
          w(m,i,j) = 0.e0
        enddo
        do nu0 = 1,nd-1
          nu = ldir2(nu0,mu)
          do m = 1,nvect
            m5=me2(mo1(m,mu),nu)
            w(m,i,j) = w(m,i,j)
     c               + w1(m,i,1,mu,nu)*dconjg(u(m,j,1,nu))
     c               + w1(m,i,2,mu,nu)*dconjg(u(m,j,2,nu))
     c               + w1(m,i,3,mu,nu)*dconjg(u(m,j,3,nu))
     c               + dconjg(u(m5,1,i,nu))*w2(m,1,j,mu,nu)
     c               + dconjg(u(m5,2,i,nu))*w2(m,2,j,mu,nu)
     c               + dconjg(u(m5,3,i,nu))*w2(m,3,j,mu,nu)
          enddo
        enddo
      enddo
      enddo
c - - - - - - - u(1,1)-u(2,2) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=1,2
      do j=1,2
      do m=1,nvect
        w3(m,i,j)=u(m,i,1,mu)*w(m,1,j)
     c           +u(m,i,2,mu)*w(m,2,j)
     c           +u(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww11=w3(m,1,1)+dconjg(w3(m,2,2))
        ww12=w3(m,1,2)-dconjg(w3(m,2,1))
        ww21=w3(m,2,1)-dconjg(w3(m,1,2))
        ww22=w3(m,2,2)+dconjg(w3(m,1,1))
        dxi2(m)=dsqrt(dble(ww11*ww22-ww12*ww21))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 29 m = 1,nvect
          if(jhit(m).eq.1) goto 29
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 29
            a0(m)=1.e0-delta(m)
            jhit(m)=1
  29    continue
      enddo

      call rndpr3(2)

      do m=1,nvect
        if (jhit(m).eq.1)then
          hitp=hitp+1.e0
          c0 =dble (w3(m,1,1)+w3(m,2,2))/dxi2(m)
          c1 =dimag(w3(m,1,2)+w3(m,2,1))/dxi2(m)
          c2 =dble (w3(m,1,2)-w3(m,2,1))/dxi2(m)
          c3 =dimag(w3(m,1,1)-w3(m,2,2))/dxi2(m)
          rad = 1.e0 - a0(m)**2
          a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
          rad2 = dsqrt(dabs(rad-a3**2))
          theta = pi2*rm2(m)
          a1 = rad2*dcos(theta)
          a2 = rad2*dsin(theta)
          d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
          d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
          d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
          d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
          y11=dcmplx( d0, d3)
          y12=dcmplx( d2, d1)
          y21=dcmplx(-d2, d1)
          y22=dcmplx( d0,-d3)
          u11=u(m,1,1,mu)
          u12=u(m,1,2,mu)
          u13=u(m,1,3,mu)
          u21=u(m,2,1,mu)
          u22=u(m,2,2,mu)
          u23=u(m,2,3,mu)
          u(m,1,1,mu)=y11*u11+y12*u21
          u(m,1,2,mu)=y11*u12+y12*u22
          u(m,1,3,mu)=y11*u13+y12*u23
          u(m,2,1,mu)=y21*u11+y22*u21
          u(m,2,2,mu)=y21*u12+y22*u22
          u(m,2,3,mu)=y21*u13+y22*u23
        endif
      enddo

c - - - - - - - u(2,2)-u(3,3) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=2,3
      do j=2,3
      do m=1,nvect
        w3(m,i,j)=u(m,i,1,mu)*w(m,1,j)
     c           +u(m,i,2,mu)*w(m,2,j)
     c           +u(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww22=w3(m,2,2)+dconjg(w3(m,3,3))
        ww23=w3(m,2,3)-dconjg(w3(m,3,2))
        ww32=w3(m,3,2)-dconjg(w3(m,2,3))
        ww33=w3(m,3,3)+dconjg(w3(m,2,2))
        dxi2(m)=dsqrt(dble(ww22*ww33-ww23*ww32))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 129 m = 1,nvect
          if(jhit(m).eq.1) goto 129
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 129
            a0(m)=1.e0-delta(m)
            jhit(m)=1
 129    continue
      enddo

      call rndpr3(2)

      do m=1,nvect
        if (jhit(m).eq.1)then
          hitp=hitp+1.e0
          c0 =dble (w3(m,2,2)+w3(m,3,3))/dxi2(m)
          c1 =dimag(w3(m,2,3)+w3(m,3,2))/dxi2(m)
          c2 =dble (w3(m,2,3)-w3(m,3,2))/dxi2(m)
          c3 =dimag(w3(m,2,2)-w3(m,3,3))/dxi2(m)
          rad = 1.e0 - a0(m)**2
          a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
          rad2 = dsqrt(dabs(rad-a3**2))
          theta = pi2*rm2(m)
          a1 = rad2*dcos(theta)
          a2 = rad2*dsin(theta)
          d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
          d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
          d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
          d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
          y22=dcmplx( d0, d3)
          y23=dcmplx( d2, d1)
          y32=dcmplx(-d2, d1)
          y33=dcmplx( d0,-d3)
          u21=u(m,2,1,mu)
          u22=u(m,2,2,mu)
          u23=u(m,2,3,mu)
          u31=u(m,3,1,mu)
          u32=u(m,3,2,mu)
          u33=u(m,3,3,mu)
          u(m,2,1,mu)=y22*u21+y23*u31
          u(m,2,2,mu)=y22*u22+y23*u32
          u(m,2,3,mu)=y22*u23+y23*u33
          u(m,3,1,mu)=y32*u21+y33*u31
          u(m,3,2,mu)=y32*u22+y33*u32
          u(m,3,3,mu)=y32*u23+y33*u33
        endif
      enddo

c - - - - - - - u(1,1)-u(3,3) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=1,3,2
      do j=1,3,2
      do m=1,nvect
        w3(m,i,j)=u(m,i,1,mu)*w(m,1,j)
     c           +u(m,i,2,mu)*w(m,2,j)
     c           +u(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww11=w3(m,1,1)+dconjg(w3(m,3,3))
        ww13=w3(m,1,3)-dconjg(w3(m,3,1))
        ww31=w3(m,3,1)-dconjg(w3(m,1,3))
        ww33=w3(m,3,3)+dconjg(w3(m,1,1))
        dxi2(m)=dsqrt(dble(ww11*ww33-ww13*ww31))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 229 m = 1,nvect
          if(jhit(m).eq.1) goto 229
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 229
            a0(m)=1.e0-delta(m)
            jhit(m)=1
 229    continue
      enddo

      call rndpr3(2)

      do m=1,nvect
        if (jhit(m).eq.1)then
          hitp=hitp+1.e0
          c0 =dble (w3(m,1,1)+w3(m,3,3))/dxi2(m)
          c1 =dimag(w3(m,1,3)+w3(m,3,1))/dxi2(m)
          c2 =dble (w3(m,1,3)-w3(m,3,1))/dxi2(m)
          c3 =dimag(w3(m,1,1)-w3(m,3,3))/dxi2(m)
          rad = 1.e0 - a0(m)**2
          a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
          rad2 = dsqrt(dabs(rad-a3**2))
          theta = pi2*rm2(m)
          a1 = rad2*dcos(theta)
          a2 = rad2*dsin(theta)
          d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
          d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
          d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
          d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
          y11=dcmplx( d0, d3)
          y13=dcmplx( d2, d1)
          y31=dcmplx(-d2, d1)
          y33=dcmplx( d0,-d3)
          u11=u(m,1,1,mu)
          u12=u(m,1,2,mu)
          u13=u(m,1,3,mu)
          u31=u(m,3,1,mu)
          u32=u(m,3,2,mu)
          u33=u(m,3,3,mu)
          u(m,1,1,mu)=y11*u11+y13*u31
          u(m,1,2,mu)=y11*u12+y13*u32
          u(m,1,3,mu)=y11*u13+y13*u33
          u(m,3,1,mu)=y31*u11+y33*u31
          u(m,3,2,mu)=y31*u12+y33*u32
          u(m,3,3,mu)=y31*u13+y33*u33
        endif
      enddo


      do m=1,nvect
        znrm1=u(m,1,1,mu)*dconjg(u(m,1,1,mu))
     c       +u(m,1,2,mu)*dconjg(u(m,1,2,mu))
     c       +u(m,1,3,mu)*dconjg(u(m,1,3,mu))
*        znrm1=cdsqrt(znrm1)
        znrm1=sqrt(znrm1)
        u(m,1,1,mu)=u(m,1,1,mu)/znrm1
        u(m,1,2,mu)=u(m,1,2,mu)/znrm1
        u(m,1,3,mu)=u(m,1,3,mu)/znrm1
        zipd21=u(m,2,1,mu)*dconjg(u(m,1,1,mu))
     c        +u(m,2,2,mu)*dconjg(u(m,1,2,mu))
     c        +u(m,2,3,mu)*dconjg(u(m,1,3,mu))
        u(m,2,1,mu)=u(m,2,1,mu)-zipd21*u(m,1,1,mu)
        u(m,2,2,mu)=u(m,2,2,mu)-zipd21*u(m,1,2,mu)
        u(m,2,3,mu)=u(m,2,3,mu)-zipd21*u(m,1,3,mu)
        znrm2=u(m,2,1,mu)*dconjg(u(m,2,1,mu))
     c       +u(m,2,2,mu)*dconjg(u(m,2,2,mu))
     c       +u(m,2,3,mu)*dconjg(u(m,2,3,mu))
*        znrm2=cdsqrt(znrm2)
        znrm2=sqrt(znrm2)
        u(m,2,1,mu)=u(m,2,1,mu)/znrm2
        u(m,2,2,mu)=u(m,2,2,mu)/znrm2
        u(m,2,3,mu)=u(m,2,3,mu)/znrm2
        zipd31=u(m,3,1,mu)*dconjg(u(m,1,1,mu))
     c        +u(m,3,2,mu)*dconjg(u(m,1,2,mu))
     c        +u(m,3,3,mu)*dconjg(u(m,1,3,mu))
        zipd32=u(m,3,1,mu)*dconjg(u(m,2,1,mu))
     c        +u(m,3,2,mu)*dconjg(u(m,2,2,mu))
     c        +u(m,3,3,mu)*dconjg(u(m,2,3,mu))
        u(m,3,1,mu)=u(m,3,1,mu)-zipd32*u(m,2,1,mu)-zipd31*u(m,1,1,mu)
        u(m,3,2,mu)=u(m,3,2,mu)-zipd32*u(m,2,2,mu)-zipd31*u(m,1,2,mu)
        u(m,3,3,mu)=u(m,3,3,mu)-zipd32*u(m,2,3,mu)-zipd31*u(m,1,3,mu)
        znrm3=u(m,3,1,mu)*dconjg(u(m,3,1,mu))
     c       +u(m,3,2,mu)*dconjg(u(m,3,2,mu))
     c       +u(m,3,3,mu)*dconjg(u(m,3,3,mu))
*        znrm3=cdsqrt(znrm3)
        znrm3=sqrt(znrm3)
        u(m,3,1,mu)=u(m,3,1,mu)/znrm3
        u(m,3,2,mu)=u(m,3,2,mu)/znrm3
        u(m,3,3,mu)=u(m,3,3,mu)/znrm3
      enddo

      do i=1,3
      do m=1,nvect
        ztr=ztr+w(m,i,1)*u(m,1,i,mu)
     c         +w(m,i,2)*u(m,2,i,mu)
     c         +w(m,i,3)*u(m,3,i,mu)
      enddo
      enddo

  20  continue


c   ------- v  -------------------------------------------
      do nu = 2,nd
      do mu = 1,nu-1
      do j  = 1,3
      do i  = 1,3
        do m = 1,nvect
          m2 = me1(m,mu)
          m4 = me1(m,nu)
          m6 = me2(m,nu)
          w1(m,i,j,mu,nu) = u(m2,i,1,nu)*dconjg(u(m4,j,1,mu))
     c                    + u(m2,i,2,nu)*dconjg(u(m4,j,2,mu))
     c                    + u(m2,i,3,nu)*dconjg(u(m4,j,3,mu))
          w2(m,i,j,mu,nu) = dconjg(u(m6,1,i,mu))*u(m6,1,j,nu)
     c                    + dconjg(u(m6,2,i,mu))*u(m6,2,j,nu)
     c                    + dconjg(u(m6,3,i,mu))*u(m6,3,j,nu)
        enddo
      enddo
      enddo
      enddo
      enddo

      do nu = 1,nd-1
      do mu = nu+1,nd
      do i  = 1,3
      do j  = 1,3
        do m = 1,nvect
          m5=mo2(me1(m,mu),nu)
          w1(m,i,j,mu,nu) = dconjg(w1(m,j,i,nu,mu))
          w3(m,i,j)       = dconjg(w2(m5,j,i,nu,mu))
        enddo
        do m = 1,nvect
          w2(m,i,j,mu,nu) = w3(m,i,j)
        enddo
      enddo
      enddo
      enddo
      enddo

      do 320 mu = 1,nd

      do j  = 1,3
      do i  = 1,3
        do m = 1,nvect
          w(m,i,j) = 0.e0
        enddo
        do nu0 = 1,nd-1
          nu = ldir2(nu0,mu)
          do m = 1,nvect
            m5=mo2(me1(m,mu),nu)
            w(m,i,j) = w(m,i,j)
     c               + w1(m,i,1,mu,nu)*dconjg(v(m,j,1,nu))
     c               + w1(m,i,2,mu,nu)*dconjg(v(m,j,2,nu))
     c               + w1(m,i,3,mu,nu)*dconjg(v(m,j,3,nu))
     c               + dconjg(v(m5,1,i,nu))*w2(m,1,j,mu,nu)
     c               + dconjg(v(m5,2,i,nu))*w2(m,2,j,mu,nu)
     c               + dconjg(v(m5,3,i,nu))*w2(m,3,j,mu,nu)
          enddo
        enddo
      enddo
      enddo
c - - - - - - - v(1,1)-v(2,2) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=1,2
      do j=1,2
      do m=1,nvect
        w3(m,i,j)=v(m,i,1,mu)*w(m,1,j)
     c           +v(m,i,2,mu)*w(m,2,j)
     c           +v(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww11=w3(m,1,1)+dconjg(w3(m,2,2))
        ww12=w3(m,1,2)-dconjg(w3(m,2,1))
        ww21=w3(m,2,1)-dconjg(w3(m,1,2))
        ww22=w3(m,2,2)+dconjg(w3(m,1,1))
        dxi2(m)=dsqrt(dble(ww11*ww22-ww12*ww21))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 329 m = 1,nvect
          if(jhit(m).eq.1) goto 329
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 329
            a0(m)=1.e0-delta(m)
            jhit(m)=1
 329    continue
      enddo
c
      call rndpr3(2)
c
      do m=1,nvect
        if (jhit(m).eq.1)then
          hitp=hitp+1.e0
          c0 =dble (w3(m,1,1)+w3(m,2,2))/dxi2(m)
          c1 =dimag(w3(m,1,2)+w3(m,2,1))/dxi2(m)
          c2 =dble (w3(m,1,2)-w3(m,2,1))/dxi2(m)
          c3 =dimag(w3(m,1,1)-w3(m,2,2))/dxi2(m)
          rad = 1.e0 - a0(m)**2
          a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
          rad2 = dsqrt(dabs(rad-a3**2))
          theta = pi2*rm2(m)
          a1 = rad2*dcos(theta)
          a2 = rad2*dsin(theta)
          d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
          d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
          d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
          d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
          y11=dcmplx( d0, d3)
          y12=dcmplx( d2, d1)
          y21=dcmplx(-d2, d1)
          y22=dcmplx( d0,-d3)
          v11=v(m,1,1,mu)
          v12=v(m,1,2,mu)
          v13=v(m,1,3,mu)
          v21=v(m,2,1,mu)
          v22=v(m,2,2,mu)
          v23=v(m,2,3,mu)
          v(m,1,1,mu)=y11*v11+y12*v21
          v(m,1,2,mu)=y11*v12+y12*v22
          v(m,1,3,mu)=y11*v13+y12*v23
          v(m,2,1,mu)=y21*v11+y22*v21
          v(m,2,2,mu)=y21*v12+y22*v22
          v(m,2,3,mu)=y21*v13+y22*v23
        endif
      enddo

c - - - - - - - v(2,2)-v(3,3) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=2,3
      do j=2,3
      do m=1,nvect
        w3(m,i,j)=v(m,i,1,mu)*w(m,1,j)
     c           +v(m,i,2,mu)*w(m,2,j)
     c           +v(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww22=w3(m,2,2)+dconjg(w3(m,3,3))
        ww23=w3(m,2,3)-dconjg(w3(m,3,2))
        ww32=w3(m,3,2)-dconjg(w3(m,2,3))
        ww33=w3(m,3,3)+dconjg(w3(m,2,2))
        dxi2(m)=dsqrt(dble(ww22*ww33-ww23*ww32))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 429 m = 1,nvect
          if(jhit(m).eq.1) goto 429
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 429
            a0(m)=1.e0-delta(m)
            jhit(m)=1
 429    continue
      enddo

      call rndpr3(2)

      do m=1,nvect
        if (jhit(m).eq.1)then
        hitp=hitp+1.e0
        c0 =dble (w3(m,2,2)+w3(m,3,3))/dxi2(m)
        c1 =dimag(w3(m,2,3)+w3(m,3,2))/dxi2(m)
        c2 =dble (w3(m,2,3)-w3(m,3,2))/dxi2(m)
        c3 =dimag(w3(m,2,2)-w3(m,3,3))/dxi2(m)
        rad = 1.e0 - a0(m)**2
        a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
        rad2 = dsqrt(dabs(rad-a3**2))
        theta = pi2*rm2(m)
        a1 = rad2*dcos(theta)
        a2 = rad2*dsin(theta)
        d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
        d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
        d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
        d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
        y22=dcmplx( d0, d3)
        y23=dcmplx( d2, d1)
        y32=dcmplx(-d2, d1)
        y33=dcmplx( d0,-d3)
        v21=v(m,2,1,mu)
        v22=v(m,2,2,mu)
        v23=v(m,2,3,mu)
        v31=v(m,3,1,mu)
        v32=v(m,3,2,mu)
        v33=v(m,3,3,mu)
        v(m,2,1,mu)=y22*v21+y23*v31
        v(m,2,2,mu)=y22*v22+y23*v32
        v(m,2,3,mu)=y22*v23+y23*v33
        v(m,3,1,mu)=y32*v21+y33*v31
        v(m,3,2,mu)=y32*v22+y33*v32
        v(m,3,3,mu)=y32*v23+y33*v33
        endif
      enddo

c - - - - - - - v(1,1)-v(3,3) - - - - - - - - - - - - - -
      do m=1,nvect
        jhit(m)=0
      enddo

      do i=1,3,2
      do j=1,3,2
      do m=1,nvect
        w3(m,i,j)=v(m,i,1,mu)*w(m,1,j)
     c           +v(m,i,2,mu)*w(m,2,j)
     c           +v(m,i,3,mu)*w(m,3,j)
      enddo
      enddo
      enddo

      do m=1,nvect
        ww11=w3(m,1,1)+dconjg(w3(m,3,3))
        ww13=w3(m,1,3)-dconjg(w3(m,3,1))
        ww31=w3(m,3,1)-dconjg(w3(m,1,3))
        ww33=w3(m,3,3)+dconjg(w3(m,1,1))
        dxi2(m)=dsqrt(dble(ww11*ww33-ww13*ww31))
      enddo

      do ihit1=1,nhit
        call rndpr3(4)
        do 529 m = 1,nvect
          if(jhit(m).eq.1) goto 529
          delta(m)=-(dlog(rm1(m))+dlog(rm2(m))*dcos(pi2*rm3(m))**2)
     c             /(dxi2(m)*b3)
          if (rm4(m)**2.gt.1.e0-delta(m)/2.e0) goto 529
            a0(m)=1.e0-delta(m)
            jhit(m)=1
 529    continue
      enddo

      call rndpr3(2)

      do m=1,nvect
        if (jhit(m).eq.1)then
          hitp=hitp+1.e0
          c0 =dble (w3(m,1,1)+w3(m,3,3))/dxi2(m)
          c1 =dimag(w3(m,1,3)+w3(m,3,1))/dxi2(m)
          c2 =dble (w3(m,1,3)-w3(m,3,1))/dxi2(m)
          c3 =dimag(w3(m,1,1)-w3(m,3,3))/dxi2(m)
          rad = 1.e0 - a0(m)**2
          a3 = dsqrt(rad)*(2.e0*rm1(m)-1.e0)
          rad2 = dsqrt(dabs(rad-a3**2))
          theta = pi2*rm2(m)
          a1 = rad2*dcos(theta)
          a2 = rad2*dsin(theta)
          d0 = a0(m)*c0+a1*c1+a2*c2+a3*c3
          d1 =-a0(m)*c1+a1*c0+a2*c3-a3*c2
          d2 =-a0(m)*c2+a2*c0+a3*c1-a1*c3
          d3 =-a0(m)*c3+a3*c0+a1*c2-a2*c1
          y11=dcmplx( d0, d3)
          y13=dcmplx( d2, d1)
          y31=dcmplx(-d2, d1)
          y33=dcmplx( d0,-d3)
          v11=v(m,1,1,mu)
          v12=v(m,1,2,mu)
          v13=v(m,1,3,mu)
          v31=v(m,3,1,mu)
          v32=v(m,3,2,mu)
          v33=v(m,3,3,mu)
          v(m,1,1,mu)=y11*v11+y13*v31
          v(m,1,2,mu)=y11*v12+y13*v32
          v(m,1,3,mu)=y11*v13+y13*v33
          v(m,3,1,mu)=y31*v11+y33*v31
          v(m,3,2,mu)=y31*v12+y33*v32
          v(m,3,3,mu)=y31*v13+y33*v33
        endif
      enddo


      do m=1,nvect
        znrm1=v(m,1,1,mu)*dconjg(v(m,1,1,mu))
     c       +v(m,1,2,mu)*dconjg(v(m,1,2,mu))
     c       +v(m,1,3,mu)*dconjg(v(m,1,3,mu))
*        znrm1=cdsqrt(znrm1)
        znrm1=sqrt(znrm1)
        v(m,1,1,mu)=v(m,1,1,mu)/znrm1
        v(m,1,2,mu)=v(m,1,2,mu)/znrm1
        v(m,1,3,mu)=v(m,1,3,mu)/znrm1
        zipd21=v(m,2,1,mu)*dconjg(v(m,1,1,mu))
     c        +v(m,2,2,mu)*dconjg(v(m,1,2,mu))
     c        +v(m,2,3,mu)*dconjg(v(m,1,3,mu))
        v(m,2,1,mu)=v(m,2,1,mu)-zipd21*v(m,1,1,mu)
        v(m,2,2,mu)=v(m,2,2,mu)-zipd21*v(m,1,2,mu)
        v(m,2,3,mu)=v(m,2,3,mu)-zipd21*v(m,1,3,mu)
        znrm2=v(m,2,1,mu)*dconjg(v(m,2,1,mu))
     c       +v(m,2,2,mu)*dconjg(v(m,2,2,mu))
     c       +v(m,2,3,mu)*dconjg(v(m,2,3,mu))
*        znrm2=cdsqrt(znrm2)
        znrm2=sqrt(znrm2)
        v(m,2,1,mu)=v(m,2,1,mu)/znrm2
        v(m,2,2,mu)=v(m,2,2,mu)/znrm2
        v(m,2,3,mu)=v(m,2,3,mu)/znrm2
        zipd31=v(m,3,1,mu)*dconjg(v(m,1,1,mu))
     c        +v(m,3,2,mu)*dconjg(v(m,1,2,mu))
     c        +v(m,3,3,mu)*dconjg(v(m,1,3,mu))
        zipd32=v(m,3,1,mu)*dconjg(v(m,2,1,mu))
     c        +v(m,3,2,mu)*dconjg(v(m,2,2,mu))
     c        +v(m,3,3,mu)*dconjg(v(m,2,3,mu))
        v(m,3,1,mu)=v(m,3,1,mu)-zipd32*v(m,2,1,mu)-zipd31*v(m,1,1,mu)
        v(m,3,2,mu)=v(m,3,2,mu)-zipd32*v(m,2,2,mu)-zipd31*v(m,1,2,mu)
        v(m,3,3,mu)=v(m,3,3,mu)-zipd32*v(m,2,3,mu)-zipd31*v(m,1,3,mu)
        znrm3=v(m,3,1,mu)*dconjg(v(m,3,1,mu))
     c       +v(m,3,2,mu)*dconjg(v(m,3,2,mu))
     c       +v(m,3,3,mu)*dconjg(v(m,3,3,mu))
*        znrm3=cdsqrt(znrm3)
        znrm3=sqrt(znrm3)
        v(m,3,1,mu)=v(m,3,1,mu)/znrm3
        v(m,3,2,mu)=v(m,3,2,mu)/znrm3
        v(m,3,3,mu)=v(m,3,3,mu)/znrm3
      enddo

      do i=1,3
      do m=1,nvect
        ztr=ztr+w(m,i,1)*v(m,1,i,mu)
     c         +w(m,i,2)*v(m,2,i,mu)
     c         +w(m,i,3)*v(m,3,i,mu)
      enddo
      enddo

 320  continue


 1000 continue

      trace=dreal(ztr)/dble(3*6*nlink)
      hitp =hitp/dble(3*nlink*iter)

      return
      end
c*********************************************************************c
c*********************************************************************c
      subroutine mafix1(ngfmax,igf_fin,omega,epsrzz)
c
c     * for SU(3) Lattice on VP.
c     --------------------------
c
c     * A configuration /var2/ is fixed in the maximally abelian
c       gauge iteratively.
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
c         /var2/ : the thermalized configuration generated by Monte.
c         /tabl/ : the list vectors which are returned by Dirct.
c         ngfmax : the number of maximum sweeps.
c         omega  : over-relaxation parameter.
c         epsrzz : The gauge-fixing sweep is continued until
c                  Rzz is smaller than Epsrzz. Rzz is the square
c                  of the off diagonal part of the operator
c                  to be diagonalized.
c
c     * Output arrays and variables
c         igf    : the number of gauge-fixing sweeps.
c         /wrk1/ : the gauge-fixed configuration.
c
c     * programmed by S.Kitahara.
c     * modified by S.Kitahara. (98.10.15)
c----------------------------------------------------------------------

      include'paravp3'
      common /var2/  u(nvect,3,3,nd),v(nvect,3,3,nd)
      common /wrk1/  x(nvect,3,3,nd),y(nvect,3,3,nd)
      common /tb1/  mioe(nsite), mioo(nsite)
      common /tabl /mo1(nvect,nd), mo2(nvect,nd),
     &              me1(nvect,nd), me2(nvect,nd),
     &              ldir2(nd-1,nd)

      common /wrk2/ zx(nvect,3,3,nd),zy(nvect,3,3,nd),wg(nvect,3,3)
     c             ,qdi(nvect),qod(nvect),qre(nvect),qim(nvect)
      common /functional/ rmax
c     dimension     zx(nvect,3,3,nd),zy(nvect,3,3,nd),wg(nvect,3,3)
c    c             ,qdi(nvect),qod(nvect),qre(nvect),qim(nvect)


      do l=1,nd
      do i=1,3
      do j=1,3
      do m=1,nvect
          x(m,i,j,l)= u(m,i,j,l)
          y(m,i,j,l)= v(m,i,j,l)
      enddo
      enddo
      enddo
      enddo


      do 1000 igf=1,ngfmax

c     odd

      do m=1,nvect
        qdi(m)=0.e0
        qod(m)=0.e0
        qre(m)=0.e0
        qim(m)=0.e0
      enddo

      do l=1,nd
      do m=1,nvect
        qdi(m)=qdi(m)+dreal( x(m,1,1,l)*dconjg(x(m,1,1,l))
     c                      +x(m,2,2,l)*dconjg(x(m,2,2,l)) )
        qod(m)=qod(m)+dreal( x(m,2,1,l)*dconjg(x(m,2,1,l))
     c                      +x(m,1,2,l)*dconjg(x(m,1,2,l)) )
        zqri= x(m,1,1,l)*dconjg(x(m,2,1,l))
     c       -x(m,1,2,l)*dconjg(x(m,2,2,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do l=1,nd
      do m=1,nvect
        m0=mo2(m,l)
        qdi(m)=qdi(m)+dreal( y(m0,1,1,l)*dconjg(y(m0,1,1,l))
     c                      +y(m0,2,2,l)*dconjg(y(m0,2,2,l)) )
        qod(m)=qod(m)+dreal( y(m0,2,1,l)*dconjg(y(m0,2,1,l))
     c                      +y(m0,1,2,l)*dconjg(y(m0,1,2,l)) )
        zqri= y(m0,1,2,l)*dconjg(y(m0,1,1,l))
     c       -y(m0,2,2,l)*dconjg(y(m0,2,1,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do m=1,nvect
        ak1 = (qdi(m)-qod(m))*0.5e0
        ak2 = qre(m)**2+qim(m)**2
        ak  = ak1+sqrt(ak1**2 + ak2)
        tmp1= sqrt( ak**2+ak2 )
        q0 = ak/tmp1
        q1 = qim(m)/tmp1
        q2 = qre(m)/tmp1
c       if(q0.ne.1.e0)then
        if(abs(q0-1.).gt.1.d-11)then
          angl=acos(q0)
          tmp =sin( angl*omega )/sqrt(1.e0-q0*q0)
          q0 =cos( angl*omega )
          q1 =q1*tmp
          q2 =q2*tmp
        else
c        stop
          q0=1.
          q1=0.
          q2=0.
        endif
        wg(m,1,1)= dcmplx( q0)
        wg(m,1,2)= dcmplx( q2, q1)
        wg(m,2,1)= dcmplx(-q2, q1)
        wg(m,2,2)= dcmplx( q0)
      enddo

      do l=1,nd
      do i=1,2
      do j=1,3
      do m=1,nvect
        zx(m,i,j,l)=x(m,i,j,l)
        zy(m,j,i,l)=y(m,j,i,l)
      enddo
      enddo
      enddo
      enddo

      do l=1,nd
      do i=1,2
      do j=1,3
      do m=1,nvect
        m0=mo2(m,l)
        x(m,i,j,l)=wg(m,i,1)*zx(m,1,j,l)
     c            +wg(m,i,2)*zx(m,2,j,l)
        y(m0,j,i,l)=zy(m0,j,1,l)*dconjg(wg(m,i,1))
     c             +zy(m0,j,2,l)*dconjg(wg(m,i,2))
      enddo
      enddo
      enddo
      enddo


      do m=1,nvect
        qdi(m)=0.e0
        qod(m)=0.e0
        qre(m)=0.e0
        qim(m)=0.e0
      enddo

      do l=1,nd
      do m=1,nvect
        qdi(m)=qdi(m)+dreal( x(m,2,2,l)*dconjg(x(m,2,2,l))
     c                      +x(m,3,3,l)*dconjg(x(m,3,3,l)) )
        qod(m)=qod(m)+dreal( x(m,3,2,l)*dconjg(x(m,3,2,l))
     c                      +x(m,2,3,l)*dconjg(x(m,2,3,l)) )
        zqri= x(m,2,2,l)*dconjg(x(m,3,2,l))
     c       -x(m,2,3,l)*dconjg(x(m,3,3,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do l=1,nd
      do m=1,nvect
        m0=mo2(m,l)
        qdi(m)=qdi(m)+dreal( y(m0,2,2,l)*dconjg(y(m0,2,2,l))
     c                      +y(m0,3,3,l)*dconjg(y(m0,3,3,l)) )
        qod(m)=qod(m)+dreal( y(m0,3,2,l)*dconjg(y(m0,3,2,l))
     c                      +y(m0,2,3,l)*dconjg(y(m0,2,3,l)) )
        zqri= y(m0,2,3,l)*dconjg(y(m0,2,2,l))
     c       -y(m0,3,3,l)*dconjg(y(m0,3,2,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do m=1,nvect
        ak1 = (qdi(m)-qod(m))*0.5e0
        ak2 = qre(m)**2+qim(m)**2
        ak  = ak1+sqrt(ak1**2 + ak2)
        tmp1= sqrt( ak**2+ak2 )
        q0 = ak/tmp1
        q1 = qim(m)/tmp1
        q2 = qre(m)/tmp1
c       if(q0.ne.1.e0)then
        if(abs(q0-1.).gt.1.d-11)then
          angl=acos(q0)
          tmp =sin( angl*omega )/sqrt(1.e0-q0*q0)
          q0 =cos( angl*omega )
          q1 =q1*tmp
          q2 =q2*tmp
        else
c         stop
          q0=1.
          q1=0.
          q2=0.
        endif
        wg(m,2,2)= dcmplx( q0)
        wg(m,2,3)= dcmplx( q2, q1)
        wg(m,3,2)= dcmplx(-q2, q1)
        wg(m,3,3)= dcmplx( q0)
      enddo

      do l=1,nd
      do i=2,3
      do j=1,3
      do m=1,nvect
        zx(m,i,j,l)=x(m,i,j,l)
        zy(m,j,i,l)=y(m,j,i,l)
      enddo
      enddo
      enddo
      enddo

      do l=1,nd
      do i=2,3
      do j=1,3
      do m=1,nvect
        m0=mo2(m,l)
        x(m,i,j,l)=wg(m,i,2)*zx(m,2,j,l)
     c            +wg(m,i,3)*zx(m,3,j,l)
        y(m0,j,i,l)=zy(m0,j,2,l)*dconjg(wg(m,i,2))
     c             +zy(m0,j,3,l)*dconjg(wg(m,i,3))
      enddo
      enddo
      enddo
      enddo


      do m=1,nvect
        qdi(m)=0.e0
        qod(m)=0.e0
        qre(m)=0.e0
        qim(m)=0.e0
      enddo

      do l=1,nd
      do m=1,nvect
        qdi(m)=qdi(m)+dreal( x(m,1,1,l)*dconjg(x(m,1,1,l))
     c                      +x(m,3,3,l)*dconjg(x(m,3,3,l)) )
        qod(m)=qod(m)+dreal( x(m,3,1,l)*dconjg(x(m,3,1,l))
     c                      +x(m,1,3,l)*dconjg(x(m,1,3,l)) )
        zqri= x(m,1,1,l)*dconjg(x(m,3,1,l))
     c       -x(m,1,3,l)*dconjg(x(m,3,3,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do l=1,nd
      do m=1,nvect
        m0=mo2(m,l)
        qdi(m)=qdi(m)+dreal( y(m0,1,1,l)*dconjg(y(m0,1,1,l))
     c                      +y(m0,3,3,l)*dconjg(y(m0,3,3,l)) )
        qod(m)=qod(m)+dreal( y(m0,3,1,l)*dconjg(y(m0,3,1,l))
     c                      +y(m0,1,3,l)*dconjg(y(m0,1,3,l)) )
        zqri= y(m0,1,3,l)*dconjg(y(m0,1,1,l))
     c       -y(m0,3,3,l)*dconjg(y(m0,3,1,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do m=1,nvect
        ak1 = (qdi(m)-qod(m))*0.5e0
        ak2 = qre(m)**2+qim(m)**2
        ak  = ak1+sqrt(ak1**2 + ak2)
        tmp1= sqrt( ak**2+ak2 )
        q0 = ak/tmp1
        q1 = qim(m)/tmp1
        q2 = qre(m)/tmp1
c       if(q0.ne.1.e0)then
        if(abs(q0-1.).gt.1.d-11)then
          angl=acos(q0)
          tmp =sin( angl*omega )/sqrt(1.e0-q0*q0)
          q0 =cos( angl*omega )
          q1 =q1*tmp
          q2 =q2*tmp
        else
c         stop
          q0=1.
          q1=0.
          q2=0.
        endif
        wg(m,1,1)= q0
        wg(m,1,3)= dcmplx( q2, q1)
        wg(m,3,1)= dcmplx(-q2, q1)
        wg(m,3,3)= q0
      enddo

      do l=1,nd
      do i=1,3,2
      do j=1,3
      do m=1,nvect
        zx(m,i,j,l)=x(m,i,j,l)
        zy(m,j,i,l)=y(m,j,i,l)
      enddo
      enddo
      enddo
      enddo

      do l=1,nd
      do i=1,3,2
      do j=1,3
      do m=1,nvect
        m0=mo2(m,l)
        x(m,i,j,l)=wg(m,i,1)*zx(m,1,j,l)
     c            +wg(m,i,3)*zx(m,3,j,l)
        y(m0,j,i,l)=zy(m0,j,1,l)*dconjg(wg(m,i,1))
     c             +zy(m0,j,3,l)*dconjg(wg(m,i,3))
      enddo
      enddo
      enddo
      enddo

c     even

      do m=1,nvect
        qdi(m)=0.e0
        qod(m)=0.e0
        qre(m)=0.e0
        qim(m)=0.e0
      enddo

      do l=1,nd
      do m=1,nvect
        qdi(m)=qdi(m)+dreal( y(m,1,1,l)*dconjg(y(m,1,1,l))
     c                      +y(m,2,2,l)*dconjg(y(m,2,2,l)) )
        qod(m)=qod(m)+dreal( y(m,2,1,l)*dconjg(y(m,2,1,l))
     c                      +y(m,1,2,l)*dconjg(y(m,1,2,l)) )
        zqri= y(m,1,1,l)*dconjg(y(m,2,1,l))
     c       -y(m,1,2,l)*dconjg(y(m,2,2,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do l=1,nd
      do m=1,nvect
        m0=me2(m,l)
        qdi(m)=qdi(m)+dreal( x(m0,1,1,l)*dconjg(x(m0,1,1,l))
     c                      +x(m0,2,2,l)*dconjg(x(m0,2,2,l)) )
        qod(m)=qod(m)+dreal( x(m0,2,1,l)*dconjg(x(m0,2,1,l))
     c                      +x(m0,1,2,l)*dconjg(x(m0,1,2,l)) )
        zqri= x(m0,1,2,l)*dconjg(x(m0,1,1,l))
     c       -x(m0,2,2,l)*dconjg(x(m0,2,1,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do m=1,nvect
        ak1 = (qdi(m)-qod(m))*0.5e0
        ak2 = qre(m)**2+qim(m)**2
        ak  = ak1+sqrt(ak1**2 + ak2)
        tmp1= sqrt( ak**2+ak2 )
        q0 = ak/tmp1
        q1 = qim(m)/tmp1
        q2 = qre(m)/tmp1
c       if(q0.ne.1.e0)then
        if(abs(q0-1.).gt.1.d-11)then
          angl=acos(q0)
          tmp =sin( angl*omega )/sqrt(1.e0-q0*q0)
          q0 =cos( angl*omega )
          q1 =q1*tmp
          q2 =q2*tmp
        else
c         stop
          q0=1.
          q1=0.
          q2=0.
        endif
        wg(m,1,1)= dcmplx( q0)
        wg(m,1,2)= dcmplx( q2, q1)
        wg(m,2,1)= dcmplx(-q2, q1)
        wg(m,2,2)= dcmplx( q0)
      enddo

      do l=1,nd
      do i=1,2
      do j=1,3
      do m=1,nvect
        zy(m,i,j,l)=y(m,i,j,l)
        zx(m,j,i,l)=x(m,j,i,l)
      enddo
      enddo
      enddo
      enddo

      do l=1,nd
      do i=1,2
      do j=1,3
      do m=1,nvect
        m0=me2(m,l)
        y(m,i,j,l)=wg(m,i,1)*zy(m,1,j,l)
     c            +wg(m,i,2)*zy(m,2,j,l)
        x(m0,j,i,l)=zx(m0,j,1,l)*dconjg(wg(m,i,1))
     c             +zx(m0,j,2,l)*dconjg(wg(m,i,2))
      enddo
      enddo
      enddo
      enddo


      do m=1,nvect
        qdi(m)=0.e0
        qod(m)=0.e0
        qre(m)=0.e0
        qim(m)=0.e0
      enddo

      do l=1,nd
      do m=1,nvect
        qdi(m)=qdi(m)+dreal( y(m,2,2,l)*dconjg(y(m,2,2,l))
     c                      +y(m,3,3,l)*dconjg(y(m,3,3,l)) )
        qod(m)=qod(m)+dreal( y(m,3,2,l)*dconjg(y(m,3,2,l))
     c                      +y(m,2,3,l)*dconjg(y(m,2,3,l)) )
        zqri= y(m,2,2,l)*dconjg(y(m,3,2,l))
     c       -y(m,2,3,l)*dconjg(y(m,3,3,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do l=1,nd
      do m=1,nvect
        m0=me2(m,l)
        qdi(m)=qdi(m)+dreal( x(m0,2,2,l)*dconjg(x(m0,2,2,l))
     c                      +x(m0,3,3,l)*dconjg(x(m0,3,3,l)) )
        qod(m)=qod(m)+dreal( x(m0,3,2,l)*dconjg(x(m0,3,2,l))
     c                      +x(m0,2,3,l)*dconjg(x(m0,2,3,l)) )
        zqri= x(m0,2,3,l)*dconjg(x(m0,2,2,l))
     c       -x(m0,3,3,l)*dconjg(x(m0,3,2,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do m=1,nvect
        ak1 = (qdi(m)-qod(m))*0.5e0
        ak2 = qre(m)**2+qim(m)**2
        ak  = ak1+sqrt(ak1**2 + ak2)
        tmp1= sqrt( ak**2+ak2 )
        q0 = ak/tmp1
        q1 = qim(m)/tmp1
        q2 = qre(m)/tmp1
c       if(q0.ne.1.e0)then
        if(abs(q0-1.).gt.1.d-11)then
          angl=acos(q0)
          tmp =sin( angl*omega )/sqrt(1.e0-q0*q0)
          q0 =cos( angl*omega )
          q1 =q1*tmp
          q2 =q2*tmp
        else
c         stop
          q0=1.
          q1=0.
          q2=0.
        endif
        wg(m,2,2)= dcmplx( q0)
        wg(m,2,3)= dcmplx( q2, q1)
        wg(m,3,2)= dcmplx(-q2, q1)
        wg(m,3,3)= dcmplx( q0)
      enddo

      do l=1,nd
      do i=2,3
      do j=1,3
      do m=1,nvect
        zy(m,i,j,l)=y(m,i,j,l)
        zx(m,j,i,l)=x(m,j,i,l)
      enddo
      enddo
      enddo
      enddo

      do l=1,nd
      do i=2,3
      do j=1,3
      do m=1,nvect
        m0=me2(m,l)
        y(m,i,j,l)=wg(m,i,2)*zy(m,2,j,l)
     c            +wg(m,i,3)*zy(m,3,j,l)
        x(m0,j,i,l)=zx(m0,j,2,l)*dconjg(wg(m,i,2))
     c             +zx(m0,j,3,l)*dconjg(wg(m,i,3))
      enddo
      enddo
      enddo
      enddo


      do m=1,nvect
        qdi(m)=0.e0
        qod(m)=0.e0
        qre(m)=0.e0
        qim(m)=0.e0
      enddo

      do l=1,nd
      do m=1,nvect
        qdi(m)=qdi(m)+dreal( y(m,1,1,l)*dconjg(y(m,1,1,l))
     c                      +y(m,3,3,l)*dconjg(y(m,3,3,l)) )
        qod(m)=qod(m)+dreal( y(m,3,1,l)*dconjg(y(m,3,1,l))
     c                      +y(m,1,3,l)*dconjg(y(m,1,3,l)) )
        zqri= y(m,1,1,l)*dconjg(y(m,3,1,l))
     c       -y(m,1,3,l)*dconjg(y(m,3,3,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do l=1,nd
      do m=1,nvect
        m0=me2(m,l)
        qdi(m)=qdi(m)+dreal( x(m0,1,1,l)*dconjg(x(m0,1,1,l))
     c                      +x(m0,3,3,l)*dconjg(x(m0,3,3,l)) )
        qod(m)=qod(m)+dreal( x(m0,3,1,l)*dconjg(x(m0,3,1,l))
     c                      +x(m0,1,3,l)*dconjg(x(m0,1,3,l)) )
        zqri= x(m0,1,3,l)*dconjg(x(m0,1,1,l))
     c       -x(m0,3,3,l)*dconjg(x(m0,3,1,l))
        qre(m)=qre(m)+dreal( zqri)
        qim(m)=qim(m)+dimag( zqri)
      enddo
      enddo

      do m=1,nvect
        ak1 = (qdi(m)-qod(m))*0.5e0
        ak2 = qre(m)**2+qim(m)**2
        ak  = ak1+sqrt(ak1**2 + ak2)
        tmp1= sqrt( ak**2+ak2 )
        q0 = ak/tmp1
        q1 = qim(m)/tmp1
        q2 = qre(m)/tmp1
c       if(q0.ne.1.e0)then
        if(abs(q0-1.).gt.1.d-11)then
          angl=acos(q0)
          tmp =sin( angl*omega )/sqrt(1.e0-q0*q0)
          q0 =cos( angl*omega )
          q1 =q1*tmp
          q2 =q2*tmp
        else
c         stop
          q0=1.
          q1=0.
          q2=0.
        endif
        wg(m,1,1)= dcmplx( q0)
        wg(m,1,3)= dcmplx( q2, q1)
        wg(m,3,1)= dcmplx(-q2, q1)
        wg(m,3,3)= dcmplx( q0)
      enddo

      do l=1,nd
      do i=1,3,2
      do j=1,3
      do m=1,nvect
        zy(m,i,j,l)=y(m,i,j,l)
        zx(m,j,i,l)=x(m,j,i,l)
      enddo
      enddo
      enddo
      enddo

      do l=1,nd
      do i=1,3,2
      do j=1,3
      do m=1,nvect
        m0=me2(m,l)
        y(m,i,j,l)=wg(m,i,1)*zy(m,1,j,l)
     c            +wg(m,i,3)*zy(m,3,j,l)
        x(m0,j,i,l)=zx(m0,j,1,l)*dconjg(wg(m,i,1))
     c             +zx(m0,j,3,l)*dconjg(wg(m,i,3))
      enddo
      enddo
      enddo
      enddo

c     normalization

      do l=1,nd
      do m=1,nvect
        znrm1=x(m,1,1,l)*dconjg(x(m,1,1,l))
     c       +x(m,1,2,l)*dconjg(x(m,1,2,l))
     c       +x(m,1,3,l)*dconjg(x(m,1,3,l))
*        znrm1=cdsqrt(znrm1)
        znrm1=sqrt(znrm1)
        x(m,1,1,l)=x(m,1,1,l)/znrm1
        x(m,1,2,l)=x(m,1,2,l)/znrm1
        x(m,1,3,l)=x(m,1,3,l)/znrm1
        zipd21=x(m,2,1,l)*dconjg(x(m,1,1,l))
     c        +x(m,2,2,l)*dconjg(x(m,1,2,l))
     c        +x(m,2,3,l)*dconjg(x(m,1,3,l))
        x(m,2,1,l)=x(m,2,1,l)-zipd21*x(m,1,1,l)
        x(m,2,2,l)=x(m,2,2,l)-zipd21*x(m,1,2,l)
        x(m,2,3,l)=x(m,2,3,l)-zipd21*x(m,1,3,l)
        znrm2=x(m,2,1,l)*dconjg(x(m,2,1,l))
     c       +x(m,2,2,l)*dconjg(x(m,2,2,l))
     c       +x(m,2,3,l)*dconjg(x(m,2,3,l))
*        znrm2=cdsqrt(znrm2)
        znrm2=sqrt(znrm2)
        x(m,2,1,l)=x(m,2,1,l)/znrm2
        x(m,2,2,l)=x(m,2,2,l)/znrm2
        x(m,2,3,l)=x(m,2,3,l)/znrm2
        zipd31=x(m,3,1,l)*dconjg(x(m,1,1,l))
     c        +x(m,3,2,l)*dconjg(x(m,1,2,l))
     c        +x(m,3,3,l)*dconjg(x(m,1,3,l))
        zipd32=x(m,3,1,l)*dconjg(x(m,2,1,l))
     c        +x(m,3,2,l)*dconjg(x(m,2,2,l))
     c        +x(m,3,3,l)*dconjg(x(m,2,3,l))
        x(m,3,1,l)=x(m,3,1,l)-zipd32*x(m,2,1,l)-zipd31*x(m,1,1,l)
        x(m,3,2,l)=x(m,3,2,l)-zipd32*x(m,2,2,l)-zipd31*x(m,1,2,l)
        x(m,3,3,l)=x(m,3,3,l)-zipd32*x(m,2,3,l)-zipd31*x(m,1,3,l)
        znrm3=x(m,3,1,l)*dconjg(x(m,3,1,l))
     c       +x(m,3,2,l)*dconjg(x(m,3,2,l))
     c       +x(m,3,3,l)*dconjg(x(m,3,3,l))
*        znrm3=cdsqrt(znrm3)
        znrm3=sqrt(znrm3)
        x(m,3,1,l)=x(m,3,1,l)/znrm3
        x(m,3,2,l)=x(m,3,2,l)/znrm3
        x(m,3,3,l)=x(m,3,3,l)/znrm3
      enddo
      enddo

      do l=1,nd
      do m=1,nvect
        znrm1=y(m,1,1,l)*dconjg(y(m,1,1,l))
     c       +y(m,1,2,l)*dconjg(y(m,1,2,l))
     c       +y(m,1,3,l)*dconjg(y(m,1,3,l))
*        znrm1=cdsqrt(znrm1)
        znrm1=sqrt(znrm1)
        y(m,1,1,l)=y(m,1,1,l)/znrm1
        y(m,1,2,l)=y(m,1,2,l)/znrm1
        y(m,1,3,l)=y(m,1,3,l)/znrm1
        zipd21=y(m,2,1,l)*dconjg(y(m,1,1,l))
     c        +y(m,2,2,l)*dconjg(y(m,1,2,l))
     c        +y(m,2,3,l)*dconjg(y(m,1,3,l))
        y(m,2,1,l)=y(m,2,1,l)-zipd21*y(m,1,1,l)
        y(m,2,2,l)=y(m,2,2,l)-zipd21*y(m,1,2,l)
        y(m,2,3,l)=y(m,2,3,l)-zipd21*y(m,1,3,l)
        znrm2=y(m,2,1,l)*dconjg(y(m,2,1,l))
     c       +y(m,2,2,l)*dconjg(y(m,2,2,l))
     c       +y(m,2,3,l)*dconjg(y(m,2,3,l))
*        znrm2=cdsqrt(znrm2)
        znrm2=sqrt(znrm2)
        y(m,2,1,l)=y(m,2,1,l)/znrm2
        y(m,2,2,l)=y(m,2,2,l)/znrm2
        y(m,2,3,l)=y(m,2,3,l)/znrm2
        zipd31=y(m,3,1,l)*dconjg(y(m,1,1,l))
     c        +y(m,3,2,l)*dconjg(y(m,1,2,l))
     c        +y(m,3,3,l)*dconjg(y(m,1,3,l))
        zipd32=y(m,3,1,l)*dconjg(y(m,2,1,l))
     c        +y(m,3,2,l)*dconjg(y(m,2,2,l))
     c        +y(m,3,3,l)*dconjg(y(m,2,3,l))
        y(m,3,1,l)=y(m,3,1,l)-zipd32*y(m,2,1,l)-zipd31*y(m,1,1,l)
        y(m,3,2,l)=y(m,3,2,l)-zipd32*y(m,2,2,l)-zipd31*y(m,1,2,l)
        y(m,3,3,l)=y(m,3,3,l)-zipd32*y(m,2,3,l)-zipd31*y(m,1,3,l)
        znrm3=y(m,3,1,l)*dconjg(y(m,3,1,l))
     c       +y(m,3,2,l)*dconjg(y(m,3,2,l))
     c       +y(m,3,3,l)*dconjg(y(m,3,3,l))
*        znrm3=cdsqrt(znrm3)
        znrm3=sqrt(znrm3)
        y(m,3,1,l)=y(m,3,1,l)/znrm3
        y(m,3,2,l)=y(m,3,2,l)/znrm3
        y(m,3,3,l)=y(m,3,3,l)/znrm3
      enddo
      enddo

c------------------------ check ----------------
      rmax=0.e0
      do l=1,nd
      do i=1,3
      do m=1,nvect
        rmax=rmax+dreal(x(m,i,i,l)*dconjg(x(m,i,i,l))
     c                 +y(m,i,i,l)*dconjg(y(m,i,i,l)))
      enddo
      enddo
      enddo
      rmax=rmax/(8.d0*nsite) - 0.5d0

      rzz=0.e0
      do m=1,nvect
        wg(m,1,1)=(0.e0,0.e0)
        wg(m,2,2)=(0.e0,0.e0)
        wg(m,3,3)=(0.e0,0.e0)
      enddo

      do l=1,nd
      do m=1,nvect
        m0=mo2(m,l)
        wg(m,1,1)=wg(m,1,1)-x(m ,1,1,l)*dconjg(x(m ,2,1,l))
     c                     +x(m ,1,2,l)*dconjg(x(m ,2,2,l))
     c                     -y(m0,1,2,l)*dconjg(y(m0,1,1,l))
     c                     +y(m0,2,2,l)*dconjg(y(m0,2,1,l))
        wg(m,2,2)=wg(m,2,2)-x(m ,1,1,l)*dconjg(x(m ,3,1,l))
     c                     +x(m ,1,3,l)*dconjg(x(m ,3,3,l))
     c                     -y(m0,1,3,l)*dconjg(y(m0,1,1,l))
     c                     +y(m0,3,3,l)*dconjg(y(m0,3,1,l))
        wg(m,3,3)=wg(m,3,3)-x(m ,2,2,l)*dconjg(x(m ,3,2,l))
     c                     +x(m ,2,3,l)*dconjg(x(m ,3,3,l))
     c                     -y(m0,2,3,l)*dconjg(y(m0,2,2,l))
     c                     +y(m0,3,3,l)*dconjg(y(m0,3,2,l))
      enddo
      enddo

      do i=1,3
      do m=1,nvect
*        rzz=rzz+cdabs(wg(m,i,i))**2
        rzz=rzz+abs(wg(m,i,i))**2
      enddo
      enddo

      do m=1,nvect
        wg(m,1,1)=(0.e0,0.e0)
        wg(m,2,2)=(0.e0,0.e0)
        wg(m,3,3)=(0.e0,0.e0)
      enddo

      do l=1,nd
      do m=1,nvect
        m0=me2(m,l)
        wg(m,1,1)=wg(m,1,1)-y(m ,1,1,l)*dconjg(y(m ,2,1,l))
     c                     +y(m ,1,2,l)*dconjg(y(m ,2,2,l))
     c                     -x(m0,1,2,l)*dconjg(x(m0,1,1,l))
     c                     +x(m0,2,2,l)*dconjg(x(m0,2,1,l))
        wg(m,2,2)=wg(m,2,2)-y(m ,1,1,l)*dconjg(y(m ,3,1,l))
     c                     +y(m ,1,3,l)*dconjg(y(m ,3,3,l))
     c                     -x(m0,1,3,l)*dconjg(x(m0,1,1,l))
     c                     +x(m0,3,3,l)*dconjg(x(m0,3,1,l))
        wg(m,3,3)=wg(m,3,3)-y(m ,2,2,l)*dconjg(y(m ,3,2,l))
     c                     +y(m ,2,3,l)*dconjg(y(m ,3,3,l))
     c                     -x(m0,2,3,l)*dconjg(x(m0,2,2,l))
     c                     +x(m0,3,3,l)*dconjg(x(m0,3,2,l))
      enddo
      enddo

      do i=1,3
      do m=1,nvect
*        rzz=rzz+cdabs(wg(m,i,i))**2
        rzz=rzz+abs(wg(m,i,i))**2
      enddo
      enddo
      rzz=rzz/dfloat(nsite)
c------------------------ check ----------------

c       if(igf.le.10.or.igf.eq.50.or.mod(igf,100).eq.0)
c    c   write(*,*)igf,rzz,rmax
      write(2,*)igf,rzz,rmax,omega,ngfmax,x(1,1,1,1)

      if(rzz.le.epsrzz.or.igf.ge.ngfmax) goto 2000
 1000 continue
      igf_fin=igf-1
 2000 continue
      igf_fin=igf
      write(4,*)igf_fin,rzz,rmax

      return
      end
c*********************************************************************c
      subroutine dcuxy
c
c     * for SU(3) Lattice on VP.
c     --------------------------
c
c     * This routine converts /wrk1/ on an even-odd lattice into
c       /var3/ on an ordinary lattice for the measurement of
c       observables.
c
c     * Before this routine is called, you must once call
c         Dir.
c
c     * Include file
c         paravp3
c
c     * Input arrays
c         /wrk1/ : a configuration on an even-odd lattice.
c         /tb1/  : list vectors calculated by Dir.
c
c     * Output arrays
c         /var3/ : the converted configuration.
c
c     * programmed by S.Kitahara.
c     * modified by S.Kitahara. (98.10.18)
c----------------------------------------------------------------------

      include'paravp3'

      common /wrk1/ u(nvect,3,3,nd),v(nvect,3,3,nd)
      common /var3/ z(nsite,3,3,nd)
      common /tb1/  mioe(nsite), mioo(nsite)

      do l=1,nd
      do i=1,3
      do j=1,3
      do m=1,nsite
        if(mioe(m).eq.1) then
          z(m,i,j,l)= u(mioo(m),i,j,l)
        else
          z(m,i,j,l)= v(mioo(m),i,j,l)
        endif
      enddo
      enddo
      enddo
      enddo
      return
      end
c*********************************************************************c
      subroutine dcuxy2
c
c     * for SU(3) Lattice on VP.
c     --------------------------
c
c     * This routine converts /var2/ on an even-odd lattice into
c       /var3/ on an ordinary lattice
c       MA gauge fixing will be applied on /var3/ config.
c
c     * Before this routine is called, you must once call
c         Dir.
c
c     * Include file
c         paravp3
c
c     * Input arrays
c         /var2/ : a configuration on an even-odd lattice.
c         /tb1/  : list vectors calculated by Dir.
c
c     * Output arrays
c         /var3/ : the converted configuration.
c
c     * programmed by S.Kitahara.
c     * modified by V.Bornyakov (08.12.00)
c----------------------------------------------------------------------

      include'paravp3'

      common /var2/ u(nvect,3,3,nd),v(nvect,3,3,nd)
      common /var3/ z(nsite,3,3,nd)
      common /tb1/  mioe(nsite), mioo(nsite)

      do l=1,nd
      do i=1,3
      do j=1,3
      do m=1,nsite
        if(mioe(m).eq.1) then
          z(m,i,j,l)= u(mioo(m),i,j,l)
        else
          z(m,i,j,l)= v(mioo(m),i,j,l)
        endif
      enddo
      enddo
      enddo
      enddo
      return
      end
c*********************************************************************c
      subroutine dcuxy3
c
c     * for SU(3) Lattice on VP.
c     --------------------------
c
c     * This routine converts /var3/ on  ordinary lattice
c       into /var2/ an even-odd lattice
c       MA gauge fixing will be applied to /var2/ config.
c
c     * Before this routine is called, you must once call
c         Dir.
c
c     * Include file
c         paravp3
c
c     * Input arrays
c         /var3/ : the converted configuration.
c         /tb1/  : list vectors calculated by Dir.
c
c     * Output arrays
c         /var2/ : a configuration on an even-odd lattice.
c
c     * programmed by S.Kitahara.
c     * modified by V.Bornyakov (11.12.00)
c----------------------------------------------------------------------

      include'paravp3'

      common /var2/ u(nvect,3,3,nd),v(nvect,3,3,nd)
      common /var3/ z(nsite,3,3,nd)
      common /tb1/  mioe(nsite), mioo(nsite)

      do l=1,nd
      do i=1,3
      do j=1,3
      do m=1,nsite
        if(mioe(m).eq.1) then
          u(mioo(m),i,j,l) = z(m,i,j,l)
        else
          v(mioo(m),i,j,l) = z(m,i,j,l)
        endif
      enddo
      enddo
      enddo
      enddo
      return
      end
c*********************************************************************c
c**************************end of dcuxy3******************************c
      subroutine bqcdread(num)
c
c     Read trajectory from bqcd program output
c     ----------------------------------------
c
c     * Before this routine is called, you must once call
c         Dir.
c
c     * Include file
c         paravp3
c
c     * Input arrays
c         /U1/ : read the output produced by program QCD.pl
c     the typical call is:
c     perl  QCD.pl -c -F 4 -o data.tes4 -f  bqcd.000.1.1.00010.*.u
c    data.tes4 is the output file,
c    -c (-b)  (do not) transform big-endian to little-endian
c
c         the trajectories bqcd.000.1.1.00010.*.u are
c          generated by bqcd
c
c         /tb1/  : list vectors calculated by Dir.
c
c     * Output arrays
c         /var3/ : a configuration on the usual lattice.
c
c     * programmed by V.Bornyakov and M.I.Polikarpov (22.12.01)
c----------------------------------------------------------------------

      include'paravp3'

      character*3 num

      common /var3/ z(nsite,3,3,nd)
      common /tb1/  mioe(nsite), mioo(nsite)
c     common /wrk1/ u1(n4,n3,n2,n1,4,3,3 )
      dimension u1(n4,n3,n2,n1,4,3,3 )
c
c  read links U
      print *,'start read bqcd link fields'
      open(5,file='/home/vborn/SU3/CONFIG/bqcd.000.1.1.0'//num//'0.lat',
     *      form='unformatted')
      read(5) u1
      close(5)
      print *,'finished read bqcd link fields'
c     write(2,*) ((U1(1,1,1,1,1,i,j),i=1,3),j=1,3)
c
      do i1 = 1,n1
         do i2 = 1,n2
            do i3 = 1,n3
               do i4 = 1,n4
c        ix = n4*n3*n2*(i1-1) + n4*n3*(i2-1) + n4*(i3-1) + i4
        ns = i1+n1*(i2-1)+n1*n2*(i3-1)+n1*n2*n3*(i4-1)
                  do mu =1,4
                  do i=1,3
                  do j=1,3
        z(ns,i,j,mu) = u1(i4,i3,i2,i1,mu,i,j)
                  end do
                  end do
                  end do
               end do
            end do
          end do
       end do
      return
      end
c*********************************************************************c
c************************** end of bqcdread **************************c

c*********************************************************************c
      subroutine dir
c
c     * for SU(2) and SU(3) Lattice on VP.
c     ------------------------------------
c
c     * Compute list vectors under periodic boundary condition.
c       For the measurement of observables, the list vectors computed
c       here are on the ordinary lattice.
c
c
c     ........       nvect           ........       nsite
c     _     _     _
c     7  7  8  8  9  9               13 14 15 16 17 18
c        _     _     _
c     4  4  5  5  6  6      ===>>    7  8  9  10 11 12
c     _     _     _
c     1  1  2  2  3  3               1  2  3  4  5  6
c
c     even-odd site numbers.         ordinary site numbers.
c     Bar{n} means odd site.
c
c
c     * List vectors Mioe and Mioo are used in order to convert
c       arrays on even-odd lattices into those on ordinary lattices.
c         mioe(n) = 1   (n means odd site)
c                 = 0   (n means even site), n: ordinary site number.
c
c         mioo(n) = even-odd site number   , n: ordinary site number.
c
c     * Im and nm give the neighboring site numbers on ordinary
c       lattice.
c         im(n,mu,i) = # of (n + (i*mu))   , n: ordinary site number
c         nm(n,mu)   = # of (n - mu)       , n: ordinary site number
c
c     * Include file
c         paravp3
c
c     * Output arrays and variables
c         /tb1/
c         /sdb2/
c
c     * Work arrays
c         /sdb1/
c
c     * programmed by KEK group.
c     * modified by S.Kitahara. (96.8.2)
c----------------------------------------------------------------------

      include'paravp3'

      common /tb1  /mioe(nsite), mioo(nsite)
      common /sdb1 /ipower(nd),isize(nd),is(nd),mio(nvect),mie(nvect)
      common /sdb2 /im(nsite,nd,n0), nm(nsite,nd)

      isize(1)=n1
      isize(2)=n2
      isize(3)=n3
      isize(4)=n4
      ipower(1)= 1
      ipower(2)= n1
      ipower(3)= n1*n2
      ipower(4)= n1*n2*n3

      do 1 i4=1,n4
        is(4)= i4
        m4   = (i4-1)*ipower(4)
        mm4  = (i4-1)*n0*n2*n3
      do 1 i3=1,n3
        is(3)= i3
        m3   = m4 +(i3-1)*ipower(3)
        mm3  = mm4+(i3-1)*n0*n2
      do 1 i2=1,n2
        is(2)= i2
        m2   = m3 +(i2-1)*ipower(2)
        mm2  = mm3+(i2-1)*n0
        ioe  = mod(i2+i3+i4,2)
      do 1 i1=1,n0
        is(1)= i1
        m1   = m2 +2*i1
        mm1  = mm2+i1
        mio(mm1) = m1-ioe
        mie(mm1) = m1+ioe-1
    1 continue

      do 2 m=1,nvect
        mioe(mio(m)) = 1
        mioo(mio(m)) = m
        mioe(mie(m)) = 0
        mioo(mie(m)) = m
    2 continue

      do 3 i4=1,n4
        is(4)= i4
        m4   = (i4-1)*ipower(4)
      do 3 i3=1,n3
        is(3)= i3
        m3   = m4 +(i3-1)*ipower(3)
      do 3 i2=1,n2
        is(2)= i2
        m2   = m3 +(i2-1)*ipower(2)
      do 3 i1=1,n1
        is(1)= i1
        m1   = m2 +i1
      do 3 l=1,nd
        nm(m1,l)=m1+(mod(is(l)-2+isize(l),isize(l))+1-is(l))*ipower(l)
      do 3 k=1,n0
        im(m1,l,k) = m1+(mod(is(l)+k-1,isize(l))+1-is(l))*ipower(l)
    3 continue
      return
      end
c*********************************************************************c
      subroutine dirct
c
c     * for SU(2) and SU(3) Lattice on VP.
c     ------------------------------------
c
c     * Compute list vectors under periodic boundary condition.
c       For vector processing, the entire lattice is divided into
c       two parts: even lattice sites and odd lattice sites.
c
c     * The list vectors are mainly used by Monte and Mafix1.
c
c
c     ........       nvect
c     _     _     _
c     7  7  8  8  9  9
c        _     _     _
c     4  4  5  5  6  6
c     _     _     _
c     1  1  2  2  3  3
c
c     even-odd site numbers.
c     Bar{n} means odd site.
c
c
c     * mo1,mo2,me1 and me2 give the neighboring site numbers.
c         mo1(n,mu) = # of (n + mu)   (n: odd,  n+mu: even)
c         mo2(n,mu) = # of (n - mu)   (n: odd,  n-mu: even)
c         me1(n,mu) = # of (n + mu)   (n: even, n+mu: odd )
c         me2(n,mu) = # of (n - mu)   (n: even, n-mu: odd )
c
c     * ldir2(i,mu) (i=1,2 and 3) gives a direction perpendicular to mu.
c       ldir2(i,1),ldir2(i,2),ldir2(i,3) and ldir2(i,4) are
c       perpendicular to each other.
c
c     * Include file
c         paravp3
c
c     * Output arrays and variables
c         /tabl/
c
c     * Work arrays
c         /vard/
c
c     * programmed by KEK group.
c     * modified by S.Kitahara. (96.8.2)
c----------------------------------------------------------------------

      include'paravp3'

      common /tabl /mo1(nvect,nd), mo2(nvect,nd),
     &              me1(nvect,nd), me2(nvect,nd),
     &              ldir2(nd-1,nd)
      common /vard /mcycle(0:n0+1), mup(50,2:nd), mdn(50,2:nd),
     &              ipower(2:nd),   is(2:nd),     isize(nd)

      do i =1,3
      do mu=1,nd
        ldir2(i,mu) = mod(i+mu-1,nd) + 1
      enddo
      enddo

      isize(1)   =  n1
      isize(2)   =  n2
      isize(3)   =  n3
      isize(4)   =  n4

      ipower(2)  =  n0
      ipower(3)  =  n0*n2
      ipower(4)  =  n0*n2*n3

      do 1 m  = 0,n0+1
    1    mcycle(m) = mod(m-1+n0,n0)+1

      do 2 i = 2,nd
      do 2 m = 1,50
         mup(m,i)  = mod(m,isize(i))+1
    2    mdn(m,i)  = mod(m-2+isize(i),isize(i))+1
      do 3 i4  = 1,n4
         is(4) = i4
         m4    = (i4-1)*ipower(4)
      do 3 i3  = 1,n3
         is(3) = i3
         m3    = m4+(i3-1)*ipower(3)
      do 3 i2  = 1,n2
         is(2) = i2
         m2    = m3+(i2-1)*ipower(2)
         ioe   = mod(i2+i3+i4,2)
      do 3 i1  = 1,n0
         m1    = m2+i1
         mo1(m1,1) = m1 + mcycle(i1+1-ioe)-i1
         mo2(m1,1) = m1 + mcycle(i1  -ioe)-i1
         me1(m1,1) = m1 + mcycle(i1  +ioe)-i1
         me2(m1,1) = m1 + mcycle(i1-1+ioe)-i1
       do 11 i = 2,nd
           mo1(m1,i) = m1 + (mup(is(i),i)-is(i))*ipower(i)
           me1(m1,i) = mo1(m1,i)
           mo2(m1,i) = m1 + (mdn(is(i),i)-is(i))*ipower(i)
           me2(m1,i) = mo2(m1,i)
   11  continue
    3  continue
       return
       end
c*********************************************************************c
      subroutine inits(ints,pi2)
c
c     * for SU(3) Lattice on VP.
c     --------------------------
c
c     * This routine returns an initial configuration.
c
c     * This routine calls
c         Rndpr3.
c
c     * Include file
c         paravp3
c
c     * Input parameters and constants
c         ints   : 1 = cold start. 0 = hot start.
c         pi2    : 2 times pi (6.28318...).
c
c     * Output arrays and variables
c         /var2/ : generated initial configuration.
c
c     * Work arrays
c         /ran1/ : random numbers generated by Rndpr3.
c
c     * programmed by S.Kitahara.
c     * modified by S.Kitahara. (98.10.14)
c----------------------------------------------------------------------

      include'paravp3'

      common /var2/ u(nvect,3,3,nd),v(nvect,3,3,nd)
      common /ran1/ rm1(nvect),rm2(nvect),rm3(nvect),rm4(nvect)

      if(ints.eq.0)then
        do i = 1,3
        do j = 1,3
        do mu = 1,nd
        do m = 1,nvect
          u(m,i,j,mu) = (0.e0,0.e0)
          v(m,i,j,mu) = (0.e0,0.e0)
        enddo
        enddo
        enddo
        enddo

        do i = 1,3
        do mu = 1,nd
        do m = 1,nvect
          u(m,i,i,mu) = (1.e0,0.e0)
          v(m,i,i,mu) = (1.e0,0.e0)
        enddo
        enddo
        enddo

      elseif(ints.eq.1)then

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y11=dcmplx( a0, a3)
            y12=dcmplx( a2, a1)
            y21=dcmplx(-a2, a1)
            y22=dcmplx( a0,-a3)
            u(m,1,1,mu)=y11
            u(m,1,2,mu)=y12
            u(m,2,1,mu)=y21
            u(m,2,2,mu)=y22
          enddo
        enddo

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y22=dcmplx( a0, a3)
            y23=dcmplx( a2, a1)
            y32=dcmplx(-a2, a1)
            y33=dcmplx( a0,-a3)
            u12=u(m,1,2,mu)
            u22=u(m,2,2,mu)
            u(m,1,2,mu)=u12*y22
            u(m,1,3,mu)=u12*y23
            u(m,2,2,mu)=u22*y22
            u(m,2,3,mu)=u22*y23
            u(m,3,2,mu)=y32
            u(m,3,3,mu)=y33
          enddo
        enddo

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y11=dcmplx( a0, a3)
            y13=dcmplx( a2, a1)
            y31=dcmplx(-a2, a1)
            y33=dcmplx( a0,-a3)
            u11=u(m,1,1,mu)
            u13=u(m,1,3,mu)
            u21=u(m,2,1,mu)
            u23=u(m,2,3,mu)
            u33=u(m,3,3,mu)
            u(m,1,1,mu)=u11*y11+u13*y31
            u(m,1,3,mu)=u11*y13+u13*y33
            u(m,2,1,mu)=u21*y11+u23*y31
            u(m,2,3,mu)=u21*y13+u23*y33
            u(m,3,1,mu)=u33*y31
            u(m,3,3,mu)=u33*y33
          enddo
        enddo

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y11=dcmplx( a0, a3)
            y12=dcmplx( a2, a1)
            y21=dcmplx(-a2, a1)
            y22=dcmplx( a0,-a3)
            v(m,1,1,mu)=y11
            v(m,1,2,mu)=y12
            v(m,2,1,mu)=y21
            v(m,2,2,mu)=y22
          enddo
        enddo

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y22=dcmplx( a0, a3)
            y23=dcmplx( a2, a1)
            y32=dcmplx(-a2, a1)
            y33=dcmplx( a0,-a3)
            v12=v(m,1,2,mu)
            v22=v(m,2,2,mu)
            v(m,1,2,mu)=v12*y22
            v(m,1,3,mu)=v12*y23
            v(m,2,2,mu)=v22*y22
            v(m,2,3,mu)=v22*y23
            v(m,3,2,mu)=y32
            v(m,3,3,mu)=y33
          enddo
        enddo

        do mu=1,nd
          call rndpr3(4)
          do m=1,nvect
            a0  =rm1(m)
            rad = 1.e0 - a0**2
            a3  = dsqrt(rad)*(2.e0*rm2(m)-1.e0)
            rad2 = dsqrt(dabs(rad-a3**2))
            theta= pi2*rm3(m)
            a1 = rad2*dcos(theta)
            a2 = rad2*dsin(theta)
            y11=dcmplx( a0, a3)
            y13=dcmplx( a2, a1)
            y31=dcmplx(-a2, a1)
            y33=dcmplx( a0,-a3)
            v11=v(m,1,1,mu)
            v13=v(m,1,3,mu)
            v21=v(m,2,1,mu)
            v23=v(m,2,3,mu)
            v33=v(m,3,3,mu)
            v(m,1,1,mu)=v11*y11+v13*y31
            v(m,1,3,mu)=v11*y13+v13*y33
            v(m,2,1,mu)=v21*y11+v23*y31
            v(m,2,3,mu)=v21*y13+v23*y33
            v(m,3,1,mu)=v33*y31
            v(m,3,3,mu)=v33*y33
          enddo
        enddo
      endif

      return
      end
c*********************************************************************c
      subroutine rndpr3(nrn)
c
c     * for SU(2) and SU(3) Lattice on VP.
c     ------------------------------------
c
c     * random number generator ( m-series,p=250,q=103 )
c
c     * Before this routine is called, you must once call
c         Seed.
c
c     * Include file
c         paravp3
c
c     * Input arrays and parameters
c         nrn    : must be 2 or 4.
c                  Rm1 and Rm2 are changed(nrn = 2).
c                  Rm1, Rm2, Rm3 and Rm4 are changed(nrn = 4).
c         /rndk/ : Initial value of K1 and K2 are 0 and 147
c                  respectively which you must give elsewhere.
c         /rndlp/: The array, Irand is computed by Seed.
c
c     * Output arrays
c         /ran1/
c
c     * Work arrays
c         /randd/
c
c     * programmed by J.Makino
c     * modified   by S.Hioki
c     * modified by S.Kitahara. (96.8.2)
c----------------------------------------------------------------------

      include'paravp3'

      parameter( nkei=256 )

      common /rndlp/irand(nkei,250)
      common /randd /rn(nsite*2+nkei)
      common /rndk /k1,k2
      common /ran1 /rm1(nvect),rm2(nvect),rm3(nvect),rm4(nvect)


c     ic = 2 ** 31 - 1
c     rc = 1.0/ic
      ic = 2 ** 30
      rc = 1.0/(ic*2.-1.)

c --  main production

      mkei=(nvect*nrn-1)/nkei
      do 10  j0 = 0 , mkei
        j1=nkei*j0
        k1=mod(k1,250)+1
        k2=mod(k2,250)+1
      do 20 j = 1 , nkei
        irr           = ieor( irand(j,k1),irand(j,k2) )
        rn(j+j1)      = dble( irr )*rc
        irand(j,k1)= irr
   20 continue
   10 continue

      if(nrn.eq.2) then
        do m=1,nvect
          rm1(m)=rn(m)
          rm2(m)=rn(m+nvect)
        enddo
      else
        do m=1,nvect
          rm1(m)=rn(m)
          rm2(m)=rn(m+nvect)
          rm3(m)=rn(m+nvect2)
          rm4(m)=rn(m+nvect3)
        enddo
      endif

      return
      end
c*********************************************************************c
      subroutine seed
c
c     * for SU(2) and SU(3) Lattice on VP.
c     ------------------------------------
c
c     * random number seed preparation.
c
c     * Output arrays
c         /rndlp/
c
c     * programmed by J.Makino
c----------------------------------------------------------------------

      parameter( nkei=256 )
      common /rndlp/irand(nkei,250)
      dimension irand0(500),l(0:499)
      data ix/1774315169/
      do 110 i1=1,250
      do 110 i2=1,nkei
        irand(i2,i1)=0
  110 continue
      do 10 i=1,250
        ix=iand(ix*48828125,2147483647)
        irand0(i)=ix
   10 continue
      do 20 i=251,500
        irand0(i)=ieor(irand0(i-250),irand0(i-103))
   20 continue
      do 900 k=1,nkei
        do 100 i=0,499
          l(i)=0
  100   continue
        if(k.le.249) then
          l(k)=1
        else
          l(k-250)=1
          l(k-103)=1
        endif
        do 500 j=242,1,-1
          do 200 i=249,0,-1
            l(2*i+1)=0
            l(2*i)=l(i)
  200     continue
          do 300 i=498,250,-1
            l(i-250)=ieor(l(i-250),l(i))
            l(i-103)=ieor(l(i-103),l(i))
  300     continue
  500   continue
        do 700 j=1,250
          do 600 i=0,249
            irand(k,j)=ieor(irand(k,j),l(i)*irand0(i+j))
  600     continue
  700   continue
  900 continue
      return
      end
