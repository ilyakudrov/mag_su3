c*********************************************************************c
c*********************************************************************c
      subroutine plaq
c
c     * for SU(3) Lattice on VP.
c     --------------------------
c
c     * Calculate SU(3) plaquette
c
c     * Before this routine is called, you must once call
c         Dir.
c
c     * Include file
c         paravp3
c
c     * Input arrays
c         /var3/ : the SU(3) configuration generated by Monte
c                  and converted by Dcuxy.
c         /sdb2/ : list vectors calculated by Dir.
c
c     * Output variables
c       plq - average plaquette
c     * Work arrays
c         /wrk5/
c
c     * programmed by S.Kitahara
c     * modified by S.Kitahara. (98.10.15)
c----------------------------------------------------------------------

      include'paravp3'
      parameter (nplaq=nsite*6)

      common /sdb2 / im(nsite,nd,n1),nm(nsite,nd)
      common /var3 / z(nsite,3,3,nd)
      common /wrk1 /zl(nsite,3,3),zk(nsite,3,3)
     c             ,zw(nsite,3,3),zzw(nsite,3,3)

      dimension fwt(n0,n0)

      plq=0.
      plqt=0.
      plqs=0.

      do k=1,nd
c     do k=nd,nd          ! only time-like WL
      do l=1,k-1

        do ii=1,3
        do jj=1,3
        do m =1,nsite
          zl(m,ii,jj)=z(m,ii,jj,l)
          zk(m,ii,jj)=z(m,ii,jj,k)
        enddo
        enddo
        enddo

          do ii=1,3
          do jj=1,3
          do m =1,nsite
            mj=im(m,k,1)
            zw(m,ii,jj) =zk(m,ii,1)*zl(mj,1,jj)
     c                  +zk(m,ii,2)*zl(mj,2,jj)
     c                  +zk(m,ii,3)*zl(mj,3,jj)
            mi=im(m,l,1)
            zzw(m,ii,jj)=zl(m,ii,1)*zk(mi,1,jj)
     c                  +zl(m,ii,2)*zk(mi,2,jj)
     c                  +zl(m,ii,3)*zk(mi,3,jj)
          enddo
          enddo
          enddo
          ztr=(0.e0,0.e0)
          do ii=1,3
          do m =1,nsite
            ztr   =ztr
     c            +zzw(m,ii,1)*dconjg(zw(m,ii,1))
     c            +zzw(m,ii,2)*dconjg(zw(m,ii,2))
     c            +zzw(m,ii,3)*dconjg(zw(m,ii,3))
          enddo
          enddo
c         fwt(i,j)=fwt(i,j)+dble(ztr)
          if(k.eq.nd) then
            plqt=plqt + dble(ztr)
          else
            plqs=plqs + dble(ztr)
          endif
          plq    =plq    +dble(ztr)

      enddo
      enddo

        plq=plq/(nplaq*3)        ! 1/3 for trace
        plqt=2.*plqt/(nplaq*3.)        ! 1/3 for trace
        plqs=2.*plqs/(nplaq*3.)        ! 1/3 for trace

        print *, plq,plqt,plqs
        write(10,*) plq,plqt,plqs
      return
      end
c*********************************************************************c