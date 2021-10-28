C******************************************************************************
C*****************************************************************************
      double precision function rnd(dmy)

c     reial*8 rn,rndin,rndout,dmy
c     integer*4 il,ih,ml,mh,ial,iah,ix,iyl,ixl,ixh,kl,kh

      data mh,ml/8188,5197/,iah,ial/1731,1419/,
     *ih,il/67108864,8192/, ixh,ixl/0,0/,rn/1.4901161d-08/

      iyl=ml*ixl+ial
      ixh=ml*ixh+mh*ixl+iah
      ixl=iyl
      ixh=mod(ixh,il)*il
      ix=mod(ixl+ixh,ih)
      ixl=mod(ix,il)
      ixh=(ix-ixl)/il
      rnd=ix*rn

      return
C*****************************************************************************
C*****************************************************************************
      entry rndin(kh,kl)

      kh=iabs(kh)
      kh=mod(kh,il)

      kl=iabs(kl)
      kl=mod(kl,il)

      ixh=kh
      ixl=kl

      rndin=0.d0
      return
C***********************************************************************
C***********************************************************************
      entry rndout(kh,kl)

      kh=ixh
      kl=ixl

      rndout=0.d0
      return
      end
C**********************************************************************
