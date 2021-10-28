c*********************************************************************c
      subroutine save_ab(num,irnd) 
c 
c     * for SU(3) Lattice on VP.
c     --------------------------
c to store abel. config into file 
c     * A common block, /abli/ is the abelian  configuration.
c
c     * Include file
c         paravp3
c
c     * Input arrays and parameters
c     *   /abli/ is the abelian  configuration.
c
c     * by V.Bornyakov 01.2002 
c----------------------------------------------------------------------
      include'paravp3'
c     character*3 num
      character*5 num
      character*2 nr(10),nrn
      data nr/'01','02','03','04','05','06','07','08','09','10'/
      common/abli/  ab(nsite,3,nd)         ! ab. link 

      nrn=nr(irnd)
      open(11,status='unknown',form='unformatted',
c    *     file='su3_ab.b5900.'//num//'.lat'//nrn)
     *     file='bqcd_ab.b5p25.k13605'//num//'.lat')
      do l=1,nd
      do i=1,2
      do m=1,nsite
        write(11) ab(m,i,l) 
      enddo
      enddo
      enddo
      close(11)
 2000 return
      end
