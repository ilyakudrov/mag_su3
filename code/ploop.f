      SUBROUTINE PLOOP(icolor) 
C              VERSION 01.2001
C TO COMPUTE ABELIAN POLYAKOV LOOP IN ALL DIRECTIONS
C***********************************************************************
      include 'paravp3'  
      parameter(lengs=n1,lengt=n4)
c     include 'su2mon.cmn'

      PARAMETER(NP1=4,NP2=lengs*NP1,NP3=lengs*NP2,NP4=lengs*NP3)
      PARAMETER(LP1=NP1*(1-lengs),LP2=NP2*(1-lengs),LP3=NP3*(1-lengs),
     *          LP4=NP4*(1-lengt))
      PARAMETER(NST3=LENGS**2*LENGT)
      PARAMETER(NST4=LENGS**3)

      common /ablink/ ab(nsite,nd,3)         ! ab. link
      common /sm_abl/ PH(nlink)              ! ab. link smeared


c     complex pll1,...
c     common /ploops/ pll1(lengs,lengs,lengt),....
      complex  STR(4),plo

      do id=1,4
      do is=1,nsite
        il=(is-1)*4+id
        PH(il)=ab(is,id,icolor)
      enddo
      enddo
      do i=1,4
       STR(i)=cmplx(0.,0.)
      enddo   
C************************************************************
C PLOOP IN THE 4th DIRECTION
C************************************************************
      DO 1 I3=1,LENGS
      DO 1 I2=1,LENGS
      DO 1 I1=1,LENGS

       S0=0.
       NS=(I1-1)*NP1+(I2-1)*NP2+(I3-1)*NP3

      DO 4 I4=1,LENGT

       IL=NS+(I4-1)*NP4+4
       S0=S0+PH(IL)

 4    CONTINUE

      plo=cmplx(cos(S0),sin(S0))
c     pll(i1,i2,i3)=plo
      STR(4)=STR(4)+plo

 1    CONTINUE

      STR(4)=STR(4)/NST4
C************************************************************
C PLOOP IN THE 3rd DIRECTION
C************************************************************
      DO 11 I4=1,LENGT
      DO 11 I2=1,LENGS
      DO 11 I1=1,LENGS

       S0=0.
       NS=(I1-1)*NP1+(I2-1)*NP2+(I4-1)*NP4

      DO 13 I3=1,LENGS

       IL=NS+(I3-1)*NP3+3
       S0=S0+PH(IL)

 13    CONTINUE

      plo=cmplx(cos(S0),sin(S0))
c     pll(i1,i2,i4)=plo
      STR(3)=STR(3)+plo

 11    CONTINUE
      STR(3)=STR(3)/NST3
C************************************************************
C PLOOP IN THE 2nd DIRECTION
C************************************************************
      DO 21 I4=1,LENGT
      DO 21 I3=1,LENGS
      DO 21 I1=1,LENGS
       S0=0.
       NS=(I1-1)*NP1+(I3-1)*NP3+(I4-1)*NP4
      DO 22 I2=1,LENGS
       IL=NS+(I2-1)*NP2+2
       S0=S0+PH(IL)
 22    CONTINUE

       plo=cmplx(cos(S0),sin(S0))
c      pll(i1,i2,i4)=plo
       STR(2)=STR(2)+plo
 21   CONTINUE
      STR(2)=STR(2)/NST3
C************************************************************
C PLOOP IN THE 1st DIRECTION
C************************************************************
      DO 32 I4=1,LENGT
      DO 32 I3=1,LENGS
      DO 32 I2=1,LENGS
       S0=0.
       NS=(I2-1)*NP2+(I3-1)*NP3+(I4-1)*NP4
      DO 31 I1=1,LENGS
       IL=NS+(I1-1)*NP1+1
       S0=S0+PH(IL)
 31    CONTINUE

       plo=cmplx(cos(S0),sin(S0))
c      pll(i1,i2,i4)=plo
       STR(1)=STR(1)+plo
 32   CONTINUE
      STR(1)=STR(1)/NST3
      
      write(11,*) STR(4)
      write(11,*) (STR(1)+STR(2)+STR(3))/3.
      RETURN
      END
C**********************************************************************
