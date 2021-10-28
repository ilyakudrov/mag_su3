C***********************************************************************
      SUBROUTINE SMEAR_ABEL(EP,icolor)
C**********************************************************************
C                   VERSION: 25.04.1993
C THE SMEARING IS DUE TO REPLACE U(X,MU) BY
C                                +          +
C EP*U(X,MU)+SUM{U(X,NU)U(X+NU,MU)U(X+MU,NU)+
C               +U(X-NU,NU)U(X-NU,MU)U(X+MU-NU,NU)}
C
C**********************************************************************

      include 'paravp3'   
c     include 'su2mon.cmn'
c     PARAMETER(EP=0.)
      parameter(nsmr=10)
      parameter(lengs=n1,lengt=n4)
      PARAMETER(MP1=4,MP2=lengs*MP1,MP3=lengs*MP2,MP4=lengs*MP3)    
      PARAMETER(KP1=MP1*(1-lengs),KP2=MP2*(1-lengs),KP3=MP3*(1-lengs),
     *          KP4=MP4*(1-lengt))

      common /ablink/ ab(nsite,nd,3)         ! ab. link
      common /sm_abl/ PH(nlink)              ! ab. link smeared
      DIMENSION NP(4),NN(4)

      do id=1,4     
      do is=1,nsite     
        il=(is-1)*4+id
        PH(il)=ab(is,id,icolor)
      enddo
      enddo

      do ism=1,nsmr

      NS0=-4
C***********************************************************************
C       THE CYCLES OVER  ALL SITES ARE HERE:
C***********************************************************************
      DO 504 I4=1,LENGT
      NP(4)=MP4
      NN(4)=MP4
      IF (I4.EQ.LENGT) NP(4)=KP4
      IF (I4.EQ.1)     NN(4)=KP4

      DO 503 I3=1,LENGS
      NP(3)=MP3
      NN(3)=MP3
      IF (I3.EQ.LENGS) NP(3)=KP3
      IF(I3 .EQ. 1)    NN(3)=KP3

      DO 502 I2=1,LENGS
      NP(2)=MP2
      NN(2)=MP2
      IF (I2.EQ.LENGS) NP(2)=KP2
      IF (I2.EQ.1)     NN(2)=KP2

      DO 501 I1=1,LENGS
      NP(1)=MP1
      NN(1)=MP1
      IF (I1.EQ.LENGS) NP(1)=KP1
      IF (I1.EQ.1)     NN(1)=KP1

C***********************************************************************
C     HERE IS THE ENUMERATION OF SITES:
C***********************************************************************
      NS0=NS0+4
C***********************************************************************
C   SELECT LINK IN MUE-DIRECTION:
C***********************************************************************
      DO 400 MUE=1,3

      P0=0.
      P3=0.

      IL=NS0+MUE
      A0=PH(IL)
C***********************************************************************
C   COMPUTE THE STAPLE IN (MUE,NUE)-PLANE:
C***********************************************************************
      DO 20 NUE=1,3
      IF (NUE.EQ.MUE) GO TO 20

      L3=NS0+NUE
      L2=IL+NP(NUE)
      L1=L3+NP(MUE)
C***********************************************************************
C      "MULTIPLICATION" OF THREE LINKS IN "POSITIVE" DIRECTION:
C***********************************************************************
                                   !         L2
      B0=PH(L1)                    !       *-----*
      C0=PH(L2)                    !       |     |
      D0=PH(L3)                    !     L3|     |L1
                                   !       |     |
      B0=B0-C0-D0                  !     NS+-----+
      P0=P0+cos(B0)                !       |  IL |
      P3=P3-sin(B0)                !     L3|     |L1
                                   !       |     |
      L3=NS0-NN(NUE)+NUE            !       *-----*
      L2=IL-NN(NUE)                !         L2
      L1=L3+NP(MUE)

C***********************************************************************
C      MULTIPLICATION OF THREE LINKS IN "NEGATIVE" DIRECTION:
C***********************************************************************
      B0=PH(L1)  
      C0=PH(L2) 
      D0=PH(L3) 

      B0=-B0-C0+D0
      P0=P0+cos(B0) 
      P3=P3-sin(B0)  

20    CONTINUE

       IF(EP.EQ.0.) go to 100
       P0=EP*COS(A0)+P0
       P3=EP*SIN(A0)+P3

  100 continue

      FIS=ATAN2(P3,P0) 
      PH(IL)=FIS

 400  CONTINUE
 501  CONTINUE
 502  CONTINUE
 503  CONTINUE
 504  CONTINUE
      
      enddo

      RETURN
      END
C**********************************************************************
C**********************************************************************
      SUBROUTINE WILT_ABEL_A

C**********************************************************************
C         LWLMAX-maximal size of the Wilson loop,
c!! loops up to LWLMAX*LWLMAX are computed (correct for lengt>lengs) !!
C**********************************************************************
      include 'paravp3'  
c     include 'su2mon.cmn'

      parameter(lengs=n1,lengt=n4)
      PARAMETER(lwlmax=lengs/2)
      PARAMETER(dnorm=2./(6.*nsite))
      PARAMETER(MP1=4,MP2=lengs*MP1,MP3=lengs*MP2,MP4=lengs*MP3)    
      parameter(nlip=2*lengs*lengt)

c     real*8 wltab(lwlmaxt,lwlmaxt)
      real*8 wltab(lwlmax,lwlmax)

      common /sm_abl/ PH(nlink)              ! ab. link smeared

      DIMENSION NP(4),KN(4),NPP(2),SL(4)
c     DIMENSION YA(LWLMAXT,NLIP)
      DIMENSION YA(lwlmax,nlip)
      DIMENSION LST(4)
      DIMENSION itab(3)

      data itab/2,3,1/
      INTEGER TAU,RHO

        NP(1)=MP1
        NP(2)=MP2
        NP(3)=MP3
        NP(4)=MP4
	DO 21 I=1,2
 21	NPP(I)=NP(I)*2/4
C***********************************************************************
C   INITIALISATION OF WILSON LOOPS:
C***********************************************************************
c     DO 2 I=1,LWLMAXT
c     DO 2 J=1,LWLMAXT
      DO  I=1,LWLMAX
      DO  J=1,LWLMAX
       WLTAB(I,J)=0.D0
      ENDDO
      ENDDO

C***********************************************************************
C   SELECT THE LOOP'S PLANE (MUE,NUE):
C***********************************************************************
      DO 4 MUE=1,3
      NUE=4

      TAU=itab(MUE)
      RHO=6-MUE-TAU

C   THE CYCLES IN ORTO-PLANE DIRECTION:

      DO 3 IT=1,LENGS
      DO 3 IR=1,LENGS

        NSS=(IT-1)*NP(TAU)+(IR-1)*NP(RHO)

C   THE CYCLES OVER SELECTED PLANE:

      DO 1 IM=1,LENGS
      DO 1 IN=1,LENGT

C   HERE IS THE ENUMERATION OF SITES:

        NS0=(IM-1)*NP(MUE)+(IN-1)*NP(NUE) + NSS

        NLA=NS0+NUE
        MLA=NS0+MUE

        NPA=(IM-1)*NPP(1)+(IN-1)*NPP(2)
C***********************************************************************
C   FILLING PLANE ARRAY WITH LINKS OF DIFFERENT SIZE:
C***********************************************************************

C   1-LINKS:                      !     |
                                  !     |
        YA(1,NPA+1)=PH(MLA)       !     |
        YA(1,NPA+2)=PH(NLA)       !     +------*
                                  !    NPA
C   2,...LWLMAXT-LINKS:

c     DO 9 J=2,LWLMAXT
      DO 9 J=2,LWLMAX

      J1=J-1
      NO=J1
      MO=J1
      IF(IM+J1.GT.LENGS) MO=J1-LENGS
      IF(IN+J1.GT.LENGT) NO=J1-LENGT
        ILMA=MLA+MO*NP(MUE)
        ILNA=NLA+NO*NP(NUE)

C   ABELIAN LINKS:    

        A0=YA(J1,NPA+1)
        B0=PH(ILMA)
        YA(J,NPA+1)=A0+B0

        A0=YA(J1,NPA+2)
        B0=PH(ILNA)
        YA(J,NPA+2)=A0+B0

 9    CONTINUE
 1    CONTINUE
C***********************************************************************
C     CONSTRUCT WILSON LOOPS WL(L1,L2) IN (MUE,NUE)-PLANE
C***********************************************************************
      DO 6 IM=1,LENGS
      DO 6 IN=1,LENGT

C   HERE IS THE ENUMERATION OF SITES IN PLANE:

      KN(1)=(IM-1)*NPP(1)+(IN-1)*NPP(2)+1
      KN(4)=KN(1)+1

c     DO 6 LM=1,LWLMAXT
      DO 6 LM=1,LWLMAX
c     DO 6 LN=1,LWLMAXT
      DO 6 LN=1,LWLMAX

      MO=LM
      NO=LN
      IF(IM+LM.GT.LENGS) MO=LM-LENGS
      IF(IN+LN.GT.LENGT) NO=LN-LENGT
      KN(3)=KN(1)+NO*NPP(2)
      KN(2)=KN(4)+MO*NPP(1)
      
      DO  J=1,4
        IF(MOD(J,2).EQ.1)  then
          L=LM
        else
          L=LN
        endif

        SL(J)=YA(L,KN(J))
      ENDDO

      WLTAB(LM,LN)=WLTAB(LM,LN)+cos(SL(1)+SL(2)-SL(3)-SL(4))

 6    CONTINUE
 3    CONTINUE
 44   CONTINUE
 4    CONTINUE

C***********************************************************************
C     NORMALISATION OF LOOPS:
C***********************************************************************
c     DO 7 LM=1,LWLMAXT
      DO 7 LM=1,LWLMAX
c     DO 7 LN=1,LWLMAXT
      DO 7 LN=1,LWLMAX
c     WLTAB(LM,LN)=2.*WLTAB(LM,LN)/NPLAQ
      WLTAB(LM,LN)=WLTAB(LM,LN)*dnorm
 7    CONTINUE

c       write(6,*) 'WIL_T_AB'
C       write(6,*)  WLTAB
       write(13,*)  WLTAB
c      do  lm=1,lwlmaxt
c       WRITE(6,68) LM,(WLT(LM,LN),LN=1,2),(WLTAB(LM,LN),LN=1,2)
c      enddo

  68  FORMAT(/2X,'LM=',I3/2X,2F9.4/2X,2F9.4/)
      RETURN
      END
C***********************************************************
