C**********************************************************************
C**********************************************************************
      SUBROUTINE MLNEW5(icolor)
C                                                   29.07.91
c                                         updated   19.12.00
C**********************************************************************
c  this is modified program mloops tuned to consider the loop with 
c  self-intersections as one cluster
c  it computes correlation function lcrf(r) also
C******************************************************************** 
C  MONOPOLE LOOP STUDY 
C  WITH FOLLOWING OBSERVABLES:
C
C  MAXIS=MAXIMAL NUMBER OF INTERSECTIONS FOR ONE CLUSTER
C  NL=NUMBER OF LOOPS
C  NPL=PERIMETER OF A LOOP (LOOP LENGTH)
C  NT=LOOP PERIMETER SUMMED UP
C  STAT=STATICITY OF A LOOP
C  TOT=LOOP STATICITY SUMMED UP
C  IDI( )=LOOP PERIMETER DISTRIBUTION
C  ICO( )=LOOP DIAMETER CORRELATION
C  NCRO( )=BOUND. CROSS. NUMBER FOR CURRENT LOOP
C  NCROS( )= SAME FOR THE HOLE CONFIG.
C  NQSUM( )= SUM OF MONOP. CHARGES OVER CURRENT LOOP
C  NMONO( )= NUMBER OF MON. AND ANTIMON. FOR THE WHOLE CONFIG.
C  LENGS,LENGT ASSUMED TO BE EVEN
C  PROCEDURE USES TOPOLOGICALLY DEF. MAGN. CHARGE DENSITY NS(4,NLI4/16 )
C
c     include 'param.txt'
c     include 'su2mon.cmn'
      include'paravp3'

      parameter(lengs=n1)
      parameter(lengt=n4)

      PARAMETER(MAXD=3*(LENGS/2)**2+(LENGT/2)**2)
C LPMAX-MAXIMAL LOOP LENGTH FOR WHICH GCR CAN BE CALCULATED
c these parameters must be adjusted depending on LENGS
c and beta, e.g. for LENGS=12, beta=2.4 it is enough
c to take LPMAX=10000
c     PARAMETER(LPMAX=10000)  ! for lengs=12, beta=2.4
      PARAMETER(LPMAX=10000)  ! to be checked in RESTLOOP and CORF
      PARAMETER(MAXIS=5000)  ! for lengs=12, beta=2.4, to be checked in RESTLOOP and CORF
c     PARAMETER(MAXIS=50000)  
c     PARAMETER(MNL=2000)
      PARAMETER(MNL=2000)
      
c new line:
c************************************************************
      common /mcurr/ cmag(nsite,nd,3)      ! mag. currents
c************************************************************

      COMMON /NSTEP/  NSTN(4),NSTL(4)
      COMMON /MLRES/  IDI1(MNL),IDI2(MNL)
      COMMON /FSITE/  NFS,IDF,IIF(4)
      COMMON /LOOPP/  NCROS(4),NMONO(4),NQSUM(4),ID(4),
     *NCONT,NCONT1,ISTAT,NPL,NSCONT(MAXIS),IDCONT(MAXIS),
     *IICONT(4,MAXIS)
c     COMMON /CRF/ LSITE(4,LPMAX),LST(LPMAX),LSO2(4),LS(4)
      COMMON /CRF/ LSITE(5,LPMAX),LST(LPMAX),LSO2(4),LS(4)
      COMMON /CLUST/ BCR(2),CHI,DFR1,DFR2,LCMAX,LCLU,FISTAT(6),
     *               FSTAT(6),NISTAT(6),NSTAT(6),LCRF(MAXD)
      common /lgclust/ lgcl(lpmax,5),nplmax

      DIMENSION II1(4),II2(4),II3(4),
     *IDQ(4),IDI(100),ICO(101),NCRO(4),NMON(4),
     *NMON1(4)
      DIMENSION NS5(4),ID5(4),II5(4,4),ID4(4)

      dimension ns(4,nsite)
      dimension nm(4)
      
C


c*************************************************
      i=icolor 
      do mu=1,nd
	nm(mu)=0
      do m=1,nsite
       if(cmag(m,mu,i).ne.0.) then
c       write(8,*) m,mu,i,cmag(m,mu,i)
	nm(mu)=nm(mu)+abs(cmag(m,mu,i))
       endif
      enddo
      enddo
       write(8,*) (nm(mu),mu=1,nd)
c*************************************************
      do idir=1,nd
      do is=1,nsite
	ns(idir,is)=cmag(is,idir,icolor)
      enddo
      enddo

      LS(1)=LENGS
      LS(2)=LENGS
      LS(3)=LENGS
      LS(4)=LENGT
C
      DO 3 J=1,4
 3    LSO2(J)=LS(J)/2

C INITIAL VALUES
      DO 12 I=1,6
      NSTAT(I)=0
      NISTAT(I)=0
      FSTAT(I)=0.
   12 FISTAT(I)=0.

      DO 14 J=1,LPMAX
      LST(J)=0
c     DO 14 I=1,4
      DO 14 I=1,5
  14  LSITE(I,J)=0

      DO 4 I=1,4
      NCRO(I)=0
      NMON1(I)=0
  4   NMON(I)=0

      NSTN(1)=1
      NSTL(1)=LS(1)

      DO 5 J=1,3
      NSTN(J+1)=NSTN(J)*LS(J)
 5    NSTL(J+1)=NSTL(J)*LS(J+1)

      DO 6 I=1,MAXD
   6  LCRF(I)=0

      NSITE1=0
      NT=0
      NT1=0
      NL=0
      NPL=0
      LCLU=0
      LCMAX=0
      CHI=0.D0
      nplmax=0
      ncontmax=0

      DO 1 I1=1,50
      IDI1(I1)=0
   1  IDI2(I1)=0
      K1=5
      K2=0

      DO 111 I1=1,100
  111  IDI(I1)=0
C**********************************************
C LOOP OVER ALL SITES TO FIND ONE OF THE LOOP POINTS

      DO 500 I4=1,LENGT
      II1(4)=I4
      DO 400 I3=1,LENGS
      II1(3)=I3
      DO 399 I2=1,LENGS
      II1(2)=I2
      DO 398 I1=1,LENGS
      II1(1)=I1

      NCLOSE=0
      NSITE1=NSITE1+1

C LOOP OVER ALL DIRECTIONS
      DO 300 ID1=1,4
C NSITE2,II2( ),ID2-CURRENT SITE,ITS COORDINATES AND DIRECTION OF LOOP
C INITIAL VALUES
      NSITE2=NSITE1
      ID2=ID1

      DO 10 J=1,4
      II2(J)=II1(J)
      ID(J)=0
      IDQ(J)=0
      NCROS(J)=0
      NMONO(J)=0
      NQSUM(J)=0
 10   CONTINUE
      NUCR=0
      NCONT=0
      NCONT1=0
      ICONT=0
      LCL=1

C  SELECT START POINT OF A MONOPOLE LOOP
C                        <   =  >
      IF(NS(ID1,NSITE1)) 20,300,30
C SIGN OF THE FIRST MONOP. CHARGE
 20   ISG=-1
      GO TO 50
 30   ISG=+1
C  NPL-LOOP LENGTH , ISTAT-STATICITY
 50   NPL=0
      ISTAT=0

c     WRITE(6,706) NL+1
C  FIRST STEP DONE ALWAYS IN NEGATIV DIRECTION ID2=ID1
 52   CONTINUE
      
C TO STORE LOOP POINT
      LST(NPL+1)=NSITE2
      DO 45 I5=1,4
  45  LSITE(I5,NPL+1)=II2(I5)
       LSITE(5,NPL+1)=ID2

c     WRITE(6,703) (II2(I),I=1,4),ID2

      II2(ID2)=II2(ID2)-1
      NSTEP=-NSTN(ID2)
      IF(II2(ID2) .NE. 0) GO TO 55
      II2(ID2)=LS(ID2)
      NSTEP=NSTEP+NSTL(ID2)
 55   CONTINUE

C TO ANNIGILATE MONOP. AT CURRENT SITE
      NS(ID2,NSITE2)=NS(ID2,NSITE2)-ISG

      NQSUM(ID2)=NQSUM(ID2)+ISG
      IQS=IABS(NQSUM(ID2))
      IF(MOD(IQS,LS(ID2)).EQ.0.AND.IQS.NE.0) THEN
      NQSUM(ID2)=0
      NCROS(ID2)=NCROS(ID2)+1
      ENDIF

      NMONO(ID2)=NMONO(ID2)+1
      NPL=NPL+1
      NSITE2=NSITE2+NSTEP

 60   CONTINUE

C TO FIND OUT NEXT LOOP SITE
      L1=0
      L2=0
      DO 80 ID3=1,4

C END OF LOOP CONDITION
      IF(ID3 .EQ. ID1 .AND. NSITE2 .EQ. NSITE1 .AND. NS(ID3,NSITE2)
     *.EQ. 0)  NCLOSE=1

      IF(ISG*NS(ID3,NSITE2)) 65,65,62
 62   L1=L1+1
      ID4(L1)=ID3

C TO FIND OUT MONOP. AT NEIGHBOUR SITE
 65   II3(ID3)=II2(ID3)+1
      NSTEP=NSTN(ID3)
      IF(II2(ID3) .NE. LS(ID3)) GO TO 70
      II3(ID3)=1
      NSTEP=NSTEP-NSTL(ID3)
 70   CONTINUE

      NSITE3=NSITE2+NSTEP
      IF(ISG*NS(ID3,NSITE3)) 82,80,80

 82   L2=L2+1
      ID5(L2)=ID3
      NS5(L2)=NSITE3

      DO 84 I14=1,4
      II5(I14,L2)=II2(I14)
      IF(I14.EQ.ID3) II5(I14,L2)=II3(I14)
 84   CONTINUE

 80   CONTINUE
C*********************************************************************
      IF(L1 .EQ. 0) GO TO 63
      IF(ID4(1).EQ.ID2) ISTAT=ISTAT+1
      ID2=ID4(1)

C IF THERE IS THE SECOND WAY TO CONTIN. LOOP WITH SAME KIND OF MON.
      IF(L1.EQ.2) THEN
       NCONT=NCONT+1
       NSCONT(NCONT)=NSITE2
       IDCONT(NCONT)=ID4(2)
       DO 66 I6=1,4
  66   IICONT(I6,NCONT)=II2(I6)
       GO TO 64
      ENDIF

  63  CONTINUE

      IF(L2.EQ.0 ) GO TO 64

      IF(L1.GE.1.AND.L2.GE.1) THEN
       NCONT=NCONT+1
       NSCONT(NCONT)=NS5(1)
       IDCONT(NCONT)=ID5(1)
       DO 86 I6=1,4
  86   IICONT(I6,NCONT)=II5(I6,1)
       GO TO 64
      ENDIF

      IF(L2.EQ.2) THEN
       NCONT=NCONT+1
       NSCONT(NCONT)=NS5(2)
       IDCONT(NCONT)=ID5(2)
       DO 87 I6=1,4
  87   IICONT(I6,NCONT)=II5(I6,2)
      ENDIF
C*****************************************************************
C TO STORE LOOP POINT
c     LST(NPL+1)=NSITE2
c     DO 46 I5=1,4
c 46  LSITE(I5,NPL+1)=II2(I5)

      NSITE2=NS5(1)
      IF(ID5(1) .EQ. ID2) ISTAT=ISTAT+1
      ID2=ID5(1)
      II2(ID2)=II5(ID2,1)
C TO STORE LOOP POINT
      LST(NPL+1)=NSITE2
      DO 46 I5=1,4
  46  LSITE(I5,NPL+1)=II2(I5)
      LSITE(5,NPL+1)=ID2

c     WRITE(6,703) (II2(I),I=1,4),ID2

C TO ANNIGILATE MONOP. AT THIS SITE
      NS(ID2,NSITE2)=NS(ID2,NSITE2)+ISG

      NQSUM(ID2)=NQSUM(ID2)-ISG
      IQS=IABS(NQSUM(ID2))
      IF(MOD(IQS,LS(ID2)).EQ.0.AND.IQS.NE.0) THEN
      NQSUM(ID2)=0
      NCROS(ID2)=NCROS(ID2)+1
      ENDIF

      NMONO(ID2)=NMONO(ID2)+1
      NPL=NPL+1

      GO TO 60

  64  CONTINUE
      IF(L1.EQ.0.AND.L2.EQ.0) GO TO 250

C NEXT MONOP. IS OF THE SAME KIND AS PREVIOS
      GO TO 52
C
 250  CONTINUE
      IF(ID3 .EQ. ID2) ISTAT=ISTAT+1

C      WRITE(6,711) NCONT
      NCONT1=NCONT
C****************************************************************
C TO CONTINUE LOOP WITH SELF-INTERSECTIONS
 260  CONTINUE
C      
      IF(NCONT1.EQ.0) GO TO 350
      ICONT=ICONT+1
C
C COORDINATES OF THE FIRST SITE FOR SUBR. RESTLOOP:
      NFS=NSCONT(ICONT)
      IDF=IDCONT(ICONT)
      DO 97 I6=1,4
  97  IIF(I6)=IICONT(I6,ICONT)

C IF LOOP HAVE WALKED THROUGH THIS POINT ALREADY 
      IF(NS(IDF,NFS).EQ.0) THEN
      NCONT1=NCONT1-1
      GO TO 260
      ENDIF

C***********************************************************
C TO COMPUTE LOOP CONTINUATION IF LOOP HAS INTERSECTIONS
C      WRITE(6,349)NFS,IDF,IIF
      CALL RESTLOOP
C***********************************************************
      NCONT1=NCONT1-1
      GO TO 260

C
 350  CONTINUE
C***********************************************************
C TO COMPUTE CORRELATION FUNCTION
c     CALL CORF
C***********************************************************

      NL=NL+1
      NT=NT+NPL

C TO CALCULATE STATICITY FOR DIFFR. TYPES OF CLUSTERS
C***********************************************************
      ICRS1=0
      DO 249 K=1,3
      ICRS1=ICRS1+NCROS(K)
  249 CONTINUE
      ICRS2=NCROS(4)

      STAT=FLOAT(ISTAT)/NPL
C ALL CLUSTERS
      FISTAT(1)=FISTAT(1)+ISTAT
      FSTAT(1)=FSTAT(1)+STAT
      
C ALL WRAP. CLUSTERS
      IF(ICRS1.NE.0.OR.ICRS2.NE.0) THEN
      FISTAT(2)=FISTAT(2)+ISTAT
      FSTAT(2)=FSTAT(2)+STAT
      NISTAT(2)=NISTAT(2)+NPL
      NSTAT(2)=NSTAT(2)+1
      ENDIF

C ALL WRAP. IN TIME-DIR. CLUSTERS
      IF(ICRS2.NE.0) THEN
      FISTAT(3)=FISTAT(3)+ISTAT
      FSTAT(3)=FSTAT(3)+STAT
      NISTAT(3)=NISTAT(3)+NPL
      NSTAT(3)=NSTAT(3)+1
      ENDIF

C ALL WRAP. IN TIME-DIR. ONLY CLUSTERS
      IF(ICRS2.NE.0.AND.ICRS1.EQ.0) THEN
      FISTAT(4)=FISTAT(4)+ISTAT
      FSTAT(4)=FSTAT(4)+STAT
      NISTAT(4)=NISTAT(4)+NPL
      NSTAT(4)=NSTAT(4)+1
      ENDIF

C ALL WRAP. IN TIME-DIR. CLUSTERS WITH LENGTH GREATER THEN 4
      IF(ICRS2.NE.0.AND.ICRS1.EQ.0.AND.NPL.GT.4) THEN
      FISTAT(5)=FISTAT(5)+ISTAT
      FSTAT(5)=FSTAT(5)+STAT
      NISTAT(5)=NISTAT(5)+NPL
      NSTAT(5)=NSTAT(5)+1
      ENDIF

C ALL WRAP. IN SPACE-DIR. ONLY CLUSTERS
      IF(ICRS1.NE.0.AND.ICRS2.EQ.0) THEN
      FISTAT(6)=FISTAT(6)+ISTAT
      FSTAT(6)=FSTAT(6)+STAT
      NISTAT(6)=NISTAT(6)+NPL
      NSTAT(6)=NSTAT(6)+1
      ENDIF
C***********************************************************

C      WRITE(6,700) NPL
C      WRITE(6,710) NQSUM
C      WRITE(6,711) NCONT

C*****************************************************************
C TO CALCULATE GIVEN CLUSTER SIZE
      DO 48 I=2,NPL                 ! this calculation of LCL is probably
      DO 49 J=1,I-1                 ! incorrect  !!!!!
      IF(LST(I).EQ.LST(J)) GO TO 48
  49  CONTINUE
      LCL=LCL+1    
  48  CONTINUE

      IF(LCL.GT.LCMAX) THEN
      LCMAX=LCL
      LLMAX=NPL
      ENDIF

      LCLU=LCLU+LCL
      CHI=CHI+DFLOAT(LCL**2)

      if(npl.gt.nplmax) then
        nplmax=npl
        ncontmax=ncont
        do k=1,5
        do j=1,npl
          lgcl(j,k)=lsite(k,j)
        enddo
        enddo
      endif
c      WRITE(18,712) LCL
c      WRITE(18,*) 'nplmax=',nplmax,' ncontmax=',ncontmax

c      DO 47 J=1,NPL
c     LSI=(LSITE(1,J)-1)*NSTN(1)+(LSITE(2,J)-1)*NSTN(2)+
c    +(LSITE(3,J)-1)*NSTN(3)+(LSITE(4,J)-1)*NSTN(4)+1
c
c  47  WRITE(18,733) (LSITE(I7,J),I7=1,5),LST(J),LSI
C******************************************************
      IT1=0
      DO 251 K=1,4
      IF(IABS(NCROS(K)).EQ.0) THEN
      IT1=IT1+1
      ENDIF
  251 CONTINUE

      IF(IT1.EQ.4) THEN
      NT1=NT1+NPL
      DO 252 K=1,4
 252  NMON1(K)=NMON1(K)+NMONO(K)
      ENDIF
C******************************************************
C
C      WRITE(6,720) NCROS
      DO 310 I=1,4
      NCRO(I)=NCRO(I)+NCROS(I)
 310  NMON(I)=NMON(I)+NMONO(I)

C******************************************************
      ICRS=0
      DO 311 J1=1,4
 311  ICRS=ICRS+NCROS(J1)

      IF(ICRS.NE.0) GO TO 320
      IF(NPL .LE. 12) THEN
      J=NPL/2-1
      IDI1(J)=IDI1(J)+1
      ELSE
      K1=K1+1
      IDI1(K1)=NPL
      ENDIF

  320 CONTINUE
      IF(ICRS.EQ.0) GO TO 321
      K2=K2+1
      IDI2(K2)=NPL
  321 CONTINUE
C******************************************************
  

 300  CONTINUE
 398  CONTINUE
 399  CONTINUE
 400  CONTINUE
 500  CONTINUE
       WRITE(18,*) nplmax,ncontmax

C  CLUSTERS FEATURES
C*******************************************************************
C CHI - SUSCEPTIBILITY IS DEFINED AS IN PH.REV.LETT.,63(1989)2169
C DFR - FRACTAL DIMENSION IS DEFINED AS IN ITEP-71-90
C*******************************************************************
      IF(LCLU.NE.0) THEN

      CHI=(CHI-LCMAX**2)/LCLU
      DFR1=LLMAX/FLOAT(LCMAX)
      DFR2=NT/FLOAT(LCLU)

      FSTAT(1)=FSTAT(1)/NL
      FISTAT(1)=FISTAT(1)/NT
      DO 599 I=2,6
      IF(NSTAT(I).NE.0)THEN
      FSTAT(I)=FSTAT(I)/NSTAT(I)
      FISTAT(I)=FISTAT(I)/NISTAT(I)
      ENDIF
  599 CONTINUE
  
      ELSE

      CHI=0.
      DFR1=1.
      DFR2=1.

      DO 699 I=1,6
      FSTAT(I)=0.
      FISTAT(I)=0.
  699 CONTINUE
      ENDIF

  800 CONTINUE

      BCR(1)=0.
      DO 600 I=1,3
  600 BCR(1)=BCR(1)+NCRO(I)
      BCR(2)=NCRO(4)
      BCR(1)=BCR(1)/3.

      WRITE(6,701) NCRO,NMON,NT,NT1,FSTAT(1),FISTAT(1)
C      WRITE(6,702) FSTAT,FISTAT,NSTAT,NISTAT
C      WRITE(6,725) IDI1,IDI2
       WRITE(6,715) LCRF
      WRITE(6,713) LCLU,LCMAX,CHI
      WRITE(6,714) DFR1,DFR2
      WRITE(6,716) BCR
      WRITE(17,*) BCR

c largest cluster:
c     do k1=1,nplmax
c      write(19,*) (lgcl(k1,k2),k2=1,5)
c     enddo
C
  162 FORMAT(/3X,'L1= ',I4/)
  182 FORMAT(/3X,'L2= ',I4/)
C      WRITE(6,349)NFS,IDF,IIF
  349 FORMAT(/2X,'CONTINUATION',2X,'NFS=',I5,2X,'IDF=',I2,2X,
     *'IIF=',4I3/)
  706 FORMAT(30X,'--LOOPS NUMBER=',I4,'--'/)
  700 FORMAT(3X,'LOOPS LENGTH=',I6/)
  710 FORMAT(3X,'DIFFER. BETWEEN NUMBER OF M AND AM =',4I6/)
  711 FORMAT(3X,'NUMBER OF POINTS WITH SELF-INTERS.=',I6/)
  712 FORMAT(3X,'CLUSTER SIZE=',I6/)
  713 FORMAT(3X,'CLUSTER SIZES SUMED UP=',I4,
     *3X,'MAX CLUSTER SIZE=',I4/3X,'SESCEPTIBILITY=',F9.6/)
  714 FORMAT(3X,'DIM. OF FRACT.=',2F9.4/)
  716 FORMAT(3X,'NUMBER OF BCR.=',2F9.4/)
  720 FORMAT(3X,'NUMBER OF BOUND. CROSS.=',4I6//)
  701 FORMAT(//30X,'TOTAL RESULTS FOR THIS ABEL. CONFIG.'//
     *3X,'TOTAL NUMBER OF BOUND CROSS=',4I6//
     *3X,'TOTAL NUMBER OF MONOP. AND ANTIMON.=',4I6/
     *3X,'TOTAL LOOPS LENGTH=',I6//
     *3X,'TOTAL LOOPS LENGTH WITHOUT BOUND.CROSSINGS=',I6//
     *2X,'STAT=',F6.3,2X,'INT. STAT.=',F6.3/)     
  702 FORMAT(/3X,'STATICITY:'/2X,6F6.3,3X,6F6.3/2X,6I6,3X,6I6/)
  703 FORMAT(3X,4I4,6X,I3)
  733 FORMAT(3X,5I4,6X,I5,6X,I5)
  715 FORMAT(//3X,'CORRELATION FUNCTION'//5(3X,10I6/),3X,2I6//)
  725 FORMAT(//3X,'DISTRIBUTION OF LOOPS LENGTH',/2(5(3X,10I5/)/)/)
  721 FORMAT(2X,'UNCLOSED LOOP'/)
      RETURN
      END


C****************************************************************** 
C****************************************************************** 
      SUBROUTINE RESTLOOP
C****************************************************************** 
C TO FIND LOOP CONTINUATION
C****************************************************************** 
c     include 'param.txt'
c     include 'su2mon.cmn'
      include'paravp3'
      parameter(lengs=n1)
      parameter(lengt=n4)

      PARAMETER(MAXD=3*(LENGS/2)**2+(LENGT/2)**2)
      PARAMETER(LPMAX=10000)
      PARAMETER(MAXIS=5000)

      COMMON /NSTEP/  NSTN(4),NSTL(4)
      COMMON /FSITE/  NSITE1,ID1,II1(4)
      COMMON /LOOPP/  NCROS(4),NMONO(4),NQSUM(4),ID(4),
     *NCONT,NCONT1,ISTAT,NPL,NSCONT(MAXIS),IDCONT(MAXIS),
     *IICONT(4,MAXIS)
c     COMMON /CRF/ LSITE(4,LPMAX),LST(LPMAX),LSO2(4),LS(4)
      COMMON /CRF/ LSITE(5,LPMAX),LST(LPMAX),LSO2(4),LS(4)

      DIMENSION II2(4),II3(4),
     *IDQ(4),IDI(100),ICO(101),NCRO(4),NMON(4),
     *NMON1(4)
      DIMENSION NS5(4),ID5(4),II5(4,4),ID4(4)
      dimension ns(4,nsite)
     

C
      NSITE2=NSITE1
      ID2=ID1
      NCLOSE=1

      DO 10 J=1,4
  10  II2(J)=II1(J)

      IF(NS(ID1,NSITE1)) 20,250,30
C SIGN OF THE FIRST MONOP. CHARGE
 20   ISG=-1
      GO TO 50
 30   ISG=+1
 50   CONTINUE

C  FIRST STEP DONE ALWAYS IN NEGATIV DIRECTION ID2=ID1
 52   CONTINUE

      LST(NPL+1)=NSITE2
      DO 45 I5=1,4
  45  LSITE(I5,NPL+1)=II2(I5)
       LSITE(5,NPL+1)=ID2

c     WRITE(6,703) (II2(I),I=1,4),ID2

      II2(ID2)=II2(ID2)-1
      NSTEP=-NSTN(ID2)
      IF(II2(ID2) .NE. 0) GO TO 55
      II2(ID2)=LS(ID2)
      NSTEP=NSTEP+NSTL(ID2)
 55   CONTINUE

C TO ANNIGILATE MONOP. AT CURRENT SITE
      NS(ID2,NSITE2)=NS(ID2,NSITE2)-ISG

      NQSUM(ID2)=NQSUM(ID2)+ISG
      IQS=IABS(NQSUM(ID2))
      IF(MOD(IQS,LS(ID2)).EQ.0.AND.IQS.NE.0) THEN
      NQSUM(ID2)=0
      NCROS(ID2)=NCROS(ID2)+1
      ENDIF

      NMONO(ID2)=NMONO(ID2)+1
      NPL=NPL+1
      NSITE2=NSITE2+NSTEP

 60   CONTINUE

C TO FIND OUT NEXT LOOP SITE
      L1=0
      L2=0

      DO 80 ID3=1,4

C END OF LOOP CONDITION
      IF(ID3 .EQ. ID1 .AND. NSITE2 .EQ. NSITE1 .AND. NS(ID3,NSITE2)
     *.EQ. 0)  NCLOSE=1

      IF(ISG*NS(ID3,NSITE2)) 65,65,62
 62   L1=L1+1
      ID4(L1)=ID3
C      WRITE(6,162) L1      

C      GO TO 80

C TO FIND OUT MONOP. AT NEIGHBOUR SITE
 65   II3(ID3)=II2(ID3)+1
      NSTEP=NSTN(ID3)
      IF(II2(ID3) .NE. LS(ID3)) GO TO 70
      II3(ID3)=1
      NSTEP=NSTEP-NSTL(ID3)
 70   CONTINUE

      NSITE3=NSITE2+NSTEP
      IF(ISG*NS(ID3,NSITE3)) 82,80,80

 82   L2=L2+1
      ID5(L2)=ID3
      NS5(L2)=NSITE3
C      WRITE(6,182) L2      

      DO 84 I14=1,4
      II5(I14,L2)=II2(I14)
      IF(I14.EQ.ID3) II5(I14,L2)=II3(I14)
 84   CONTINUE

 80   CONTINUE

      IF(L1 .EQ. 0) GO TO 63
      IF(ID4(1).EQ.ID2) ISTAT=ISTAT+1
      ID2=ID4(1)

C IF THERE IS THE SECOND WAY TO CONTIN. LOOP WITH SAME KIND OF MON.
      IF(L1.EQ.2) THEN
      NCONT=NCONT+1
      NCONT1=NCONT1+1
      NSCONT(NCONT)=NSITE2
      IDCONT(NCONT)=ID4(2)
      DO 66 I6=1,4
  66  IICONT(I6,NCONT)=II2(I6)

      GO TO 64
      ENDIF

  63  CONTINUE

      IF(L2.EQ.0 ) GO TO 64

      IF(L1.GE.1.AND.L2.GE.1) THEN
      NCONT=NCONT+1
      NCONT1=NCONT1+1
      NSCONT(NCONT)=NS5(1)
      IDCONT(NCONT)=ID5(1)
      DO 86 I6=1,4
  86  IICONT(I6,NCONT)=II5(I6,1)
      GO TO 64
      ENDIF

      IF(L2.EQ.2) THEN
       NCONT=NCONT+1
       NCONT1=NCONT1+1
       NSCONT(NCONT)=NS5(2)
       IDCONT(NCONT)=ID5(2)
       DO 87 I6=1,4
  87   IICONT(I6,NCONT)=II5(I6,2)
      ENDIF

c     LST(NPL+1)=NSITE2
c     DO 46 I5=1,4
c 46  LSITE(I5,NPL+1)=II2(I5)

      NSITE2=NS5(1)
      IF(ID5(1) .EQ. ID2) ISTAT=ISTAT+1
      ID2=ID5(1)
      II2(ID2)=II5(ID2,1)
      LST(NPL+1)=NSITE2
      DO 46 I5=1,4
  46  LSITE(I5,NPL+1)=II2(I5)
      LSITE(5,NPL+1)=ID2

c     WRITE(6,703) (II2(I),I=1,4),ID2

C TO ANNIGILATE MONOP. AT THIS SITE
      NS(ID2,NSITE2)=NS(ID2,NSITE2)+ISG
      NQSUM(ID2)=NQSUM(ID2)-ISG
      IQS=IABS(NQSUM(ID2))
      IF(MOD(IQS,LS(ID2)).EQ.0.AND.IQS.NE.0) THEN
      NQSUM(ID2)=0
      NCROS(ID2)=NCROS(ID2)+1
      ENDIF

      NMONO(ID2)=NMONO(ID2)+1
      NPL=NPL+1

      GO TO 60

  64  CONTINUE
      IF(L1.EQ.0.AND.L2.EQ.0) GO TO 250

C NEXT MONOP. IS OF THE SAME KIND AS PREVIOS
      GO TO 52
C
 250  CONTINUE
C      WRITE(6,711) NCONT

  703 FORMAT(3X,4I4,6X,I3)
  711 FORMAT(/3X,'NUMBER OF POINTS WITH SELF-INTERS.=',I6//)
  162 FORMAT(/3X,'L1= ',I4/)
  182 FORMAT(/3X,'L2= ',I4/)
  
      RETURN
      END
C*******************************************************************
      SUBROUTINE CORF
C TO CALCULATE CORRELATION FUNCTION LCRF(R)
C*******************************************************************
c     include 'param.txt'
c     include 'su2mon.cmn'
      include'paravp3'

      parameter(lengs=n1)
      parameter(lengt=n4)
      PARAMETER(MAXD=3*(LENGS/2)**2+(LENGT/2)**2)
      PARAMETER(LPMAX=10000)
      PARAMETER(MAXIS=5000)

c     COMMON /CRF/ LSITE(4,LPMAX),LST(LPMAX),LSO2(4),LS(4)
      COMMON /CRF/ LSITE(5,LPMAX),LST(LPMAX),LSO2(4),LS(4)
      COMMON /LOOPP/  NCROS(4),NMONO(4),NQSUM(4),ID(4),
     *NCONT,NCONT1,ISTAT,NPL,NSCONT(MAXIS),IDCONT(MAXIS),
     *IICONT(4,MAXIS)
      COMMON /CLUST/ BCR(2),CHI,DFR1,DFR2,LCMAX,LCLU,FISTAT(6),
     *               FSTAT(6),NISTAT(6),NSTAT(6),LCRF(MAXD)

      DIMENSION NCRF(MAXD)

      DO 6 I=1,MAXD
   6  NCRF(I)=0

      DO 1 I1=1,NPL-1
      DO 1 I2=I1+1,NPL

      NDI=0

      DO 3 I3=1,4
      IDI=ABS(LSITE(I3,I1)-LSITE(I3,I2))
      IF(IDI.GT.LSO2(I3)) IDI=LS(I3)-IDI
  3   NDI=NDI+IDI**2

      IF(NDI.EQ.0) THEN
      NIS=NIS+1
      GO TO 1
      ENDIF

      NCRF(NDI)=NCRF(NDI)+1
      LCRF(NDI)=LCRF(NDI)+1

  1   CONTINUE
C      WRITE(6,715) NCRF

  715 FORMAT(//3X,'CORRELATION FUNCTION'//5(3X,10I6/),3X,2I6//)

      RETURN
      END
C*************************************************************************
