      subroutine vsup(nd,npt,its,fxn,avgi,sd,chisq)
      implicit double precision (a-h,o-z)
      character *25 fname8
      logical plot
      common/bvegas/ndim,ncall,nprn,itmx
      common/plots/plot
      common/iparallel/ilabel
      common/chfile/fname8
      external fxn
      ndim=nd
      ncall=npt
      nprn = 1
      itmx=its
      plot=.false.
      call vegas_m(1,fxn,avgi,sd,chisq)
      do i=1,itmx
         call vegas_m(3,fxn,avgi,sd,chisq)
      enddo
      ncall=10*npt
      itmx=2*its
      plot=.true.
      call vegas_m(2,fxn,avgi,sd,chisq)
      do i=1,itmx
         call vegas_m(4,fxn,avgi,sd,chisq)
      enddo
      return
      end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ Plus specific vsup]      
      subroutine vsup_2(nd,npt,its,fxn,avgi,sd,chisq)
      implicit double precision (a-h,o-z)
      character *25 fname8
      logical plot
      common/bvegas/ndim,ncall,nprn,itmx
      common/plots/plot
      common/iparallel/ilabel
      common/chfile/fname8
      external fxn
      ndim=nd
      ncall=npt*10
      nprn = 1
      itmx=its
      plot=.false.
      call vegas_m(1,fxn,avgi,sd,chisq)
      do i=1,itmx
         call vegas_m(3,fxn,avgi,sd,chisq)
      enddo
      ncall=ncall*10
      itmx=its*2
      plot=.true.
      call vegas_m(2,fxn,avgi,sd,chisq)
      do i=1,itmx
         call vegas_m(4,fxn,avgi,sd,chisq)
      enddo
      return
      end

      SUBROUTINE VEGAS_M(ISTAT,FXN,AVGI,SD,CHI2A)
C
C     PERFORMS N-DIMENSIONAL MONTE CARLO INTEGRATION
C     BY G.P. LEPAGE   SEPT 1976/(REV)APR 1978
C
      IMPLICIT REAL*8 (A-H,O-Z)
      character *25 fname8
      COMMON/BVEGAS/NDIM,NCALL,NPRN,ITMX
      common/iparallel/ilabel
      common/chfile/fname8
      DIMENSION XI(50,10),D(50,10),DI(50,10),XIN(50),R(50),DT(10),X(10),
     1     KG(10),IA(10),QRAN(10)
      EXTERNAL FXN
      SAVE
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1D0/,MDS/1/
C
      IF(ISTAT.EQ.1.OR.ISTAT.EQ.2)THEN
C     
C     INITIALIZE CUMULATIVE VARIABLES 
C     
         IT=0
         SI=0D0
         SI2=0D0
         SWGT=0D0
         SCHI=0D0
         ND=NDMX
         NG=1
         IF(MDS.NE.0)THEN 
            NG=(0.5D0*NCALL)**(ONE/NDIM)
            MDS=1
            IF((2*NG-NDMX).GE.0)THEN
               MDS=-1
               NPG=NG/NDMX+1
               ND=NG/NPG
               NG=NPG*ND
            ENDIF
         ENDIF
         K=NG**NDIM
         NPG=NCALL/K
         IF(NPG.LT.2)NPG=2
         CALLS=NPG*K
         DXG=ONE/NG
         DV2G=(CALLS*DXG**NDIM)**2/NPG/NPG/(NPG-ONE)
         XND=ND
         NDM=ND-1
         DXG=DXG*XND
         XJAC=ONE/CALLS
         IF(NPRN.NE.0) WRITE(6,200) NDIM,CALLS,IT,ITMX,MDS,ND
c         IF(NPRN.NE.0) WRITE(ilabel,200) NDIM,CALLS,IT,ITMX,MDS,ND

      ENDIF
      IF(ISTAT.EQ.1)THEN
C
C   CONSTRUCT UNIFORM GRID  
C
         RC=ONE/XND
         DO J=1,NDIM
            XI(1,J)=ONE
            K=0
            XN=0D0
            DR=0D0
            I=0
 4          K=K+1
            DR=DR+ONE
            XO=XN
            XN=XI(K,J)
 5          IF(RC.GT.DR) GO TO 4
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR
            IF(I.LT.NDM) GO TO 5
            DO  I=1,NDM
               XI(I,J)=XIN(I)
            ENDDO
            XI(ND,J)=ONE
         ENDDO
         NDO=ND
         RETURN
C
      ELSEIF(ISTAT.EQ.2)THEN
C
C   RESCALE REFINED GRID TO NEW ND VALUE - PRESERVE BIN DENSITY
C
         IF(ND.NE.NDO)THEN
            RC=NDO/XND
            DO J=1,NDIM
               K=0
               XN=0D0
               DR=0D0
               I=0
 6             K=K+1
               DR=DR+ONE
               XO=XN
               XN=XI(K,J)
 7             IF(RC.GT.DR) GO TO 6
               I=I+1
               DR=DR-RC
               XIN(I)=XN-(XN-XO)*DR
               IF(I.LT.NDM) GO TO 7
               DO  I=1,NDM
                  XI(I,J)=XIN(I)
               ENDDO
               XI(ND,J)=ONE
            ENDDO
         ENDIF
         RETURN
C     
      ELSEIF(ISTAT.EQ.3.OR.ISTAT.EQ.4)THEN
C     
C     MAIN INTEGRATION LOOP
C         
         IT=IT+1
         TI =0D0
         TSI=0D0
         DO J=1,NDIM
            KG(J)=1
            DO I=1,ND
               D (I,J)=0D0
               DI(I,J)=0D0
            ENDDO
         ENDDO
C     
 11      FB=0D0
         F2B=0D0
         K=0
 12      K=K+1
         CALL BRM48(QRAN,NDIM)
         WGT=XJAC
         DO J=1,NDIM
            XN=(KG(J)-QRAN(J))*DXG+ONE
            IA(J)=XN
            IF(IA(J).GT.1)THEN
               XO=XI(IA(J),J)-XI(IA(J)-1,J)
               RC=XI(IA(J)-1,J)+(XN-IA(J))*XO
            ELSE
               XO=XI(IA(J),J)
               RC=(XN-IA(J))*XO
            ENDIF
            X(J)=RC
            WGT=WGT*XO*XND
         ENDDO
C     
         F=FXN(X,WGT/ITMX)*WGT
         F2=F*F
         FB=FB+F
         F2B=F2B+F2
         DO J=1,NDIM
            DI(IA(J),J)=DI(IA(J),J)+F
            IF(MDS.GE.0) D(IA(J),J)=D(IA(J),J)+F2
         ENDDO
         IF(K.LT.NPG) GO TO 12
C     
         F2B=DSQRT(F2B*NPG)
         F2B=(F2B-FB)*(F2B+FB)
         TI=TI+FB
         TSI=TSI+F2B
         IF(MDS.LT.0) THEN
            DO J=1,NDIM
               D(IA(J),J)=D(IA(J),J)+F2B
            ENDDO
         ENDIF
         K=NDIM
 19      KG(K)=MOD(KG(K),NG)+1
         IF(KG(K).NE.1) GO TO 11
         K=K-1
         IF(K.GT.0) GO TO 19
C     
C     FINAL RESULTS FOR THIS ITERATION
C     
         TSI=TSI*DV2G
         TI2=TI*TI
         WGT=TI2/TSI
         SI=SI+TI*WGT
         SI2=SI2+TI2
         SWGT=SWGT+WGT
         SCHI=SCHI+TI2*WGT
         AVGI=SI/SWGT
         SD=SWGT*IT/SI2
         CHI2A=0D0
         IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
         SD=DSQRT(ONE/SD)
C     
         IF(NPRN.NE.0) THEN
            TSI=DSQRT(TSI)
            WRITE(6,201) IT,TI,TSI,AVGI,SD,CHI2A
c            WRITE(ilabel,201) IT,TI,TSI,AVGI,SD,CHI2A

         ENDIF
         IF(NPRN.LT.0) THEN
            DO J=1,NDIM
            WRITE(6,202) J,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
c               WRITE(ilabel,202) J,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
            ENDDO
         ENDIF
      ENDIF
C     
C     REFINE GRID
C     
      IF(ISTAT.EQ.3)THEN
         DO J=1,NDIM
            XO=D(1,J)
            XN=D(2,J)
            D(1,J)=0.5D0*(XO+XN)
            DT(J)=D(1,J)
            DO I=2,NDM
               D(I,J)=XO+XN
               XO=XN
               XN=D(I+1,J)
               D(I,J)=(D(I,J)+XN)/3D0
               DT(J)=DT(J)+D(I,J)
            ENDDO
            D(ND,J)=0.5D0*(XN+XO)
            DT(J)=DT(J)+D(ND,J)
         ENDDO
C     
         DO 28 J=1,NDIM
            RC=0D0
            DO 24 I=1,ND
               R(I)=0D0
               IF(D(I,J).LE.0D0) GO TO 24
               XO=DT(J)/D(I,J)
               R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
 24            RC=RC+R(I)
               RC=RC/XND
               K=0
               XN=0D0
               DR=XN
               I=K
 25            K=K+1
               DR=DR+R(K)
               XO=XN
               XN=XI(K,J)
 26            IF(RC.GT.DR) GO TO 25
               I=I+1
               DR=DR-RC
               XIN(I)=XN-(XN-XO)*DR/R(K)
               IF(I.LT.NDM) GO TO 26
               DO 27 I=1,NDM
 27            XI(I,J)=XIN(I)
 28      XI(ND,J)=ONE
         RETURN
      ENDIF
C     
 200  FORMAT(///' INPUT PARAMETERS FOR VEGAS:  NDIM=',I3,'  NCALL=',
     1    F12.0/28X,'  IT=',I5,'  ITMX=',I5
     2    /28X,'  MDS=',I3,'   ND=',I4)
201   FORMAT(i4,4x,g15.7,g11.4,g15.7,g11.4,g11.4)
202   FORMAT(' DATA FOR AXIS',I2 / ' ',6X,'X',7X,'  DELT I  ',
     1    2X,' CONV''CE  ',11X,'X',7X,'  DELT I  ',2X,' CONV''CE  '
     2   ,11X,'X',7X,'  DELT I  ',2X,' CONV''CE  ' /
     2    (' ',3G12.4,5X,3G12.4,5X,3G12.4))
      RETURN
      END

      SUBROUTINE BRM48(RVEC,LENV)
C     Double-precision version of
C Universal random number generator proposed by Marsaglia and Zaman
C in report FSU-SCRI-87-50
C        based on RANMAR, modified by F. James, to generate vectors
C        of pseudorandom numbers RVEC of length LENV, where the numbers
C        in RVEC are numbers with at least 48-bit mantissas.
C   Input and output entry points: BRM48I, BRM48O
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RM48:                                    ++
C!!!      CALL RM48 (RVEC, LEN)     returns a vector RVEC of LEN     ++
C!!!                   64-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RM48IN(I1,N1,N2)   initializes the generator from one ++
C!!!                   64-bit integer I1, and number counts N1,N2    ++
C!!!                  (for initializing, set N1=N2=0, but to restart ++
C!!!                    a previously generated sequence, use values  ++
C!!!                    output by RM48UT)                            ++
C!!!      CALL RM48UT(I1,N1,N2)   outputs the value of the original  ++
C!!!                  seed and the two number counts, to be used     ++
C!!!                  for restarting by initializing to I1 and       ++
C!!!                  skipping N2*100000000+N1 numbers.              ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C for 32-bit machines, use IMPLICIT DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character *25 fname8
      DIMENSION RVEC(*)
      COMMON/R48ST1/U(97),C,I97,J97
      common/iparallel/ilabel
      common/chfile/fname8
      PARAMETER (MODCNS=1000000000)
      SAVE CD, CM, TWOM24, NTOT, NTOT2, IJKL,TWOM49
      SAVE ZERO, ONE
      DATA NTOT,NTOT2,IJKL/-1,0,0/
C
      IF (NTOT .GE. 0)  GO TO 50
C
C        Default initialization. User has called RM48 without RM48IN.
      IJKL = 54217137
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
C
      ENTRY      BRM48I(IJKLIN, NTOTIN,NTOT2N)
C         Initializing routine for RM48, may be called before
C         generating pseudorandom numbers with RM48.   The input
C         values should be in the ranges:  0<=IJKLIN<=900 OOO OOO
C                                          0<=NTOTIN<=999 999 999
C                                          0<=NTOT2N<<999 999 999!
C To get the standard values in Marsaglia's paper, IJKLIN=54217137
C                                            NTOTIN,NTOT2N=0
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
C          always come here to initialize
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
      WRITE(6,'(A,I10,2X,2I10)') ' RM48 INITIALIZED:',IJKL,NTOT,NTOT2
c      WRITE(ilabel,'(A,I10,2X,2I10)') ' RM48 INITIALIZED:',
c     &      IJKL,NTOT,NTOT2
CCC      PRINT '(A,4I10)', '   I,J,K,L= ',I,J,K,L
      ONE = 1.
      HALF = 0.5
      ZERO = 0.
      DO 2 II= 1, 97
      S = 0.
      T = HALF
      DO 3 JJ= 1, 48
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
    3    T = HALF*T
    2 U(II) = S
      TWOM49 = T
      TWOM24 = ONE
      DO 4 I24= 1, 24
    4 TWOM24 = HALF*TWOM24
      C  =   362436.*TWOM24
      CD =  7654321.*TWOM24
      CM = 16777213.*TWOM24
      I97 = 97
      J97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT
      IF (NOW .GT. 0)  THEN
      WRITE(6,'(A,I15)') ' RM48IN SKIPPING OVER ',NOW
          DO 40 IDUM = 1, NTOT
          UNI = U(I97)-U(J97)
          IF (UNI .LT. ZERO)  UNI=UNI+ONE
          U(I97) = UNI
          I97 = I97-1
          IF (I97 .EQ. 0)  I97=97
          J97 = J97-1
          IF (J97 .EQ. 0)  J97=97
          C = C - CD
          IF (C .LT. ZERO)  C=C+CM
   40     CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED .EQ. 1)  RETURN
C
C          Normal entry to generate LENV random numbers
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI .LT. ZERO)  UNI=UNI+ONE
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. ZERO)  C=C+CM
      UNI = UNI-C
      IF (UNI .LT. ZERO) UNI=UNI+ONE
      RVEC(IVEC) = UNI
C             Replace exact zeros by 2**-49
         IF (UNI .EQ. ZERO)  THEN
            RVEC(IVEC) = TWOM49
         ENDIF
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT .GE. MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
C           Entry to output current status
      ENTRY BRM48O(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
      RETURN
      END
