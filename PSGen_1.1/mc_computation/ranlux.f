      SUBROUTINE RANLUX(RVEC,LENV)
C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C     Fortran 77 coded by F. James, 1993
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT) initializes the generator from one    ++
C!!!                   32-bit integer INT and sets Luxury Level LUX  ++
C!!!                   which is integer between zero and MAXLEV,     ++
C!!!                   or if LUX .GT. 24, it sets p=LUX directly     ++
C!!!      CALL RLUXIN(IVEC)    restarts the generator from vector    ++
C!!!                   IVEC of 25 32-bit integers (see RLUXUT)       ++
C!!!      CALL RLUXUT(IVEC)    outputs the current values of the 25  ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      PARAMETER (TWOP12=4096.)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24
      INTEGER LUXLEV
      LOGICAL NOTYET
      DATA NOTYET, LUXLEV, IN24 /.TRUE., LXDFLT, 0/
      DATA I24,J24,CARRY/24,10,0./
C
      DATA NDSKIP/ 0,  24,  73, 199, 365 /
Corresponds to p= 24   48   97  223  389
C     time factor  1    2    4    9   16   (worst case)
C                  1   1.5   2    3    5   (in practice, approx.)
C
C  NOTYET is .TRUE. if no initialization has been performed yet.
C              Default Initialization by Multiplicative Congruential

      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = 314159265
ckk>
c         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
ckk<
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
ckk>
c         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
c     +        LUXLEV,'      p =',LP
ckk<
            TWOM24 = 1.
         DO I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
         end do

         TWOM12 = TWOM24 * 4096.
         DO I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         end do
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .LE. SEEDS(14)) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO IVEC= 1, LENV
CCCCCC     This has now been corrected for the bug discovered by Luescher  ***
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY                               <<<<
CCCCCC     This has now been corrected for the bug discovered by Luescher  ***
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = I24 - 1
      IF (I24 .EQ. 0)  I24 = 24
      J24 = J24 - 1
      IF (J24 .EQ. 0)  J24 = 24
      RVEC(IVEC) = UNI
C  small numbers (with less than 12 "significant" bits) are "padded".
      IF (UNI .LT. TWOM12)  THEN
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C        and zero is forbidden in case someone takes a logarithm
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
      ENDIF
C        Skipping to luxury.  As proposed by Martin Luscher.
      IN24 = IN24 + 1
      IF (IN24 .EQ. 24)  THEN
         IN24 = 0
         DO ISK= 1, NSKIP
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
CCCCCC         The above has been corrected
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = I24 - 1
         IF (I24 .EQ. 0)  I24 = 24
         J24 = J24 - 1
         IF (J24 .EQ. 0)  J24 = 24
         end do
      ENDIF
      end do
      RETURN
C
C           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
         TWOM24 = 1.
         DO I= 1, 24
          TWOM24 = TWOM24 * 0.5
         end do
         TWOM12 = TWOM24 * 4096.
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      DO I= 1, 24
       SEEDS(I) = REAL(ISDEXT(I))*TWOM24
      end do
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
        IF (LUXLEV .LE. MAXLEV) THEN
          NSKIP = NDSKIP(LUXLEV)
          WRITE (6,'(A,I2)') ' RANLUX luxury level set to ',LUXLEV
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX p-value set to ',LUXLEV
        ELSE
          NSKIP = NDSKIP(MAXLEV)
          WRITE (6,'(A,I5)') ' RANLUX illegal input luxury: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
      RETURN
C
C                    Entry to ouput seeds as integers
      ENTRY RLUXUT(ISDEXT)
      DO I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
      end do
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
      RETURN
C
C                    Entry to initialize from one integer
      ENTRY RLUXGO(LUX,INSEED)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
            WRITE (6,'(A,I7)') ' RANLUX illegal input luxury: ',LUX
         ELSE
            LUXLEV = LUX
            DO ILX= 0, MAXLEV
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
            end do
         ENDIF
      IF (LUXLEV .LE. MAXLEV)  THEN
         NSKIP = NDSKIP(LUXLEV)
         WRITE(6,'(A,I2,A,I4)')   ' RANLUX luxury level set to ',
     +        LUXLEV,'      p =', NSKIP+24
      ELSE
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX p-value set to ',LUXLEV
      ENDIF
      IN24 = 0
      IF (INSEED .NE. 0)  THEN
        JSEED = INSEED
        WRITE(6,'(A,I12)') ' RANLUX INITIALIZED FROM SEED ',INSEED
      ELSE
        JSEED = 314159265
        WRITE(6,'(A)') ' RANLUX INITIALIZED FROM DEFAULT SEED '
      ENDIF
      NOTYET = .FALSE.
      TWOM24 = 1.
         DO I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
         end do
         TWOM12 = TWOM24 * 4096.
         DO I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         end do
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .LE. SEEDS(14)) CARRY = TWOM24
      RETURN
      END

