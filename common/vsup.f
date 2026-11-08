c#####################################################################
c  VEGAS driver: an adaptation pass followed by a frozen-grid
c  production pass.
c
c  ADAPTATION  (npt_ad x its_ad, ISTAT=3)  refines the grid after every
c  iteration. Nothing is kept: ISTAT=2 below resets the cumulative
c  variables, so the quoted result comes from the production pass alone.
c
c  PRODUCTION  (npt_pr x its_pr, ISTAT=4)  leaves the grid alone. This
c  is the only pass that may be histogrammed -- adaptation points come
c  from a different and still-moving density -- which is what the two
c  hooks below mark out:
c
c     histo_begin_production   arms the fill and zeroes the accumulators
c     histo_end_iteration      closes one iteration
c     checkpoint_write_state   snapshots grid + RNG + accumulators
c
c  checkpoint_write_state is called unconditionally, like the flush
c  below: whether a checkpoint is actually wanted is checkpoint_mod's
c  own business, not a second flag threaded through here. Its natural
c  home is this same boundary, because an iteration is the only point
c  at which the histogram accumulators are whole.
c
c  One iteration is the sample unit for every histogram error bar, so
c  the hook has to sit exactly here, at the iteration boundary, and not
c  inside vegas_m. See analysis_mod.f90.
c
c  THE SIZES ARE ARGUMENTS, not multiples. This routine used to hardcode
c  the production pass as 10*npt points and 2*its iterations, which made
c  the statistics behind a plot something you had to read the integrator
c  to know. The run card states all four now.
c#####################################################################
      subroutine vsup(nd,npt_ad,its_ad,npt_pr,its_pr,epsrel,fxn,avgi,
     &                sd,chisq,converged,resume,start_it)
      implicit double precision (a-h,o-z)
      character *25 fname8
c  vs_done is a LOGICAL FUNCTION. Without this declaration the implicit
c  rule above would make it double precision and the exit test would
c  never compile correctly.
      logical plot,resume,converged,vs_done
c  start_it would otherwise be REAL*8 under the implicit rule above.
      integer start_it
      common/bvegas/ndim,ncall,nprn,itmx
      common/plots/plot
      common/iparallel/ilabel
      common/chfile/fname8
      external fxn
      ndim=nd
      nprn = 1
c  Set on every path before either production loop, so a caller can always
c  trust it: .true. only when the epsrel test actually fired.
      converged=.false.
c
c  RESUME skips straight to the production loop. Everything the skipped
c  calls would have set -- the grid, the running sums, the histogram
c  accumulators, the RNG position -- was restored by checkpoint_resume
c  before this routine was entered, so replaying vegas_m(1/2) or
c  histo_begin_production here would overwrite it with a fresh start.
c  See the ISTAT=2 note at the checkpoint hooks in vegas_m.
c
      if (resume) then
         ncall=npt_pr
         itmx=its_pr
         plot=.true.
         write(6,'(/a,i0,a,i0)')
     &      ' [ VEGAS ] production pass -- RESUMED at iteration ',
     &      start_it+1, ' of ', itmx
         flush(6)
         if (epsrel.gt.0d0) then
            write(6,'(a,i0,a,g11.4)')
     &         '           stop at ', itmx,
     &         ' iterations or sd/|sigma| < ', epsrel
            flush(6)
         endif
         do i=start_it+1,itmx
            call vegas_m(4,fxn,avgi,sd,chisq)
            call histo_end_iteration
            call checkpoint_write_state
            if (vs_done(epsrel,avgi,sd,chisq,i,itmx)) then
               converged=.true.
               return
            endif
         enddo
         call vs_ceiling(epsrel,avgi,sd,chisq,itmx)
         return
      endif
c
c  --- adaptation: grid refined, results discarded ---
      ncall=npt_ad
      itmx=its_ad
      plot=.false.
c
c  Every progress write below is flushed on the spot. A campaign job runs
c  with stdout redirected to a log (tools/run_seeds.sh), and a redirected
c  unit is block-buffered: without the flush the whole iteration table
c  appears only when the run ends, which is exactly when it stops being
c  useful. One flush per iteration costs nothing next to npt points.
c
      write(6,'(/a)') ' [ VEGAS ] adaptation pass -- grid refining'
      flush(6)
      call vegas_m(1,fxn,avgi,sd,chisq)
      do i=1,itmx
         call vegas_m(3,fxn,avgi,sd,chisq)
      enddo
c
c  --- production: grid frozen, histograms filled ---
      ncall=npt_pr
      itmx=its_pr
      plot=.true.
      write(6,'(/a)') ' [ VEGAS ] production pass -- grid frozen'
      if (epsrel.gt.0d0) then
         write(6,'(a,i0,a,g11.4)')
     &      '           stop at ', itmx,
     &      ' iterations or sd/|sigma| < ', epsrel
      endif
      flush(6)
      call vegas_m(2,fxn,avgi,sd,chisq)
      call histo_begin_production
      do i=1,itmx
         call vegas_m(4,fxn,avgi,sd,chisq)
         call histo_end_iteration
         call checkpoint_write_state
         if (vs_done(epsrel,avgi,sd,chisq,i,itmx)) then
            converged=.true.
            return
         endif
      enddo
      call vs_ceiling(epsrel,avgi,sd,chisq,itmx)
      return
      end

c#####################################################################
c  The production-pass stopping test.
c
c  A drop in SD alone is not convergence: with few iterations SD can
c  fluctuate low while the iteration-to-iteration spread still says the
c  grid does not describe the integrand. So a second condition is
c  required. CHI2A, which vegas_m already reports PER DEGREE OF
c  FREEDOM, must be sane: 1 is ideal and a few is tolerable, so
c  CHI2MAX=5 is loose enough not to fire on an honest run and tight
c  enough to catch a pass whose iterations genuinely disagree.
c
c  IT IS THE CUMULATIVE ITERATION COUNT IN BOTH LOOPS, the resumed one
c  included. start_it was read from '# iterations_done =', which
c  checkpoint_snapshot writes from VEGAS's own cumulative IT;
c  vegas_put_state restores that IT and analysis_read_state restores the
c  histogram nit. So the loop variable i is simultaneously the VEGAS IT,
c  the histogram nit, and the total -- no separate counter is needed.
c
c  Nothing is tested below MINIT iterations: histo1d_error returns 0 for
c  nit < 2, and CHI2A is identically 0 at IT=1 and carries a single very
c  noisy degree of freedom at IT=2.
c#####################################################################
      logical function vs_done(epsrel,avgi,sd,chisq,it,itmx)
      implicit double precision (a-h,o-z)
      integer it,itmx
      parameter (CHI2MAX=5d0)
      parameter (MINIT=3)
      vs_done=.false.
      if (epsrel.le.0d0) return
      if (it.lt.MINIT) return
      if (avgi.eq.0d0) return
      rel=sd/dabs(avgi)
      if (rel.ge.epsrel) return
      if (chisq.gt.CHI2MAX) return
      vs_done=.true.
      write(6,'(/a,i0,a,i0)')
     &   ' [ VEGAS ] CONVERGED -- stopped at iteration ', it,
     &   ' of a ceiling of ', itmx
      write(6,'(a,g11.4,a,g11.4)')
     &   '           sd/|sigma| = ', rel, '  <  epsrel = ', epsrel
      write(6,'(a,g11.4,a,g11.4)')
     &   '           chi2/dof   = ', chisq, '  <= ', CHI2MAX
      write(6,'(a)')
     &   '           checkpoint KEPT -- lower epsrel and rerun to'
      write(6,'(a)')
     &   '           continue from here.'
      flush(6)
      end
c#####################################################################
c  Reported when the ceiling was reached without meeting epsrel. Silent
c  when the feature is off, so a card with no epsrel key sees no new
c  output at all.
c#####################################################################
      subroutine vs_ceiling(epsrel,avgi,sd,chisq,itmx)
      implicit double precision (a-h,o-z)
      integer itmx
      if (epsrel.le.0d0) return
      write(6,'(/a,i0,a)')
     &   ' [ VEGAS ] CEILING reached -- ', itmx,
     &   ' iterations, epsrel NOT met'
      if (avgi.ne.0d0) then
         write(6,'(a,g11.4,a,g11.4)')
     &      '           sd/|sigma| = ', sd/dabs(avgi),
     &      '  vs epsrel = ', epsrel
      endif
      write(6,'(a,g11.4)')
     &   '           chi2/dof   = ', chisq
      write(6,'(a)')
     &   '           raise it_prod and rerun with resume = yes to'
      write(6,'(a)')
     &   '           add iterations.'
      flush(6)
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
C     MXDIM is the largest NDIM this routine can take, and it is checked
C     rather than trusted: the arrays below are FIXED, VEGAS_M never sees
C     the caller's NDIM until it is already indexing them, and overrunning
C     them corrupts the grid silently. Upstream this was 10, which is
C     exactly what  2 (x1,x2) + 8 (PSGen) + 1 (channel) = 11  needs one
C     more than. The 50 is NDMX, bins per axis, and is a different number.
      PARAMETER (MXDIM=16)
      DIMENSION XI(50,MXDIM),D(50,MXDIM),DI(50,MXDIM),XIN(50),R(50),
     1     DT(MXDIM),X(MXDIM),KG(MXDIM),IA(MXDIM),QRAN(MXDIM)
C     Dummies of the VEGAS_GET_STATE/VEGAS_PUT_STATE entries at the foot
C     of this routine. They belong here rather than beside the ENTRY
C     lines because a declaration may not follow executable statements.
      DIMENSION XI_O(50,MXDIM),XI_I(50,MXDIM)
      EXTERNAL FXN
      SAVE
      LOGICAL FIRST_CLOCK
      INTEGER*8 CLOCK_START,CLOCK_NOW,CLOCK_RATE,CLOCK_MAX
      CHARACTER*20 ELAPSED_STR
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1D0/,MDS/1/,FIRST_CLOCK/.TRUE./
      IF(NDIM.GT.MXDIM.OR.NDIM.LT.1)THEN
         WRITE(6,*)'ERROR vegas_m: NDIM=',NDIM,' outside 1..',MXDIM
         STOP 1
      ENDIF
      IF(FIRST_CLOCK)THEN
         CALL SYSTEM_CLOCK(CLOCK_START,CLOCK_RATE,CLOCK_MAX)
         FIRST_CLOCK=.FALSE.
      ENDIF
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
         IF(NPRN.NE.0) THEN
            WRITE(6,200) NDIM,CALLS,IT,ITMX,MDS,ND
            FLUSH(6)
         ENDIF
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
            CALL SYSTEM_CLOCK(CLOCK_NOW)
            IF(CLOCK_NOW.LT.CLOCK_START) CLOCK_NOW=CLOCK_NOW+CLOCK_MAX
            ELAPSED_SEC=DBLE(CLOCK_NOW-CLOCK_START)/DBLE(CLOCK_RATE)
            NDAYS=INT(ELAPSED_SEC/86400D0)
            NHOURS=INT((ELAPSED_SEC-NDAYS*86400D0)/3600D0)
            NMINS=INT((ELAPSED_SEC-NDAYS*86400D0-NHOURS*3600D0)/60D0)
            IF(NDAYS.GT.0)THEN
               WRITE(ELAPSED_STR,'(A,I0,A,I0,A,I0,A)')
     &            '[',NDAYS,'d ',NHOURS,'h ',NMINS,'m]'
            ELSEIF(NHOURS.GT.0)THEN
               WRITE(ELAPSED_STR,'(A,I0,A,I0,A)')
     &            '[',NHOURS,'h ',NMINS,'m]'
            ELSE
               WRITE(ELAPSED_STR,'(A,I0,A)') '[',NMINS,'m]'
            ENDIF
            WRITE(6,201) IT,TI,TSI,AVGI,SD,CHI2A,ELAPSED_STR
c            WRITE(ilabel,201) IT,TI,TSI,AVGI,SD,CHI2A,ELAPSED_STR
            FLUSH(6)

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
201   FORMAT(i4,4x,g15.7,g11.4,g15.7,g11.4,g11.4,2x,a)
202   FORMAT(' DATA FOR AXIS',I2 / ' ',6X,'X',7X,'  DELT I  ',
     1    2X,' CONV''CE  ',11X,'X',7X,'  DELT I  ',2X,' CONV''CE  '
     2   ,11X,'X',7X,'  DELT I  ',2X,' CONV''CE  ' /
     2    (' ',3G12.4,5X,3G12.4,5X,3G12.4))
      RETURN
C
C=====================================================================
C  CHECKPOINT HOOKS.
C
C  Everything above lives in this routine's SAVE, which is exactly why
C  it is otherwise unreachable: there is no common block and no argument
C  through which a caller could see the grid. These two ENTRY points
C  open that state up for a checkpoint and nothing else, using the same
C  trick BRM48I/BRM48O already use on the RNG a few hundred lines down.
C
C  THE DERIVED BLOCK (NG,MDS,NPG,CALLS,DXG,DV2G,XND,NDM,XJAC) IS
C  CARRIED WHOLE, and a resume must never replay ISTAT=1 or ISTAT=2 to
C  rebuild it. The ISTAT=2 branch recomputes ND from NCALL and rescales
C  XI from NDO bins to the new ND -- but never updates NDO, so ND.NE.NDO
C  permanently once production starts. Replaying it against a grid that
C  has ALREADY been rescaled would fire that rescale a second time and
C  silently deform the grid.
C
C  AVGI/SD/CHI2A are absent on purpose: they are rebuilt in full from
C  SI/SI2/SWGT/SCHI/IT at the end of the very next iteration. TI/TSI/
C  FB/F2B likewise -- per-iteration scratch, zeroed at the top of every
C  ISTAT=3/4 call, never live at the boundary a checkpoint is taken on.
C=====================================================================
      ENTRY VEGAS_GET_STATE(ND_O,NDO_O,NG_O,MDS_O,NPG_O,NDM_O,IT_O,
     &  CALLS_O,DXG_O,DV2G_O,XND_O,XJAC_O,SI_O,SI2_O,SWGT_O,SCHI_O,
     &  XI_O)
      ND_O   = ND
      NDO_O  = NDO
      NG_O   = NG
      MDS_O  = MDS
      NPG_O  = NPG
      NDM_O  = NDM
      IT_O   = IT
      CALLS_O= CALLS
      DXG_O  = DXG
      DV2G_O = DV2G
      XND_O  = XND
      XJAC_O = XJAC
      SI_O   = SI
      SI2_O  = SI2
      SWGT_O = SWGT
      SCHI_O = SCHI
      DO J=1,MXDIM
         DO I=1,50
            XI_O(I,J)=XI(I,J)
         ENDDO
      ENDDO
      RETURN
C
      ENTRY VEGAS_PUT_STATE(ND_I,NDO_I,NG_I,MDS_I,NPG_I,NDM_I,IT_I,
     &  CALLS_I,DXG_I,DV2G_I,XND_I,XJAC_I,SI_I,SI2_I,SWGT_I,SCHI_I,
     &  XI_I)
      ND   = ND_I
      NDO  = NDO_I
      NG   = NG_I
      MDS  = MDS_I
      NPG  = NPG_I
      NDM  = NDM_I
      IT   = IT_I
      CALLS= CALLS_I
      DXG  = DXG_I
      DV2G = DV2G_I
      XND  = XND_I
      XJAC = XJAC_I
      SI   = SI_I
      SI2  = SI2_I
      SWGT = SWGT_I
      SCHI = SCHI_I
      DO J=1,MXDIM
         DO I=1,50
            XI(I,J)=XI_I(I,J)
         ENDDO
      ENDDO
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
