      SUBROUTINE TRANSF(IGRAPH,DEPTH_FLAG)
C
C  THIS SUBROUTINE EXPLICITLY SOLVES THE TRANSFER EQUATION
C  FOR A SET OF NODES ON THE STAR DISK. THE RESULTS ARE:
C  AN ARRAY TABLE(WAVELENGTH) WITH SPECIFIC INTENSITIES
C  (LINE OPACITY INCLUDED) AND FC* WITH CONTINUUM INTENSITIES
C  AT BOTH ENDS OF SPECTRAL INTERVAL. THE RESULTS ARE
C  WRITTEN TO THE FILE #11, AS WELL AS THE INFORMATION ABOUT
C  THE NUMBER OF WAVELENGTHS, THE NUMBER OF NODES ON THE DISK,
C  MODEL TEMPERATURE AND GRAVITY, THE ABUNDANCE AND
C  THE WAVELENGTH RANGE.
C
C  Author: N.Piskunov
C
C  LAST UPDATE: September 13, 1993.
C
      INCLUDE 'SIZES.SYN'
      INCLUDE 'COMMONS.SYN'
      DOUBLE PRECISION EPS1
      PARAMETER (EPS1=1.D-05,EPS2=0.0002,NSTOT=6)
C
C FLUX_SCALE = exp(30.)*1.E-8
C
      PARAMETER (FLXSCL=1.0686475E5)
      LOGICAL DEPTH_FLAG
      REAL TABLE(NWSIZE),CONTRIB(NWSIZE)
      DOUBLE PRECISION WL(NWSIZE)
      DOUBLE PRECISION TAULIN(NSTOT+1),WGHTS(NSTOT)
      SAVE TAULIN,WGHTS
      DATA TAULIN/
     * 0.00000, 0.16134, 0.94600, 2.52410, 5.00660, 8.88090, 14.75900/
      DATA WGHTS/
     * 0.27191, 0.18626, 0.38898E-1, 0.28674E-2, 0.60078E-4, 0.17675E-6/
C
C  INTEGRATE TRANSFER EQUATION FOR FLUX
C
      CALL RKINTS(TAULIN,WGHTS,NSTOT,EPS1,EPS2,FCBLUE,FCRED,
     *            TABLE,CONTRIB,NWL,WL,IGRAPH)
C
      DO 1 LINE=1,NLINES
      IF(GAMRAD(LINE).GT.0) THEN
        GRLG10=LOG10(GAMRAD(LINE))
      ELSE
        GRLG10=0.0
      END IF
      IF(GAMQST(LINE).GT.0.AND.ELEMC(LINE).NE.'H') THEN
        GSLG10=LOG10(GAMQST(LINE))
      ELSE
        GSLG10=GAMQST(LINE)
      END IF
      IF(GAMVW(LINE).GT.0.AND.GAMVW(LINE).LT.10.0.AND.
     *   ELEMC(LINE).NE.'H') THEN
        GWLG10=LOG10(GAMVW(LINE))
      ELSE
        GWLG10=GAMVW(LINE)
      END IF
      CALL NOBLNK(SREFGF(LINE),J1,J2)
      WRITE(11,200) WLCENT(LINE),ELEMC(LINE),EXCIT(LINE),
     *              VTURB(LINE),LOG10(GF(LINE)),GRLG10,
     *              GSLG10,GWLG10,XLANDE(LINE),RLINE(LINE),
     *              SREFGF(LINE)(J1:J2)
 200  FORMAT(F10.4,',''',A4,''',',F8.4,',',F4.1,',',F7.3,',',
     *       F8.3,',',F8.3,',',F8.3,',',F8.3,',',F7.3,
     *          ', ''',A,'''')
   1  CONTINUE
C
C  INTEGRATE TRANSFER EQUATION FOR SPECIFIC INTENSITY AND EVALUATE
C  LIMB DARKENING FOR CONTINUUM IN A FORM I=I0*(1-U+U*MU)
C
      CALL LIMBD(EPS1,U,FCA,FCB,FCC)
C
C  WRITE RESULTS TO THE FILE
C
      WRITE(11,'(I5,F8.4,'' - number of wavelengths, limb darkening'')')
     *           NWL,U
      DO 2 IWL=1,NWL
      IF(DEPTH_FLAG) THEN
        WRITE(11,'(F10.4,F10.7,E13.6)') WL(IWL),TABLE(IWL),CONTRIB(IWL)
      ELSE
        WRITE(11,'(F10.4,F10.7)') WL(IWL),TABLE(IWL)
      END IF
   2  CONTINUE
c      WRITE(11,'(8F10.4)') (WL(IWL),IWL=1,NWL)
c      WRITE(11,'(8F10.7)') (TABLE(IWL),IWL=1,NWL)
C
C  WRITE MODEL ATMOSPHERE PARAMETERS
C
      WRITE(11,*) 'Temperature=',TEFF,',  Gravity=',GRAV
      WRITE(11,'('' FCblue='',E12.6,'' FCred='',E12.6,
     *           '' FC(0)='',E12.6)')
     *              FCBLUE*FLXSCL,FCRED*FLXSCL,FCA*FLXSCL
C
      RETURN
      END

      SUBROUTINE RKINTS(TAULIN,WGHTS,NSTOT,EPS1,EPS2,
     *                  FCBLUE,FCRED,TABLE,CONTRIB,NWL,WL,IGRAPH)
C
C  THIS SUBROUTINE CALLS SUBROUTINE FCINTG TO INTEGRATE THE EMMERGING
C  SPECIFIC INTENSITIES FOR CONTINUUM AT THE EDGES OF SPECTRAL
C  INTERVAL (RETURNED AS "FC*") AND SUBROUTINE TBINTG FOR THE LINE
C  (RETURNED AS "TABLE").
C  DELWL is the minimum stepsize in wavelength.
C
C  Author: N.Piskunov
C
C  LAST UPDATE: August 2, 1996.
C
      INCLUDE 'SIZES.SYN'
      INCLUDE 'COMMONS.SYN'
      DOUBLE PRECISION TAULIN(NSTOT),WGHTS(NSTOT-1)
      REAL TABLE(NWSIZE),CONTRIB(NWSIZE),RINTNS(LINSIZ)
      DOUBLE PRECISION WW,WL(NWSIZE),EPS1,DELWL
      PARAMETER (EPS3=5.,DELWL=1.D-3)
C
C  CALCULATE CONTINUUM FLUX FOR BOTH ENDS OF THE INTERVAL
C
      CALL FCINTG(TAULIN,WGHTS,NSTOT,FCBLUE,FCRED,EPS1)
      FCSLOP=(FCRED-FCBLUE)/(WLAST-WFIRST)
c      write(*,*) 'FCblue=',FCBLUE,',   FCred=',FCRED
c      write(*,*)
C
C  NOW WE SHELL INCLUDE LINE OPACITY INTO CONSIDERATION AND
C  START A LOOP OVER WAVELENGTHS IN A GIVEN SPECTRAL INTERVAL
C
C  FIRST WE CALCULATE FLUX AT THE BLUE END OF SPECTRAL INTERVAL
C
      WL(1)=WFIRST
      CALL TBINTG(WFIRST,TAULIN,WGHTS,NSTOT,TABLE(1),CONTRIB(1),EPS1)
      TABLE(1)=TABLE(1)/FCBLUE
C
C  ADD ONE POINT AT EACH LINE CENTER AND ONE POINT BETWEEN EACH PAIR OF LINES
C
      IWL=1
      DO 3 LINE=1,NLINES
      WW=WLCENT(LINE)
      IF(WW.GT.WFIRST.AND.WW.LT.WLAST.AND.WW-WL(IWL).GT.DELWL) THEN
        IWL=IWL+1
        WL(IWL)=WW
        CALL TBINTG(WL(IWL),TAULIN,WGHTS,NSTOT,TABLE(IWL),
     *              CONTRIB(IWL),EPS1)
        TABLE(IWL)=TABLE(IWL)/(FCSLOP*(WL(IWL)-WFIRST)+FCBLUE)
        IF(LINE.GT.1.) THEN
          DO 1 LIN=1,LINE
          IF(.NOT.MARK(LIN).AND.WRIGHT(LIN).GT.WL(IWL).AND.
     *       ALMAX(LIN)*100.LT.EPS2) WRIGHT(LIN)=WL(IWL)
   1      CONTINUE
        END IF
        IF(TABLE(IWL).GT.1.-EPS2) MARK(LINE)=.TRUE.
        IF(LINE.LT.NLINES) THEN
          DO 2 LIN=LINE+1,NLINES
          IF(.NOT.MARK(LIN).AND.ALMAX(LIN)*100.LT.EPS2)
     *       WLEFT(LIN)=WL(IWL)
   2      CONTINUE
        END IF
      END IF
      IF(WW.GT.WFIRST.AND.WW.LT.WLAST) RINTNS(LINE)=TABLE(IWL)
C
C  Central line strength
C
      CALL LINSTR(WLCENT(LINE),TAULIN,WGHTS,NSTOT,
     *            RLINE(LINE),EPS1,LINE)
      RLINE(LINE)=1.-RLINE(LINE)/(FCSLOP*(WLCENT(LINE)-WFIRST)+FCBLUE)
   3  CONTINUE
C
C  AND AT LAST DO THE SAME WITH THE RED END OF SPECTRAL INTERVAL
C
      IF(WLAST-WL(IWL).GT.DELWL) THEN
        IWL=IWL+1
        WL(IWL)=WLAST
        CALL TBINTG(WL(IWL),TAULIN,WGHTS,NSTOT,TABLE(IWL),
     *              CONTRIB(IWL),EPS1)
        TABLE(IWL)=TABLE(IWL)/FCRED
      END IF
      NWL=IWL
C
C  AND NOW ADJUST STEP SIZE OF ABS(TABLE(IWL)-TABLE(IWL-1)) IS TOO BIG
C
      IPERCO=0
      IPERC =0
      IWL=2

   4  IF(NWL.GE.NWSIZE) THEN
        WRITE(*,*) 'NOT ENOUGH NWSIZE TO OBTAIN REQUIRED ACCURACY'
        RETURN
      END IF
C
C  Linear estimate in the intermediate point
C
      FCL=0.5*(TABLE(IWL)+TABLE(IWL-1))
      DO 5 IIWL=NWL,IWL,-1
      WL(IIWL+1)=WL(IIWL)
      TABLE  (IIWL+1)=TABLE  (IIWL)
      CONTRIB(IIWL+1)=CONTRIB(IIWL)
   5  CONTINUE
      NWL=NWL+1
      WL(IWL)=(WL(IWL-1)+WL(IWL+1))*0.5D0
C
C  Calculate true value
C
      CALL TBINTG(WL(IWL),TAULIN,WGHTS,NSTOT,TABLE(IWL),
     *            CONTRIB(IWL),EPS1)
      TABLE(IWL)=TABLE(IWL)/(FCSLOP*(WL(IWL)-WFIRST)+FCBLUE)
      IF(WL(IWL)-WL(IWL-1).LT.DELWL.OR.
     *   ABS(TABLE(IWL)-FCL).LT.EPS2) THEN
        IPERC=INT(100.*(WL(IWL)-WFIRST)/(WLAST-WFIRST))
C
C  Now we will move right of the WL(IWL) and will never come back, mark
C  permanently all weak lines left of this wavelength. Unmark all
C  temporary marked lines.
C
        DO 6 LINE=1,NLINES
        IF(.NOT.MARK(LINE).AND.WRIGHT(LINE).GT.WL(IWL).AND.
     *     WLCENT(LINE).LE.WL(IWL).AND.ALMAX(LINE)*100.LT.EPS2) then
          if(IELEM(line).NE.1) WRIGHT(LINE)=WL(IWL)
        END IF
   6    CONTINUE
        IWL=IWL+2
      ELSE
        if(IPERC.gt.IPERCO) then
          IPERCO=IPERC
          write(*,200) NWL,WL(IWL),IPERC
        end if
C
C At this point we will add more points to the left, so we can
C ignore all weak lines to the right of this wavelength.
C
        DO 7 LINE=1,NLINES
        IF(.NOT.MARK(LINE).AND.WLEFT(LINE).LT.WL(IWL).AND.
     *     WLCENT(LINE).GT.WL(IWL).AND.
     *     ALMAX(LINE)*100.LT.EPS2) WLEFT(LINE)=WL(IWL)
   7    CONTINUE
      END IF
      IF(IWL.LE.NWL) GO TO 4
C
 200  format(' NWL=',I5,',  WL=',F10.4,',  Done=',I4,'%')
      RETURN
      END


      SUBROUTINE TBINTG(WL,TAULIN,WGHTS,NSTOT,FLUX,CONTRIB,EPS)
C
C  THIS SUBROUTINE INTEGRATE THE SPECIFIC INTENSITY AT WAVELENGTH "WL".
C  THE RESULT IS RETURNED AS "FLUX".
C  "TAULIN" IS THE STANDARD OPTICAL DEPTH SCALE AND "BPLANK" & "TAUNEW"
C  ARE WORKING ARRAYS TO STORE INTERMEDIATE VALUES OF PLANK FUNCTION
C  AND LINE OPTICAL DEPTH.
C
C  Author: N.Piskunov
C
C  LAST UPDATE: March 15, 1991.
C
      INCLUDE 'SIZES.SYN'
      INCLUDE 'COMMONS.SYN'
      LOGICAL FIRSTS
      DOUBLE PRECISION TAULIN(NSTOT),WGHTS(NSTOT-1),BPLANK
      DOUBLE PRECISION WL,FL,CONWL5,RATIO,TAUNEW,XBEG,XEND,EPS,FRDP
      EXTERNAL RATIO
C
C  INTEGRATE THE STANDARD WL OPTICAL DEPTH NODES FOR
C  CORRESPONDING TAULIN AT WL
C
      CONWL5=EXP(50.7649141-5.*LOG(WL))
      HNUK=1.43868E+08/WL
      TAUNEW=RHOX(1)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  here is a good place to write out model stuff
c      write(*,*) NRHOX
c      WL=WLCENT(2)
c      do 99 im=1,nrhox
c 99   W=RATIO(RHOX(im),WL,0,TEMPER)
c      stop
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  THE LOOP THROUGH THE ATMOSPHERE USING THE LOG DEPTH SCALE
C  STARTING AT THE SURFACE
C
      FL=0.D0
      FRDP=0.D0
      FIRSTS=.TRUE.
      DO 1 IM=2,NSTOT
      XBEG=TAULIN(IM-1)
      XEND=TAULIN(IM)
      CALL RKINT(XBEG,XEND,EPS,TAUNEW,RHOX(2)-RHOX(1),WL,0,
     *           FIRSTS,TEMPER,RATIO)
      FIRSTS=.FALSE.
C
C  PLANK FUNCTION AGAIN
C
      BPLANK=CONWL5/(EXP(HNUK/TEMPER)-1.)
      FL=FL+BPLANK*WGHTS(IM-1)
      FRDP=FRDP+TAUNEW*BPLANK*WGHTS(IM-1)
   1  CONTINUE
      FLUX=FL
      CONTRIB=FRDP/FL
C
      RETURN
      END

      SUBROUTINE LINSTR(WL,TAULIN,WGHTS,NSTOT,FLUX,EPS,LINE)
C
C  THIS SUBROUTINE INTEGRATE THE SPECIFIC INTENSITY AT WAVELENGTH "WL".
C  THE RESULT IS RETURNED AS "FLUX".
C  "TAULIN" IS THE STANDARD OPTICAL DEPTH SCALE AND "BPLANK" & "TAUNEW"
C  ARE WORKING ARRAYS TO STORE INTERMEDIATE VALUES OF PLANK FUNCTION
C  AND LINE OPTICAL DEPTH.
C
C  Author: N.Piskunov
C
C  LAST UPDATE: August 1, 1997.
C
      INCLUDE 'SIZES.SYN'
      INCLUDE 'COMMONS.SYN'
      LOGICAL FIRSTS
      DOUBLE PRECISION TAULIN(NSTOT),WGHTS(NSTOT-1),BPLANK
      DOUBLE PRECISION WL,FL,CONWL5,RATIO,TAUNEW,XBEG,XEND,EPS
      EXTERNAL RATIO
C
C  INTEGRATE THE STANDARD WL OPTICAL DEPTH NODES FOR
C  CORRESPONDING TAULIN AT WL
C
      CONWL5=EXP(50.7649141-5.*LOG(WL))
      HNUK=1.43868E+08/WL
      TAUNEW=RHOX(1)
C
C  THE LOOP THROUGH THE ATMOSPHERE USING THE LOG DEPTH SCALE
C  STARTING AT THE SURFACE
C
      FL=0.D0
      FIRSTS=.TRUE.
      DO 1 IM=2,NSTOT
      XBEG=TAULIN(IM-1)
      XEND=TAULIN(IM)
      CALL RKINT(XBEG,XEND,EPS,TAUNEW,RHOX(2)-RHOX(1),WL,LINE,
     *           FIRSTS,TEMPER,RATIO)
      FIRSTS=.FALSE.
C
C  PLANK FUNCTION AGAIN
C
      BPLANK=CONWL5/(EXP(HNUK/TEMPER)-1.)
      FL=FL+BPLANK*WGHTS(IM-1)
   1  CONTINUE
      FLUX=FL
C
      RETURN
      END

      SUBROUTINE FCINTG(TAULIN,WGHTS,NSTOT,FCBLUE,FCRED,EPS)
C
C  THIS SUBROUTINE INTEGRATE CONTINUUM FLUX FOR BOTH ENDS OF
C  THE SPECTRAL INTERVAL. THE RESULTS ARE RETURNED AS "FCBLUE" AND "FCRED".
C  "TAULIN" IS THE STANDARD OPTICAL DEPTH SCALE.
C
C  Author: N.Piskunov
C
C  LAST UPDATE: March 15, 1991.
C
      INCLUDE 'SIZES.SYN'
      INCLUDE 'COMMONS.SYN'
      LOGICAL FIRSTS
      DOUBLE PRECISION TAULIN(NSTOT+1),WGHTS(NSTOT),BPLANK
      DOUBLE PRECISION WLCONT,FC,CONWL5,RATIO,TAUNEW,XBEG,XEND,EPS
      EXTERNAL RATIO
C
C  INTEGRATE THE STANDARD WL OPTICAL DEPTH NODES FOR CONTINUUM
C  CORRESPONDING TO TAULIN AT WLCONT.
C
      DO 2 IWL=1,2
      IF(IWL.EQ.1) WLCONT=WFIRST
      IF(IWL.EQ.2) WLCONT=WLAST
      CONWL5=EXP(50.7649141-5.0*LOG(WLCONT))
      HNUK=1.43868E+08/WLCONT
      TAUNEW=RHOX(1)
      FC=0.D0
      FIRSTS=.TRUE.
C
C  THE LOOP THROUGH THE ATMOSPHERE USING THE LOG DEPTH SCALE
C  STARTING AT THE SURFACE
C
      DO 1 IM=2,NSTOT
      XBEG=TAULIN(IM-1)
      XEND=TAULIN(IM)
      CALL RKINT(XBEG,XEND,EPS,TAUNEW,RHOX(2)-RHOX(1),WLCONT,-IWL,
     *           FIRSTS,TEMPER,RATIO)
      FIRSTS=.FALSE.
C
C  PLANK FUNCTION AGAIN
C
      BPLANK=CONWL5/(EXP(HNUK/TEMPER)-1.)
      FC=FC+BPLANK*WGHTS(IM-1)
   1  CONTINUE
      IF(IWL.EQ.1) FCBLUE=FC
      IF(IWL.EQ.2) FCRED =FC
   2  CONTINUE
C
      RETURN
      END


      DOUBLE PRECISION FUNCTION RATIO(TAU1,WAVE,ICODE,TEMPER)
C
C  THIS FUNCTION CALCULATES OPACITY OR OPACITY RATIO (OPACWL/OPACSTD)
C  PER GRAMM OF STELLAR MATER (cm^2/gm) PER ANGSTROEM AT DEPTH #IM
C  OF THE STANDARD MODEL DEPTH SCALE. WAVELENGTH IS TAKEN EITHER FROM
C  WAVE (ICODE=0) OR FROM EDGES OF SPECTRAL INTERVAL (ICODE=1,2).
C
C  Author: N.Piskunov
C
C  LAST MAJOR UPDATES: August 5, 1996.
C  February 12, 1999 : Temperature dependent van der Waals if
C                      ALPHA and SIGMA are available and 
C                      reduced mass of perturbers by Paul Barklem 
C  July 3, 2000:       New inteface with Hydrogen line opacity routine
C
      INCLUDE 'SIZES.SYN'
      INCLUDE 'COMMONS.SYN'
C
C                   pi*e^2
C  Line opacity is: ------ * gf * N_absorb * STIM * f(nu-nu0)
C                    m*c
C
C  where the line profile f(nu) is assumed to be normalized so that:
C
C  \integ f(nu-nu0) d nu = 1
C
C  This is true for Voigt, Hydrogen and (I hope) Fano profiles.
C                                                     1
C  E.g., in case of Voigt profile f(nu-nu0)= -------------------- * H(a,v)
C                                            sqrt(pi)*del_nu_Dopp
C  where del_nu_Dopp = DNDOPL  is in Hz,
C
C  where H(a,v) is the Voigt function with normalization:
C  \integ H(a,v) d v = sqrt(pi)
C
C  Two Hydrogen line profiles are computed externally by Kurucz
C  approximation (HLINOP) or by interpolation in Stehle's tables (HTABLE)
C  and are area normalized!
C
C  Therefore the normalization factor Z=PI*e^2/(m*c) with speed
C  of light in cm/s. The net result is that Z is in cm^2/s !!!
C
C  Other constants: K  - Boltzmann's constant J/K,
C                   M0 - unit atomic mass kg (Carbon 12 scale),
C                   A0 - Bohr radius m
C
      REAL      GX,X,SIGMA,ALPHA,K,M0,A0,GAMMAF,VBAR
      PARAMETER (C=2.997925E+18,PI=3.14159265,C4PI=C*4.*PI)
      PARAMETER (SQRTPI=1.77245385,K=1.380658E-23,M0=1.660540E-27)
      PARAMETER (A0=5.29177249E-11,Z=0.026540045)
c      PARAMETER (VWFUDGE=2.5)
      PARAMETER (VWFUDGE=1.0)

      DOUBLE PRECISION WAVE,WL1,TAU,TAU1,DSMTAU
      REAL*8 WAVE2,VACAIR,AIRVAC
      REAL YABUND(LINSIZ),XMASS(LINSIZ),EXCUP(LINSIZ),VTURB2(LINSIZ),
     *     GR(LINSIZ),ENU4(LINSIZ),ENL4(LINSIZ)
      LOGICAL FIRSTL(LINSIZ)
      INTEGER IDHEL(LINSIZ)
      SAVE WL1,FIRSTL,IDHEL,YABUND,XMASS,EXCUP,VTURB2,GR,ENU4,ENL4,
     *     HNUXXX,DDWL

      DATA WL1/-1.D0/,FIRSTL/LINSIZ*.TRUE./
C
      TAU=MAX(MIN(TAU1,RHOX(NRHOX)),RHOX(1))
      IM=2
   1  IF(RHOX(IM).LT.TAU.AND.IM.LT.NRHOX) THEN
        IM=IM+1
        GO TO 1
      END IF
      TEMPER=SPLINT(IM-1,IM,RHOX,T,     SP_T,  NRHOX,TAU)
      OPCONB=SPLINT(IM-1,IM,RHOX,COPBLU,SP_CBL,NRHOX,TAU)
      OPCONR=SPLINT(IM-1,IM,RHOX,COPRED,SP_CRD,NRHOX,TAU)
      IF(MOTYPE.EQ.0) THEN
        RATIO=1.D0
      ELSE
        RATIO=SPLINT(IM-1,IM,RHOX,COPSTD,SP_CST,NRHOX,TAU)
      END IF
C
C  SELECT ICODE
C
      IF(ICODE.EQ.-1) THEN
        RATIO=RATIO/OPCONB
      ELSE IF(ICODE.EQ.-2) THEN
        RATIO=RATIO/OPCONR
      ELSE IF(ICODE.GE.0) THEN
        IF(ABS(WAVE-WL1).GT.1.0D-6) THEN
          WL1=WAVE
          HNUXXX=C*6.6256D-27/WAVE
          DDWL=(WAVE-WFIRST)/(WLAST-WFIRST)
C
C  Wavelength has changed. Clear opacity contribution array
C
          DO 2 LINE=1,NLINES
   2      ALMAX(LINE)=0.
        END IF
        OPCON=(OPCONR-OPCONB)*DDWL+OPCONB
        XNATOM=SPLINT(IM-1,IM,RHOX,XNA,SP_NA ,NRHOX,TAU)
        XNELEC=SPLINT(IM-1,IM,RHOX,XNE,SP_NE ,NRHOX,TAU)
        XXRHO =SPLINT(IM-1,IM,RHOX,RHO,SP_RHO,NRHOX,TAU)
C
C  PREPARE THINGS FOR IONIZATION FRACTIONS INTERPOLATION
C
        IF(TAU.GT.RHOX(IM).AND.IM.LT.NRHOX) THEN
          ISMTAU=NSMTAU*(IM-1)+INT((TAU-RHOX(IM))/HSMTAU(IM))+1
          DSMTAU=MOD(TAU-RHOX(IM),HSMTAU(IM))/HSMTAU(IM)
        ELSE
          ISMTAU=NSMTAU*(IM-2)+INT((TAU-RHOX(IM-1))/HSMTAU(IM-1))+1
          DSMTAU=MOD(TAU-RHOX(IM-1),HSMTAU(IM-1))/HSMTAU(IM-1)
        END IF
C
C  CALCULATE FRACTION OF H I AND HE I
C
        IF(ISMTAU.GT.(MOSIZE-1)*NSMTAU) THEN
          H1FRC =FRACT(1,3,(MOSIZE-1)*NSMTAU)
          HE1FRC=FRACT(2,4,(MOSIZE-1)*NSMTAU)
        ELSE
          H1FRC =(FRACT(1,3,ISMTAU+1)-FRACT(1,3,ISMTAU))*DSMTAU+
     +            FRACT(1,3,ISMTAU)
          HE1FRC=(FRACT(2,4,ISMTAU+1)-FRACT(2,4,ISMTAU))*DSMTAU+
     +            FRACT(2,4,ISMTAU)
        END IF
C
C  PREPARE SOME OTHER USEFULL THINGS
C
        XTK=TEMPER*1.38054E-16
        XSTIM=1.-EXP(-HNUXXX/XTK)
        VH=1.28466E+04*SQRT(TEMPER)
        GVWPRT=17.*(VH**0.6)*H1FRC
        TEMP6=(TEMPER/10000.)**(1./6.)*XNELEC
        TEMP3=(TEMPER/10000.)**0.3*(H1FRC+0.42*HE1FRC)
C
C  LOOP OVER SPECTRAL LINES
C
        ALINE=0.0
        IF(ICODE.EQ.0) THEN
          LINE1=1
          LINE2=NLINES
        ELSE
          LINE1=ICODE
          LINE2=ICODE
        END IF
        DO 3 LINE=LINE1,LINE2
        IF(MARK(LINE).OR.WAVE.LT.WLEFT(LINE).OR.
     *     WAVE.GT.WRIGHT(LINE)) GO TO 3
c        IF(IELEM(LINE).NE.1.AND..NOT.AUTOION(LINE).AND.
c     *     ABS(WAVE-WLCENT(LINE)).GT.15.0) GO TO 3
C
        WLC=WLCENT(LINE)
C
C  THE NUMBER OF ATOMS (NEW VERSION)
C
        IF(ISMTAU.LT.(MOSIZE-1)*NSMTAU) THEN
          FR=(FRACT(IELEM(LINE),ION(LINE),ISMTAU+1)-
     -        FRACT(IELEM(LINE),ION(LINE),ISMTAU))*DSMTAU+
     +        FRACT(IELEM(LINE),ION(LINE),ISMTAU)
        ELSE
          FR=FRACT(IELEM(LINE),ION(LINE),ISMTAU)
        END IF
        EFRACT=FR*EXP(-EXCIT(LINE)/(8.6171E-05*TEMPER))
C
C  Store some variables for current line that can be computed only once
C
        IF(FIRSTL(LINE)) THEN
          FIRSTL(LINE)=.FALSE.
          YABUND(LINE)=Z*ABUND(IELEM(LINE))*GF(LINE)
          XMASS(LINE)=1.66355E+24/C/C/AMASS(IELEM(LINE))
          EXCUP(LINE)=EXCIT(LINE)+1./(WLC*8065.47E-8)
C
C  VTURB is in km/s, 1.E13 converts C to km/s, so VTURB2 is dimensionless
C
          VTURB2(LINE)=1.E+26/C/C*VTURB(LINE)**2
C
C  Variables used in Unsold approximation for Stark and van der Waals damping
C
          IF(.NOT.AUTOION(LINE).AND.
     *       (GAMVW(LINE).LE.0.0.OR.GAMQST(LINE).LE.0.0)) THEN
            ENU4(LINE)=(ION(LINE)**2*13.598/
     /                 (POTION(IELEM(LINE),ION(LINE))-EXCUP(LINE)))**2
          END IF
          IF(.NOT.AUTOION(LINE).AND.GAMVW(LINE).LE.0.0) THEN
            ENL4(LINE)=(ION(LINE)**2*13.598/
     /                 (POTION(IELEM(LINE),ION(LINE))-EXCIT(LINE)))**2
          END IF
C
C  Radiative damping parameter: if not set, use classical value
C
          IF(GAMRAD(LINE).GT.0.0) THEN
            GR(LINE)=GAMRAD(LINE)
          ELSE
            GR(LINE)=0.222E+16/(WLC*WLC)
          END IF
C
C  Identify Helium lines included in Dimitrijevic & Sahal-Brechot table;
C  Stark damping for those will be computed in subroutine GAMHE
C
          IF(IELEM(LINE).EQ.2.AND..NOT.MARK(LINE)) THEN
            IDHEL(LINE)=0
            IF(INT(WLC).EQ.3819) IDHEL(LINE)= 1
            IF(INT(WLC).EQ.3867) IDHEL(LINE)= 2
            IF(INT(WLC).EQ.3871) IDHEL(LINE)= 3
            IF(INT(WLC).EQ.3888) IDHEL(LINE)= 4
            IF(INT(WLC).EQ.3926) IDHEL(LINE)= 5
            IF(INT(WLC).EQ.3964) IDHEL(LINE)= 6
            IF(INT(WLC).EQ.4009) IDHEL(LINE)= 7
            IF(INT(WLC).EQ.4120.OR.
     *         INT(WLC).EQ.4121) IDHEL(LINE)= 8
            IF(INT(WLC).EQ.4143) IDHEL(LINE)= 9
            IF(INT(WLC).EQ.4168.OR.
     *         INT(WLC).EQ.4169) IDHEL(LINE)=10
            IF(INT(WLC).EQ.4437) IDHEL(LINE)=11
            IF(INT(WLC).EQ.4471) IDHEL(LINE)=12
            IF(INT(WLC).EQ.4713) IDHEL(LINE)=13
            IF(INT(WLC).EQ.4921.OR.
     *         INT(WLC).EQ.4922) IDHEL(LINE)=14
            IF(INT(WLC).EQ.5015.OR.
     *         INT(WLC).EQ.5016) IDHEL(LINE)=15
            IF(INT(WLC).EQ.5047) IDHEL(LINE)=16
            IF(INT(WLC).EQ.5875) IDHEL(LINE)=17
            IF(INT(WLC).EQ.6678) IDHEL(LINE)=18
            IF(INT(WLC).EQ.4026) IDHEL(LINE)=19
            IF(INT(WLC).EQ.4387.OR.
     *         INT(WLC).EQ.4388) IDHEL(LINE)=20
          END IF
        END IF
C
C  This van der Waals part is written by Paul Barklem
C  Doppler broadening: DOPL is in fact delta_lambda/lambda
C  DLDOPL is delta_lambda in Angstroems
C  DNDOPL is delta_nu in Hz
C
        DOPL=SQRT(TEMPER*XMASS(LINE)+VTURB2(LINE))
        DLDOPL=WAVE*DOPL
        DNDOPL=DOPL*C/WAVE
C
C  If this happens to be a hydrogen line
C
        IF(IELEM(LINE).EQ.1) THEN
          NBLO=GAMQST(LINE)+0.1
          NBUP=GAMVW(LINE) +0.1

          HNORM=EFRACT*YABUND(LINE)*XSTIM*XNATOM/XXRHO
          CALL HLINPROF(WAVE,WLCENT(LINE),TEMPER,XNELEC,
     *                NBLO,NBUP,H1FRC,HE1FRC,DOPL,AHLIN1)
          AHLINE=HNORM*AHLIN1*WAVE*WAVE/C 
          IF(AHLIN1.LT.0.) THEN
            AHLIN2=HLINOP(WAVE,NBLO,NBUP,WLCENT(LINE),TEMPER,XNELEC,
     *                    H1FRC,HE1FRC,DOPL)
            AHLINE=HNORM*AHLIN2
          ENDIF
C
C  Store maximum contribution of hydrogen line to the total line opacity
C
          IF(AHLINE/OPCON.GT.ALMAX(LINE)) ALMAX(LINE)=AHLINE/OPCON
          ALINE=ALINE+AHLINE
          GO TO 3
        END IF
C
C  Quadratic Stark damping: if constant is set, follow D.Gray's
C  formula, otherwise use C.Cowley's modification of Unsold
C  approximation.
C
        IF(IELEM(LINE).NE.2.OR.IDHEL(LINE).EQ.0) THEN
          IF(GAMQST(LINE).GT.0.0.OR.AUTOION(LINE)) THEN
            GQST=GAMQST(LINE)*TEMP6
          ELSE
            IF(ION(LINE).EQ.1) THEN
              GQST=2.26E-07*ENU4(LINE)*XNELEC
            ELSE
              GQST=5.42E-07*ENU4(LINE)*XNELEC/(ION(LINE)+1)**2
            END IF
          END IF
        ELSE
C
C  For Helium lines from Dimitrijevich & Sahal-Brechot list compute
C  Stark broadening and wavelength shifts separately
C
          SHFT=0
          CALL GAMHE(IDHEL(LINE),TEMPER,XNELEC,XNATOM,GQST,SHFT)
          WLC=WLC+SHFT
        END IF
C
C  van der Waals damping parameter
C  
        IF(ANSTEE(LINE)) THEN
C
C  Compute the broadening by hydrogen from cross-section data which is
C  in m^2
C  Unpack the temperature dependent van der Waals parameters:
C  integer part is SIGMA and decimal part is ALPHA.
C
          SIGMA=INT(GAMVW(LINE))*A0*A0
          ALPHA=GAMVW(LINE)-INT(GAMVW(LINE))
C
C  Compute the Gamma function of X, this function is valid over the range 1<X<2
C
          X=2.-ALPHA*.5
          GX=X-1.0
          GAMMAF=1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*GX
     ;           )*GX)*GX)*GX)*GX
C
C  Compute the halfwidth per unit perturber density for vbar
C  
          GVW=(4./PI)**(ALPHA*0.5)*GAMMAF*1.E4*SIGMA
          VBAR=SQRT(8.*K*TEMPER/PI/M0*(1./1.008+1./AMASS(IELEM(LINE))))
          GVW=GVW*((VBAR/1.E4)**(1.-ALPHA))
C
C  Fullwidth given H1FRC perturbers per cm^3 with approximate HeI broadening
C  The factor of 2 comes from converting to the full width.
C  
          GVW=GVW*(H1FRC+0.42*HE1FRC)*1.E6*2.
        ELSE IF(.NOT.ANSTEE(LINE).AND.
     *          (GAMVW(LINE).GT.0.0.OR.AUTOION(LINE))) THEN
C
C  Input is log line width per unit density (rad/s cm^3)
C  
          GVW=GAMVW(LINE)*TEMP3
        ELSE
C
C  Input of zero means it is computed by Unsold theory
C  
          CW=1.61E-33*(ENU4(LINE)-ENL4(LINE))/(ION(LINE)*ION(LINE))
C
C  VWFUDGE is a fudge factor for Unsold formula found from fitting Solar
C  spectrum to be approx. 2.5.
C
          GVW=GVWPRT*CW**0.4*VWFUDGE
        END IF
cC 
cC  This is how it used to be
cC
c        IF(GAMVW(LINE).GT.0.0.OR.AUTOION(LINE)) THEN
c          GVW=GAMVW(LINE)*TEMP3
c        ELSE
c          CW=1.61E-33*(ENU4(LINE)-ENL4(LINE))/(ION(LINE)**2)
c          GVW=GVWPRT*CW**0.4*VWFUDGE
c        END IF
C
C  For autionizing line compute simple Fano profile
C  (U.Fano, Physical Review, 1961, Vol. 124, No. 6)
C  Fano Q is assumed to be passed via GAMQST(LINE)
C
c        IF(AUTOION(LINE)) THEN
c          DNU=2*C*(WAVE-WLC)/(WLC*WLC)/GR(LINE)
c          Q=GAMQST(LINE)
c          FANO=MAX(0.,(DNU*Q+1.)/(DNU*DNU+1.))*C*2/GR(LINE)/PI
c          ALINE1=EFRACT*YABUND(LINE)*XSTIM*XNATOM/XXRHO*FANO
c        ELSE
C
C  Compute total damping parameter. VWFUDGE is a fudge factor found from
C  fitting Solar spectrum to be approx. 2.5.
C
c          if(ABS(WAVE-WLCENT(LINE)).lt.1.d-4.and.IM.eq.2)
c     *       write(*,*) ANSTEE(LINE),GVW,IELEM(LINE)
          GAMTOT=GR(LINE)+GQST+GVW
          A=GAMTOT/(DNDOPL*4*PI)
          V=(WAVE-WLC)/DLDOPL
          VGT=VOIGT(A,V)/SQRTPI
C
C  Line absorption with Voigt function and total line opacity
C
          ALINE1=EFRACT*YABUND(LINE)*XSTIM*XNATOM/(XXRHO*DNDOPL)*VGT
c          write(*,*) TEMPER,ALINE1/FR/XNATOM/ABUND(IELEM(LINE)),
c          if(LINE.eq.2) then
c            write(*,*) TEMPER,OPCON,
c     *         FR*XNATOM*ABUND(IELEM(LINE)),
c     *         ALINE1/VGT/SQRTPI
c            write(*,'(F8.1,6e12.5)') TEMPER,ALINE1,VGT*SQRTPI,A,V,
c     *          DLDOPL,ABS(WAVE-WLC)
c          endif
c        END IF
        ALINE=ALINE+ALINE1
C
C  Store maximum contribution from given line to the total line opacity
C
        IF(ALINE1/OPCON.GT.ALMAX(LINE)) ALMAX(LINE)=ALINE1/OPCON
   3    CONTINUE
C
C  Compute the ratio according to model atmosphere type
C
        RATIO=RATIO/(ALINE+OPCON)
      END IF
C
      RETURN
      END

      FUNCTION CVOIGT(A,V)
C
C  Voigt function calculation: Humlicek's approximation
C  (Humlicek, J. 1982, J.Q.S.R.T. 27, 437)
C
C  COMPLEX approximation
C
      REAL*8 SAV
      COMPLEX*8 TAV,UAV,W4,V4
C
      TAV=CMPLX(A,-V)
      SAV=ABS(V)+A
      UAV=TAV*TAV
      IF(SAV.GE.15.d0) THEN
        W4=TAV*0.5641896d0/(0.5d0+UAV)
      ELSE IF(SAV.GE.5.5d0) THEN
        W4=TAV*(1.410474d0+UAV*0.5641896d0)/(0.75d0+UAV*(3.d0+UAV))
      ELSE IF(A.GE.0.195d0*V-0.176d0) THEN
        W4=(16.4955d0+TAV*(20.20933d0+TAV*(11.96482d0+
     +     TAV*(3.778987d0+TAV*0.5642236d0))))/(16.4955d0+
     +     TAV*(38.82363d0+TAV*(39.27121d0+TAV*(21.69274d0+
     +     TAV*(6.699398d0+TAV)))))
      ELSE
        W4=TAV*(36183.31d0-UAV*(3321.9905d0-UAV*(1540.787d0-
     -     UAV*(219.0313d0-UAV*(35.76683d0 -UAV*(1.320522d0-
     -     UAV*0.56419d0))))))
        V4=(32066.6d0-UAV*(24322.84d0-UAV*(9022.228d0-
     -      UAV*(2186.181d0-UAV*(364.2191d0-UAV*(61.57037d0-
     -      UAV*(1.841439d0-UAV)))))))
        W4=EXP(UAV)-W4/V4
      END IF
      CVOIGT=REAL(W4)
C
      RETURN
      END

      FUNCTION VOIGT(A,V)
C
C REAL approximation (slightly faster)
C
      TR= A
      TI=-V
      UR=A*A-V*V
      UI=-2*A*V
      SAV=ABS(V)+A
      IF(SAV.GE.15.) THEN
        UR=UR+0.5
        XX=MAX(A*A,V*V)
        TR=TR/XX
        TI=TI/XX
        UR=UR/XX
        UI=UI/XX
        VOIGT=0.5641896*(TR*UR+TI*UI)/(UR*UR+UI*UI)
      ELSE IF(SAV.GE.5.5) THEN
        X1=UR*0.5641896+1.410474
        Y1=UI*0.5641896
        XX=X1*TR-Y1*TI
        YY=X1*TI+Y1*TR
        X1=UR+3.
        Y1=UI
        UU=X1*UR-Y1*UI+0.75
        VV=X1*UI+Y1*UR
        VOIGT=(XX*UU+YY*VV)/(UU*UU+VV*VV)
      ELSE IF(A.GE.0.195*V-0.176) THEN
        X1=3.778987+TR*0.5642236
        Y1=         TI*0.5642236
        X2=X1*TR-Y1*TI+11.96482
        Y2=X1*TI+Y1*TR
        X1=X2*TR-Y2*TI+20.20933
        Y1=X2*TI+Y2*TR
        XX=X1*TR-Y1*TI+16.4955
        YY=X1*TI+Y1*TR
        X1=TR+6.699398
        Y1=TI
        X2=X1*TR-Y1*TI+21.69274
        Y2=X1*TI+Y1*TR
        X1=X2*TR-Y2*TI+39.27121
        Y1=X2*TI+Y2*TR
        X2=X1*TR-Y1*TI+38.82363
        Y2=X1*TI+Y1*TR
        UU=X2*TR-Y2*TI+16.4955
        VV=X2*TI+Y2*TR
        VOIGT=(XX*UU+YY*VV)/(UU*UU+VV*VV)
      ELSE
        X1=1.320522 -UR*0.56419
        Y1=         -UI*0.56419
        X2=35.76683 -(X1*UR-Y1*UI)
        Y2=         -(X1*UI+Y1*UR)
        X1=219.0313 -(X2*UR-Y2*UI)
        Y1=         -(X2*UI+Y2*UR)
        X2=1540.787 -(X1*UR-Y1*UI)
        Y2=         -(X1*UI+Y1*UR)
        X1=3321.9905-(X2*UR-Y2*UI)
        Y1=         -(X2*UI+Y2*UR)
        X2=36183.31 -(X1*UR-Y1*UI)
        Y2=         -(X1*UI+Y1*UR)
        XX=X2*TR-Y2*TI
        YY=X2*TI+Y2*TR
        X1=1.841439-UR
        Y1=        -UI
        X2=61.57037-(X1*UR-Y1*UI)
        Y2=        -(X1*UI+Y1*UR)
        X1=364.2191-(X2*UR-Y2*UI)
        Y1=        -(X2*UI+Y2*UR)
        X2=2186.181-(X1*UR-Y1*UI)
        Y2=        -(X1*UI+Y1*UR)
        X1=9022.228-(X2*UR-Y2*UI)
        Y1=        -(X2*UI+Y2*UR)
        X2=24322.84-(X1*UR-Y1*UI)
        Y2=        -(X1*UI+Y1*UR)
        UU=32066.6 -(X2*UR-Y2*UI)
        VV=        -(X2*UI+Y2*UR)
        VOIGT=EXP(UR)*COS(UI)-(XX*UU+YY*VV)/(UU*UU+VV*VV)
      END IF
C
      RETURN
      END

      FUNCTION VPF(A,V)
      COMPLEX Z
      SAVE
      DATA C0,C1,C2,C3,C4,C5,C6 /0.5641641,0.8718681,1.474395,
     1 -19.57862,802.4513,-4850.316,8031.468/
      DATA A0,A1,A2,A3,A4,A5,A6 /122.6079,214.3824,181.9285,93.15558,
     1 30.18014,5.912626,0.5641895/
      DATA B0,B1,B2,B3,B4,B5,B6 /122.6079,352.7306,457.3345,348.7039,
     1 170.3540,53.99291,10.47986/
      IF (A .LE. 1.0E-3 .AND. V .GE. 2.5) THEN
         V2 = V**2
         IF (V2 .GE. 100.0) THEN
            FACTOR = 0.0
         ELSE
            FACTOR = EXP(-V2)
         ENDIF
         VPF = FACTOR*(1.0 + (A**2)*(1.0 - 2.0*V2)) + (A/V2)*
     1    ((((((C6/V2 + C5)/V2 + C4)/V2 + C3)/V2 + C2)/V2 +C1)/V2 + C0)
      ELSE
         Z = CMPLX(A,-V)
         VPF = REAL(((((((A6*Z + A5)*Z + A4)*Z + A3)*Z + A2)*Z + A1)
     1    *Z + A0)/(((((((Z + B6)*Z + B5)*Z + B4)*Z + B3)*Z + B2)
     2    *Z + B1)*Z + B0))
      END IF
      RETURN
      END

      SUBROUTINE LIMBD(EPS,U,FC00,FC45,FC75)
C
C  THIS SUBROUTINE CALLS SUBROUTINE FCMU TO INTEGRATE THE EMMERGING
C  SPECIFIC INTENSITIES FOR CONTINUUM AT THE CENTER OF SPECTRAL
C  INTERVAL (RETURNED AS "FC*") AND FINDS LIMB DARKENING COEF.
C  U: FC=FC00*(1-U+U*mu) USING LINEAR REGRESSION FOR THREE ANGLES.
C
C  Author: N.Piskunov
C
C  LAST UPDATE: March 15, 1991.
C
      INCLUDE 'SIZES.SYN'
      INCLUDE 'COMMONS.SYN'
C
      DOUBLE PRECISION WLCONT,EPS
      PARAMETER (PI=3.14159265)
C
      XMU00=1.00
      WLCONT=WFIRST
      CALL FCMU(FC00,WLCONT,XMU00,EPS)
      write(*,*) WLCONT,FC00
      WLCONT=(WFIRST+WLAST)*0.5
      CALL FCMU(FC00,WLCONT,XMU00,EPS)
      write(*,*) WLCONT,FC00
      WLCONT=WLAST
      CALL FCMU(FC00,WLCONT,XMU00,EPS)
      write(*,*) WLCONT,FC00
c
c      WLCONT=(WFIRST+WLAST)*0.5
      WLCONT=WFIRST
C
C  CALCULATE CONTINUUM SPECIFIC INTENSITY AT THETA=0 (DISK CENTER)
C
      XMU00=1.00
      CALL FCMU(FC00,WLCONT,XMU00,EPS)
C
C  CALCULATE CONTINUUM SPECIFIC INTENSITY AT THETA=30
C
      XMU45=COS(PI*45./180.)
      CALL FCMU(FC45,WLCONT,XMU45,EPS)
C
C  CALCULATE CONTINUUM SPECIFIC INTENSITY AT THETA=60
C
      XMU75=COS(PI*75./180.)
      XMU75=0.5
      CALL FCMU(FC75,WLCONT,XMU75,EPS)
C
C  NOW WE SHELL LINEARLY APPROXIMATE LIMB DARKENING U: I=I0*(1-U+U*MU)
C
      A=1.-XMU45
      B=1.-XMU75
      U=(A*(1.-FC45/FC00)+B*(1.-FC75/FC00))/(A*A+B*B)
C      write(*,*) 'MU00=',XMU00,',  MU45=',XMU45,',  MU75=',XMU75
      write(*,*) 'FC00=', FC00,',  FC45=', FC45,',  FC75=', FC75
C
      RETURN
      END


      SUBROUTINE FCMU(FLUX,WLCONT,XMU,EPS)
C
C  THIS SUBROUTINE INTEGRATE THE SPECIFIC INTENSITY AT WAVELENGTH "WLCONT".
C  THE RESULT IS RETURNED AS "FLUX".
C
C  Author: N.Piskunov
C
C  LAST UPDATE: March 15, 1991.
C
      INCLUDE 'SIZES.SYN'
      INCLUDE 'COMMONS.SYN'
      PARAMETER (NSTOT=9)
      LOGICAL FIRSTS
      DOUBLE PRECISION TAULIN(11),WGHTS(10),BPLANK
      DOUBLE PRECISION WLCONT,FC,CONWL5,RATIO1,TAUNEW,XBEG,XEND,EPS
      EXTERNAL RATIO1
      SAVE TAULIN, WGHTS
      DATA TAULIN/
     1  0.000000000000, 0.137793470540, 0.729454549503, 1.808342901740,
     2  3.401433697855, 5.552496140064, 8.330152746764,11.843785837900,
     3 16.279257831378,21.996585811981,29.920697012274/
      DATA WGHTS /
     1  3.08441115765E-01, 4.01119929155E-01, 2.18068287612E-01,
     2  6.20874560987E-02, 9.50151697518E-03, 7.53008388588E-04,
     3  2.82592334960E-05, 4.24931398496E-07, 1.83956482398E-09,
     4  9.91182721961E-13/
C
C  INTEGRATE THE STANDARD WL OPTICAL DEPTH NODES FOR
C  CORRESPONDING TAULIN AT WL
C
      CONWL5=EXP(50.7649141-5.*LOG(WLCONT))
      HNUK=1.43868E+08/WLCONT
      TAUNEW=RHOX(1)
C
C  THE LOOP THROUGH THE ATMOSPHERE USING THE LOG DEPTH SCALE
C  STARTING AT THE SURFACE
C
      FC=0.D0
      FIRSTS=.TRUE.
      DO 1 IM=2,NSTOT
      XBEG=TAULIN(IM-1)*XMU
      XEND=TAULIN(IM)*XMU
      CALL RKINT(XBEG,XEND,EPS,TAUNEW,RHOX(2)-RHOX(1),WLCONT,1,
     *           FIRSTS,TEMPER,RATIO1)
      FIRSTS=.FALSE.
C
C  PLANK FUNCTION AGAIN
C
      BPLANK=CONWL5/(EXP(HNUK/TEMPER)-1.)
      FC=FC+BPLANK*WGHTS(IM-1)
   1  CONTINUE
      FLUX=FC
C
      RETURN
      END


      DOUBLE PRECISION FUNCTION RATIO1(TAU1,WL,ICODE,TEMPER)
C
C  THEIS ROUTINE DO ALMOST THE SAME AS "RATIO" EXCEPT THAT ONLY
C  CONTINUUM OPACITY AT THE CENTER OF INTERVAL IS CONSIDERED
C
C  Author: N.Piskunov
C
C  LAST UPDATE: March 15, 1991.
C
      INCLUDE 'SIZES.SYN'
      INCLUDE 'COMMONS.SYN'
      DOUBLE PRECISION WL,TAU,TAU1
C
      TAU=MAX(TAU1,RHOX(1))
      IM=2
   1  IF(RHOX(IM).LT.TAU.AND.IM.LT.NRHOX) THEN
        IM=IM+1
        GO TO 1
      END IF
      TEMPER=SPLINT(IM-1,IM,RHOX,T,     SP_T,  NRHOX,TAU)
      OPCONB=SPLINT(IM-1,IM,RHOX,COPBLU,SP_CBL,NRHOX,TAU)
      OPCONR=SPLINT(IM-1,IM,RHOX,COPRED,SP_CRD,NRHOX,TAU)
      IF(MOTYPE.EQ.0) THEN
        RATIO1=1.D0
      ELSE
        RATIO1=SPLINT(IM-1,IM,RHOX,COPSTD,SP_CST,NRHOX,TAU)
      END IF
C
      RATIO1=RATIO1*2.D0/(OPCONB+OPCONR)
C
      RETURN
      END


      SUBROUTINE RKINT(X,XEND,TOL,Y,DELTAY,WL,ICODE,FIRST,TEMPER,FUNC)
C
C  Integrate the standard optical depth Y for given values of
c  monochromatic optical depth X using 6th order
C  Runge-Kutta procedure. Integration is carried from X to XEND.
C  On input Y contains standard optical depth value corrsponding
C  to X. On return Y contains value corresponding to XEND.
C  Other input parameters:
C  TOL is the accuracy of integration, DELTAY is an estimate of
C  optical depth step, WL is the current wavelength, ICODE specifies
C  the independent depth scale of model atmosphere (mass or stadard
C  optical depth), FIRST is a flag to avoid recalculating certain
C  static variable.
C  Apart of Y, RKINT also return the temperature TEMPER in XEND.
C  FUNC is an external function which computes the ratio of opacities -
C  right part of differential equation.
C
C  Author: N.Piskunov
C
C  LAST UPDATE:  July 31, 1996.
C
      EXTERNAL FUNC
      DOUBLE PRECISION WL,WL1,Y,EPS,TINY,WK(10),FUNC,XTRIAL,
     *       YMAX,TOL,HMAX,HMIN,HMAG,HTRIAL,TEMP,EST,X,XEND,DELTAY
      DOUBLE PRECISION RK01,RK02,RK03,RK04,RK05,RK06,RK07,RK08,RK09,
     1                 RK10,RK11,RK12,RK13,RK14,RK15,RK16,RK17,RK18,
     2                 RK19,RK20,RK21,RK22,RK23,RK24,RK25,RK26,RK27,
     3                 RK28,RK29,RK30,RK31,RK32,RK33,RK34,RK35,RK36,
     4                 RK37,RK38,RK39,RK40,RK41,RK42,RK43
      PARAMETER (TINY=1.D-34,EPS =7.10706D-15)
      PARAMETER (            RK01=1.D0/6.D0,       RK02=4.D0/75.D0,
     1 RK03=16.D0/75.D0,     RK04=5.D0/6.D0,       RK05=-8.D0/3.D0,
     2 RK06=5.D0/2.D0,       RK07=-165.D0/64.D0,   RK08=55.D0/6.D0,
     3 RK09=-425.D0/64.D0,   RK10=85.D0/96.D0,     RK11=12.D0/5.D0,
     4 RK12=-8.D0,           RK13=4015.D0/612.D0,  RK14=-11.D0/36.D0,
     5 RK15=88.D0/255.D0,    RK16=-8263.D0/15000.D0,
     6 RK17=124.D0/75.D0,    RK18=-4501.D0/4760.D0,RK19=-81.D0/250.D0,
     7 RK20=2484.D0/10625.D0,RK21=3501.D0/1720.D0, RK22=-300.D0/43.D0,
     8 RK23=297275.D0/52632.D0,
     9 RK24=-319.D0/2322.D0, RK25=24068.D0/84065.D0,RK26=0.D0,
     A RK27=3850.D0/26703.D0,RK28=3.D0/40.D0,      RK29=0.D0,
     B RK30=875.D0/2244.D0,  RK31=23.D0/72.D0,
     C RK32=264.D0/1955.D0,  RK33=0.D0,
     D RK34=125.D0/11592.D0, RK35=43.D0/616.D0,    RK36=1.D0/160.D0,
     E RK37=0.D0,            RK38=125.D0/17952.D0, RK39=-1.D0/144.D0,
     F RK40=84.D0/13685.D0,  RK41=3.D0/44.D0,
     G RK42=-125.D0/11592.D0,RK43=-43.D0/616.D0)
      LOGICAL FIRST
      SAVE WL1,YMAX,WK,HINIT,HMAX,HMIN,HMAG,MXSTEP,EST,NSUCST,
     *     NSUCFL,NSTEP
C
C  Initial values
C
      IF(FIRST) THEN
        WL1=WL
        YMAX=MAX(1.D0,ABS(Y))
C
C  Calculate slope
C
        WK(1) =FUNC(Y,WL,ICODE,TEMPER)
        HINIT =DELTAY/(5.D0*WK(1))
        IF(HINIT.EQ.0.D0) THEN
          HMAG=0.5D0*HMAX*TOL**(1.D0/6.D0)
        ELSE
          HMAG=0.5D0*HINIT
        END IF
        HMAX  =2.D0
        MXSTEP=50000
C
C  Initialize EPS and HMAG
C
        WK(2)=0.D0
        EST=0.D0
        NSUCST=0
        NSUCFL=0
        NSTEP =0
      END IF
C
C  Loop through the following four stages, once for each trial step
C  until the occurrence of one of the following (A) Normal return
C  on reaching XEND in stage 4, or (B) An error stop in stage 1 or 4.
C
C  Stage 1 -- Do calculations of HMIN, HMAX, etc., and some
C  parameter checking. End up with suitable values of HMAG, XTRIAL
C  and HTRIAL for use in the integration step. Check MXFCN
C
   1  CONTINUE
C
C  Check MXSTEP
C
      IF(NSTEP.GE.MXSTEP) THEN
        WRITE(*,'('' Maximum number of steps allowed ('',I6,
     *            '') was used''/
     *            '' The problem may be stiff'')') MXSTEP
        WRITE(*,*) 'WL=',WL
        STOP
      END IF
C
C  Calculate HMIN
C
      HMIN=10.D0*MAX(TINY,EPS*MAX(ABS(XEND),ABS(X)))
C
C  Error stop if HMIN .GT. HMAX
C
      IF(HMIN.GT.HMAX) THEN
        WRITE(*,'('' HMIN ='',G12.6,'' is greater than HMAX ='',
     *             G12.6)') HMIN,HMAX
        STOP
      END IF
C
C  Calculate preliminary HMAG:
C
      IF(NSUCFL.LE.1) THEN
C
C  After a successful step, or at most one failure (IF is to avoid overflow
C  and MAX to avoid reduction by more than half)
C
        TEMP=2.D0*HMAG
        IF(TOL.LT.(2.D0/0.9D0)**6*EST)
     *                   TEMP=0.9D0*(TOL/EST)**(1.D0/6.D0)*HMAG
        HMAG=MAX(TEMP,0.5D0*HMAG)
      ELSE
C
C  After two or more successive failures
C
        HMAG=0.5D0*HMAG
      END IF
C
C  Check against HMAX and HMIN
C
      HMAG=MAX(MIN(HMAG,HMAX),HMIN)
C
C  Calculate HMAG, XTRIAL - depending on preliminary HMAG, XEND
C
      IF(HMAG.LT.ABS(XEND-X)) THEN
C
C  Do not step more than half way to XEND
C
        HMAG=MIN(HMAG,0.5D0*ABS(XEND-X))
        XTRIAL=X+SIGN(HMAG,XEND-X)
      ELSE
C
C  Hit XEND exactly
C
        HMAG=ABS(XEND-X)
        XTRIAL=XEND
      END IF
C
C  Calculate HTRIAL
C
      HTRIAL=XTRIAL-X
C
C  Stage 2 -- Calculate YTRIAL. WK(2),..., WK(8) will hold intermediate
C  results needed in stage 3. WK(9) is temporary storage until finally
C  it holds YTRIAL
C
      WK(9)=Y+HTRIAL*WK(1)*RK01
      WK(2)=FUNC(WK(9),WL,ICODE,TEMPER)
C      IF(WK(9).LT.Y) write(*,*) '2:',Y,HTRIAL,WK(9)
C
      WK(9)=Y+HTRIAL*(WK(1)*RK02+WK(2)*RK03)
      WK(3)=FUNC(WK(9),WL,ICODE,TEMPER)
C      IF(WK(9).LT.Y) write(*,*) '3:',Y,HTRIAL,WK(9)
C
      WK(9)=Y+HTRIAL*(WK(1)*RK04+WK(2)*RK05+WK(3)*RK06)
      WK(4)=FUNC(WK(9),WL,ICODE,TEMPER)
C      IF(WK(9).LT.Y) write(*,*) '4:',Y,HTRIAL,WK(9)
C
      WK(9)=Y+HTRIAL*(WK(1)*RK07+WK(2)*RK08+WK(3)*RK09+WK(4)*RK10)
      WK(5)=FUNC(WK(9),WL,ICODE,TEMPER)
C      IF(WK(9).LT.Y) write(*,*) '5:',Y,HTRIAL,WK(9)
C
      WK(9)=Y+HTRIAL*(WK(1)*RK11+WK(2)*RK12+WK(3)*RK13+WK(4)*RK14+
     +                WK(5)*RK15)
      WK(6)=FUNC(WK(9),WL,ICODE,TEMPER)
C      IF(WK(9).LT.Y) write(*,*) '6:',Y,HTRIAL,WK(9)
C
      WK(9)=Y+HTRIAL*(WK(1)*RK16+WK(2)*RK17+WK(3)*RK18+WK(4)*RK19+
     +                WK(5)*RK20)
      WK(7)=FUNC(WK(9),WL,ICODE,TEMPER)
C      IF(WK(9).LT.Y) write(*,*) '7:',Y,HTRIAL,WK(9)
C
      WK(9)=Y+HTRIAL*(WK(1)*RK21+WK(2)*RK22+WK(3)*RK23+WK(4)*RK24+
     +                WK(5)*RK25+WK(6)*RK26+WK(7)*RK27)
      WK(8)=FUNC(WK(9),WL,ICODE,TEMPER)
C      IF(WK(9).LT.Y) write(*,*) '8:',Y,HTRIAL,WK(9)
C
C  Calculate YTRIAL, the extrapolated approximation and store in WK(9)
C
      WK(9)=Y+HTRIAL*(WK(1)*RK28+WK(2)*RK29+WK(3)*RK30+WK(4)*RK31+
     +                WK(5)*RK32+WK(6)*RK33+WK(7)*RK34+WK(8)*RK35)
C
C  Stage 3 -- Calculate the error estimate EST Calculate the
C             unweighted absolute error estimate vector
C
      WK(2)=WK(1)*RK36+WK(2)*RK37+WK(3)*RK38+WK(4)*RK39+
     *      WK(5)*RK40+WK(6)*RK41+WK(7)*RK42+WK(8)*RK43
C
C  Calculate EST - the weighted max norm of WK(2)
C  EST is intended to be a measure of the error per unit step in YTRIAL
C
      TEMP=ABS(WK(2))/YMAX
      EST=TEMP*HMAG
C
C  Stage 4 -- Make decision
C
      IF(EST.LE.TOL) THEN
C
C  Step accepted, so update X, Y from XTRIAL and WK(9), also YMAX
C
        X=XTRIAL
        Y=WK(9)
        YMAX=MAX(YMAX,ABS(Y))
        NSUCST=NSUCST+1
        NSUCFL=0
        NSTEP =NSTEP+1
        WK(1)=FUNC(Y,WL,ICODE,TEMPER)
C
C  Normal return when XEND has been reached
C
        IF(X.EQ.XEND) RETURN
      ELSE
C
C  Step not accepted - Add 1 to number successive failures
C
        NSUCFL=NSUCFL+1
C
C  Error stop if HMAG .LE. HMIN
C
        IF(HMAG.LE.HMIN) THEN
          WRITE(*,'('' Unable to satisfy the error requirement''/
     *              '' TOL='',G12.6,'' may be too small'')') TOL
          STOP
        END IF
      END IF
C
C  End stage 4
C
      GO TO 1
      END

      SUBROUTINE GAMHE(IND,T,ANE,ANP,GAM,SHIFT)
C
C   NEUTRAL HELIUM STARK BROADENING PARAMETERS
C   AFTER DIMITRIJEVIC AND SAHAL-BRECHOT, 1984, J.Q.S.R.T. 31, 301
C   OR FREUDENSTEIN AND COOPER, 1978, AP.J. 224, 1079  (FOR C(IND).GT.0)
C
      DOUBLE PRECISION W(5,20),V(4,20),
     *                 SHIFTE(4,20),SHIFTP(4,20),C(20)
      SAVE W,V,SHIFTE,SHIFTP,C
C
C   ELECTRONS T= 5000   10000   20000   40000     LAMBDA
C
      DATA W /  5.990,  6.650,  6.610,  6.210,    3819.60,
     1          2.950,  3.130,  3.230,  3.300,    3867.50,
     2        109.000, 94.400, 79.500, 65.700,    3871.79,
     3          0.142,  0.166,  0.182,  0.190,    3888.65,
     4         70.700, 60.700, 50.900, 41.900,    3926.53,
     5          1.540,  1.480,  1.400,  1.290,    3964.73,
     6         41.600, 50.500, 57.400, 65.800,    4009.27,
     7          1.320,  1.350,  1.380,  1.460,    4120.80,
     8          7.830,  8.750,  8.690,  8.040,    4143.76,
     9          5.830,  6.370,  6.820,  6.990,    4168.97,
     A          2.280,  2.320,  2.360,  2.430,    4437.55,
     B          2.470,  2.200,  1.910,  1.650,    4471.50,
     C          0.588,  0.620,  0.641,  0.659,    4713.20,
     D          2.600,  2.480,  2.240,  1.960,    4921.93,
     E          0.627,  0.597,  0.568,  0.532,    5015.68,
     F          1.050,  1.090,  1.110,  1.140,    5047.74,
     G          0.277,  0.298,  0.296,  0.293,    5875.70,
     H          0.714,  0.666,  0.602,  0.538,    6678.15,
     I          3.490,  3.630,  3.470,  3.190,    4026.20,
     J          4.970,  5.100,  4.810,  4.310,    4387.93/
C
C   PROTONS   T= 5000   10000   20000   40000
C
      DATA V /  1.520,  4.540,  9.140, 10.200,
     1          0.607,  0.710,  0.802,  0.901,
     2          0.000,  0.000,  0.000,  0.000,
     3          0.0396, 0.0434, 0.0476, 0.0526,
     4          0.000,  0.000,  0.000,  0.000,
     5          0.507,  0.585,  0.665,  0.762,
     6          0.930,  1.710, 13.600, 27.200,
     7          0.288,  0.325,  0.365,  0.410,
     8          1.330,  6.800, 12.900, 14.300,
     9          1.100,  1.370,  1.560,  1.760,
     A          0.516,  0.579,  0.650,  0.730,
     B          1.520,  1.730,  1.830,  1.630,
     C          0.128,  0.143,  0.161,  0.181,
     D          2.040,  2.740,  2.950,  2.740,
     E          0.187,  0.210,  0.237,  0.270,
     F          0.231,  0.260,  0.291,  0.327,
     G          0.0591, 0.0650, 0.0719, 0.0799,
     H          0.231,  0.260,  0.295,  0.339,
     I          2.180,  3.760,  4.790,  4.560,
     J          1.860,  5.320,  7.070,  7.150/
C
C  Shifts due to electrons
C
      DATA SHIFTE/-0.698, -0.558, -0.354, -0.216,
     1          1.800,  1.930,  1.810,  1.670,
     2          8.510,  5.340,  2.560,  1.560,
     3          0.075,  0.061,  0.049,  0.035,
     4          7.130,  4.270,  1.960,  0.560,
     5         -0.459, -0.345, -0.249, -0.179,
     6         10.400, 20.700, 29.700, 38.000,
     7          0.890,  0.931,  0.851,  0.677,
     8          0.924,  0.856,  0.775,  0.656,
     9          3.120,  3.430,  3.490,  3.500,
     A          1.690,  1.600,  1.270,  0.906,
     B          0.062, -0.064, -0.015, -0.006,
     C          0.409,  0.456,  0.439,  0.349,
     D          0.436,  0.368,  0.298,  0.221,
     E         -0.236, -0.179, -0.132, -0.095,
     F          0.730,  0.745,  0.668,  0.528,
     G         -0.073, -0.040, -0.012, -0.005,
     H          0.249,  0.222,  0.180,  0.144,
     I         -0.425, -0.315, -0.209, -0.136,
     J          0.665,  0.558,  0.450,  0.336/
C
C  Shifts due to protons
C
      DATA SHIFTP/ 0.000,  0.055,  1.790,  6.100,
     1          0.243,  0.422,  0.579,  0.725,
     2          0.000,  0.000,  0.000,  0.000,
     3          0.028,  0.033,  0.039,  0.044,
     4          0.000,  0.000,  0.000,  0.000,
     5         -0.232, -0.367, -0.488, -0.602,
     6          0.000,  0.000,  0.089,  4.630,
     7          0.170,  0.234,  0.294,  0.351,
     8          0.000,  0.028,  1.540,  6.750,
     9          0.280,  0.676,  1.030,  1.340,
     A          0.465,  0.532,  0.604,  0.684,
     B          1.350,  1.560,  1.840,  2.110,
     C          0.094,  0.117,  0.139,  0.161,
     D          0.261,  1.140,  2.010,  2.650,
     E         -0.131, -0.164, -0.197, -0.231,
     F          0.158,  0.203,  0.246,  0.288,
     G         -0.045, -0.052, -0.060, -0.069,
     H          0.171,  0.211,  0.250,  0.292,
     I          0.002,  0.544,  2.200,  3.680,
     J          0.001,  0.359,  2.770,  5.140/
C
      DATA C /2*0.,1.83E-4,0.,1.13E-4,5*0.,1.6E-4,9*0./
      PARAMETER (TT1=3.699, TT2=4., TT3=4.301, TT4=4.602)
C
      IF(W(4,IND).EQ.0.) THEN
        WRITE(*,*) ' NO DATA AVAILABLE FOR He I',W(5,IND)
        STOP
      END IF
      IF(W(1,IND).NE.0.0) THEN
C
C CUBIC INTERPOLATION OVER T=5000,10000,20000,40000 IN LOG SCALE
C
        TLG=LOG10(T)
        IF(TLG.LE.TT3) THEN
          J=3
          TJ =(TT3-TT2)*(TT3-TT1)*(TT2-TT1)
          TJ0=(TLG-TT1)*(TLG-TT2)*(TT2-TT1)/TJ
          TJ1=(TLG-TT1)*(TT3-TLG)*(TT3-TT1)/TJ
          TJ2=(TLG-TT2)*(TLG-TT3)*(TT3-TT2)/TJ
        ELSE
          J=4
          TJ =(TT4-TT3)*(TT4-TT2)*(TT3-TT2)
          TJ0=(TLG-TT2)*(TLG-TT3)*(TT3-TT2)/TJ
          TJ1=(TLG-TT2)*(TT4-TLG)*(TT4-TT2)/TJ
          TJ2=(TLG-TT3)*(TLG-TT4)*(TT4-TT3)/TJ
        END IF
        GAM=((TJ0*W(J,IND)+TJ1*W(J-1,IND)+TJ2*W(J-2,IND))*ANE
     *      +(TJ0*V(J,IND)+TJ1*V(J-1,IND)+TJ2*V(J-2,IND))*ANP)
     *      *1.884E3/W(5,IND)**2
        IF(GAM.LT.0.) GAM=0.
C
        SHIFT=(TJ0*SHIFTE(J,IND)+TJ1*SHIFTE(J-1,IND)+
     *         TJ2*SHIFTE(J-2,IND))*(ANE/1.E16) +
     *        (TJ0*SHIFTP(J,IND)+TJ1*SHIFTP(J-1,IND)+
     *         TJ2*SHIFTP(J-2,IND))*(ANP/1.E16)
C
      ELSE
        GAM=C(IND)*T**0.16667*ANE
        SHIFT=0
      END IF
C
      RETURN
      END
