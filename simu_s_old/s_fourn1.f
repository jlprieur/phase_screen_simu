      SUBROUTINEFOURN1(DATA,DL,DM,ISIGN)
      IMPLICITNONE
      CHARACTERID_X*(36)
      DATAID_X/'!@(#) four2.f    23/08/91 v3.01 >'/
      REALDATA(*)
      INTEGERDL,DM,ISIGN
      REAL*8WR,WI,WPR,WPI,WTEMP,THETA
      integeri1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2
      integerip1,ip2,ip3,itot,k1,k2,n,nprev,nrem,ntot
      realfac,tempi,tempr
      if(isign.eq.0)return
      NTOT=DL*DM
      if(isign.gt.0)then
      fac=1.0/float(ntot)
      DO102itot=1,2*ntot
      data(itot)=fac*data(itot)
  102 CONTINUE
      endif
      NPREV=1
      N=DL
      DO18IDIM=1,2
      NREM=NTOT/(N*NPREV)
      IP1=2*NPREV
      IP2=IP1*N
      IP3=IP2*NREM
      I2REV=1
      DO14I2=1,IP2,IP1
      IF(I2.LT.I2REV)THEN
      DO13I1=I2,I2+IP1-2,2
      DO12I3=I1,IP3,IP2
      I3REV=I2REV+I3-I2
      TEMPR=DATA(I3)
      TEMPI=DATA(I3+1)
      DATA(I3)=DATA(I3REV)
      DATA(I3+1)=DATA(I3REV+1)
      DATA(I3REV)=TEMPR
      DATA(I3REV+1)=TEMPI
   12 CONTINUE
   13 CONTINUE
      ENDIF
      IBIT=IP2/2
    1 IF((IBIT.GE.IP1).AND.(I2REV.GT.IBIT))THEN
      I2REV=I2REV-IBIT
      IBIT=IBIT/2
      GOTO1
      ENDIF
      I2REV=I2REV+IBIT
   14 CONTINUE
      IFP1=IP1
    2 IF(IFP1.LT.IP2)THEN
      IFP2=2*IFP1
      THETA=ISIGN*6.28318530717959D0/(IFP2/IP1)
      WPR=-2.D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.D0
      WI=0.D0
      DO17I3=1,IFP1,IP1
      DO16I1=I3,I3+IP1-2,2
      DO15I2=I1,IP3,IFP2
      K1=I2
      K2=K1+IFP1
      TEMPR=SNGL(WR)*DATA(K2)-SNGL(WI)*DATA(K2+1)
      TEMPI=SNGL(WR)*DATA(K2+1)+SNGL(WI)*DATA(K2)
      DATA(K2)=DATA(K1)-TEMPR
      DATA(K2+1)=DATA(K1+1)-TEMPI
      DATA(K1)=DATA(K1)+TEMPR
      DATA(K1+1)=DATA(K1+1)+TEMPI
   15 CONTINUE
   16 CONTINUE
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
   17 CONTINUE
      IFP1=IFP2
      GOTO2
      ENDIF
      NPREV=N*NPREV
      N=DM
   18 CONTINUE
      RETURN
      END
