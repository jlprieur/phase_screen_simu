       SUBROUTINE FOURN1(DATA,DL,DM,ISIGN)
C++*************************************************************
C De Vos version (FOUR2)
C ---  -----------------------------------------------------------
C      Four2 makes a fast fourier transform of the complex data
C      in data (data(i) is real part, data(i+1) is complex part.
C      The axes have length dl and dm. Set isign to either -1
C      or 1 (inv), in the last case, data is scaled by 1/(dl*dm).
C ---  -----------------------------------------------------------
C
       IMPLICIT NONE
        CHARACTER ID_X*(36)
        DATA      ID_X/'!@(#) four2.f    23/08/91 v3.01 >'/
       REAL    DATA(*)
       INTEGER DL,DM,ISIGN
       REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
       integer i1, i2, i2rev, i3, i3rev, ibit, idim, ifp1, ifp2
       integer ip1, ip2, ip3, itot, k1, k2, n, nprev, nrem, ntot
       real    fac, tempi, tempr
c
       if (isign.eq.0) return
c
       NTOT=DL*DM
       if (isign.gt.0) then
          fac = 1.0/float(ntot)
          DO 102 itot = 1, 2*ntot
             data(itot) = fac * data(itot)
  102     CONTINUE
       endif
c
       NPREV=1
       N=DL
       DO 18 IDIM=1,2
*        N=NN(IDIM)
         NREM=NTOT/(N*NPREV)
         IP1=2*NPREV
         IP2=IP1*N
         IP3=IP2*NREM
         I2REV=1
         DO 14 I2=1,IP2,IP1
           IF(I2.LT.I2REV)THEN
             DO 13 I1=I2,I2+IP1-2,2
               DO 12 I3=I1,IP3,IP2
                 I3REV=I2REV+I3-I2
                 TEMPR=DATA(I3)
                 TEMPI=DATA(I3+1)
                 DATA(I3)=DATA(I3REV)
                 DATA(I3+1)=DATA(I3REV+1)
                 DATA(I3REV)=TEMPR
                 DATA(I3REV+1)=TEMPI
12             CONTINUE
13           CONTINUE
           ENDIF
           IBIT=IP2/2
1          IF ((IBIT.GE.IP1).AND.(I2REV.GT.IBIT)) THEN
             I2REV=I2REV-IBIT
             IBIT=IBIT/2
           GO TO 1
           ENDIF
           I2REV=I2REV+IBIT
14       CONTINUE
         IFP1=IP1
2        IF(IFP1.LT.IP2)THEN
           IFP2=2*IFP1
           THETA=ISIGN*6.28318530717959D0/(IFP2/IP1)
           WPR=-2.D0*DSIN(0.5D0*THETA)**2
           WPI=DSIN(THETA)
           WR=1.D0
           WI=0.D0
           DO 17 I3=1,IFP1,IP1
             DO 16 I1=I3,I3+IP1-2,2
               DO 15 I2=I1,IP3,IFP2
                 K1=I2
                 K2=K1+IFP1
                 TEMPR=SNGL(WR)*DATA(K2)-SNGL(WI)*DATA(K2+1)
                 TEMPI=SNGL(WR)*DATA(K2+1)+SNGL(WI)*DATA(K2)
                 DATA(K2)=DATA(K1)-TEMPR
                 DATA(K2+1)=DATA(K1+1)-TEMPI
                 DATA(K1)=DATA(K1)+TEMPR
                 DATA(K1+1)=DATA(K1+1)+TEMPI
15             CONTINUE
16           CONTINUE
             WTEMP=WR
             WR=WR*WPR-WI*WPI+WR
             WI=WI*WPR+WTEMP*WPI+WI
17         CONTINUE
           IFP1=IFP2
         GO TO 2
         ENDIF
         NPREV=N*NPREV
         N=DM
18     CONTINUE
       RETURN
       END
