C++******************************************************************
C OBJECT
C To create an image for simulations
C Written initially by A. Lannes
C
C Syntax:    OBJET NX,NY FILE
C Example:   objet 64,64 test.bdf
C
C JLP
C Version:  18-05-92
C--******************************************************************
	PROGRAM OBJECT
	PARAMETER(NMAX=40,NCLMAX=256)
	CHARACTER NAME*40, COMMENTS*80
	INTEGER X(NMAX),Y(NMAX),Z(NMAX)
	REAL A(NMAX),B(NMAX),SUM
	REAL OBJET(NCLMAX,NCLMAX)
	CALL JLP_BEGIN
	CALL JLP_INQUIFMT
 
	WRITE(6,11)
11	FORMAT(' Program OBJECT, Version 18-05-92',/,
     1	' To create a test object for simulations')
	WRITE(6,12)
12	FORMAT(' Size of the output image (less than 256,256):')
	READ(5,*) NC,NL
 
	PI=3.14159265358979323846
 
C DEFINITION OBJET
C POSITION : X,Y ; Z : PARAM GAUSS ; A : AMPLITUDE
	NO = 26
 
C POINT 1 FOND CONTINU
        I=1
	X(I)=0
	Y(I)=0
	Z(I)=30
	A(I)=0.5
 
C POINT 2
        I=2
	X(I)=-8.
	Y(I)=-8.
	Z(I)=4
	A(I)=1
 
C Can be zeroed if needed:
C POINT 3
        I=3
	X(I)=8.
	Y(I)=8
	Z(I)=4
	A(I)=0.1
 
C Big spot on the middle:
C POINT 4
        I=4
	X(I)=-7
	Y(I)=8
	Z(I)=10
	A(I)=2.2
 
C Can be zeroed if needed:
C POINT 5
        I=5
	X(I)=-10
	Y(I)=-16
	Z(I)=4
	A(I)=0.1
 
C Two points on right:
C POINT 6
        I=6
	X(I)=26
	Y(I)=22
	Z(I)=4
	A(I)=2.
 
C POINT 13
        I=13
	X(I)=12
	Y(I)=24
	Z(I)=4
	A(I)=1.9
 
C POINT 7
        I=7
	X(I)= -4
	Y(I)=-7
	Z(I)=4
	A(I)=1.6
 
C POINT 8
        I=8
	X(I)= -24
	Y(I)=4
	Z(I)=4
	A(I)=1.2
 
C POINT 9
        I=9
	X(I)=3
	Y(I)=14
	Z(I)=4
	A(I)=0.5
 
C 3 points on top:
C POINT 10
        I=10
	X(I)=-13
	Y(I)=32
	Z(I)=4.5
	A(I)=2.2
 
C POINT 21
        I=21
	X(I)=-4
	Y(I)=25
	Z(I)=4.5
	A(I)=2.2
 
C POINT 22
        I=22
	X(I)=2
	Y(I)=34
	Z(I)=4.5
	A(I)=2.2

C POINT 12
        I=12
	X(I)=12
	Y(I)=-28
	Z(I)=4
	A(I)=1.8
 
C POINT 14
        I=14
	X(I)=-34
	Y(I)=20
	Z(I)=6
	A(I)=1.7
 
C Can be zeroed if needed: 
C POINT 15
        I=15
	X(I)=-22
	Y(I)=5
	Z(I)=5
	A(I)=.9
 
C POINT 16
        I=16
	X(I)=-24
	Y(I)=8
	Z(I)=4
	A(I)=.5
 
C POINT 17
        I=17
	X(I)=-34
	Y(I)=20
	Z(I)=5
	A(I)=1.
 
C Big spot on the bottom:
C POINT 18
        I=18
	X(I)=-14
	Y(I)=-22
	Z(I)=10
	A(I)=2.3
 
C Can be zeroed if needed:
C POINT 19
        I=19
	X(I)=-12
	Y(I)=14
	Z(I)=4
	A(I)=0.1
 
C POINT 20
        I=20
	X(I)=4
	Y(I)=-30
	Z(I)=4
	A(I)=2.0
 
C Points to fill the gaps:
C POINT 23
        I=23
	X(I)=-20
	Y(I)=+20
	Z(I)=5
	A(I)=2.0
 
C POINT 24
        I=24
	X(I)=0
	Y(I)=-20
	Z(I)=4
	A(I)=2.0
 
C POINT 25
        I=25
	X(I)=-20
	Y(I)=-10
	Z(I)=4
	A(I)=2.0
 
 
C Two big spots on the right 
C POINT 11
        I=11
	X(I)=+17
	Y(I)=+8
	Z(I)=10
	A(I)=2.1
 
C POINT 26
        I=26
	X(I)=+14
	Y(I)=-13
	Z(I)=10
	A(I)=2.2
 
C RECENTRAGE ET CALCUL DES PARAMETRES
	IXC=(NC/2)+1
	IYC=(NL/2)+1
	DO 1 N=1,NO	
	  X(N) = IXC+X(N)*FLOAT(NC)/128.
	  Y(N) = IYC+Y(N)*FLOAT(NL)/128.
	  B(N) = -PI/FLOAT(Z(N)**2)
1	CONTINUE	
 
C INTENSITE OBJET (GLOBAL) CENTRE
C OBTENUE LIGNE PAR LIGNE, POINT PAR POINT
 
        SUM=0.
	DO 2 J=1,NL
	DO 3 I=1,NC
	XO=0.
C EN CHAQUE POINT DU RESEAU LA CONTRIBUTION DE CHAQUE OBJET
	   DO 4 N=1,NO
	   TT=FLOAT((I-X(N))**2+(J-Y(N))**2)
	   XO = XO + A(N)*EXP(B(N)*TT)
4	   CONTINUE
 
	OBJET(I,J)= XO
        SUM=SUM+XO
3	CONTINUE
2	CONTINUE
 
C Normalization to one:
	DO 22 J=1,NL
	 DO 23 I=1,NC
           OBJET(I,J)=OBJET(I,J)/SUM  
23       CONTINUE
22      CONTINUE

C ECRITURE DU TABLEAU OBJET DANS FICHIER OBJET.MID
	COMMENTS=' Simulated object'
	NAME=' '
	CALL JLP_WRITEIMAG(OBJET,NC,NL,NCLMAX,NAME,COMMENTS)
	CALL JLP_END
	STOP
	END
