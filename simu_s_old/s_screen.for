        program s_pscorg
*******************************************************************************
*
*       is a phase screen simulator which produces artificial fringe
*       or speckle patterns and obtain integrated power spectrum or(/and)
*       optionally specified bispectral components.
*       The atmospheric disturbance is given by the Kolmogorov spectrum.
*       Seeing condition, aperture configuration and the object ( at this
*       point, only multiple point sources ) and the usage of the phase screen
*       must be specified by the user in an ascii input file.
*
*       Following is a sample input file. (free format, ! indicates a
*       comment line.)
*
*       ! seed for random number generator (large odd integer)
*       Random
*       7597551
*       ! outer scale and inner scale of the turbulence(m) - recommended values
*       Turbulence
*       10.0  0.01
*       ! central wavelength(m), fractional bandwidth, sampling points
*       Wavelength
*       0.65e-6 0.016 3
*       ! source: # of point sources
*       !         intensity(arbitrary unit), position offset(x,y) in arcsecond
*       Source
*       2
*       1.0   0.0  0.0
*       0.5   0.02 0.02
*       ! # of circular apertures (including voids)
*       ! configuraton: center of aperture(x,y), radius(m) and type
*       ! 1 is an aperture,0 is a void. Voids should follow ordinary apertures.
*       ! Example: annular mask
*       Aperture
*       2
*       0.0 0.0 2.0  1
*       0.0 0.0 1.9  0
*       ! Multiple subsets are taken out of one phase screen.
*       ! Equally spaced sampling is specified by
*       ! x_min, x_max, step size in x(m)
*       ! y_min, y_max, step size in y(m)
*       Phase_screen
*       1.25 8.75  0.14
*       1.25 8.75  2.5
*       ! total number of frames to be processed
*       Frame
*       100
*       ------ as an option `no atmosphere' can be specified by ------
*       ! obtaining beam pattern
*       No_disturbance
*       --------------------------------------------------------------
*
*       Written in standard FORTRAN and expected to run both on VMS and UNIX
*       simply by replacing FFT routines and random number generator.
*       On UNIX all tempory arrays are to be dynamically allocated by 'calloc()'
*       Source code should be called pscreen.for on VMS, pscreen.f on UNIX.
*       On a vector machine, effeciency of loops should be reconsidered.
*
*       This is a re-organized version of the pscsn.for developed on Tacos.
*       The program is planned to be as users oriented as possible for the
*       future development of Optical/IR interferometry.
*
*       version 1.0    Dec 25, 1988 (Merry Christmas!)     Tadashi Nakajima
*       Today's resolution : `I will have real holidays next year and live
*                             a human life.'
*       version 2.0    Dec 30, 1988   Now sharp cut-offs of fluctuation power
*       spectrum at kapmin and kapmax are introduced to investigatate the
*       relative significance of the diffractive effect and refractive effect.
*       For instance, the refractive effect can be seen by setting l_o = 0.1,
*       while the diffractive effect can be investigated by setting l_0 = 0.4.
*******************************************************************************
*------ declarations ---------*
*       constants
        real*4    pi, twopi
 
*       integrated structure constant of the phase fluctuation
        real*4 lcn2
 
*       size of phase screen and resolution (m/pix)
        real*4 l_s, l_p
        parameter( l_s = 10.24 )
        parameter( l_p = 0.02  )
 
*  size of phase screen
        integer*4 nsise
        parameter( nsize = 256  )
 
*       wave number, spatial frequencies and scale length of tubulence
        real kk, kapmax, kapmin, l_0, dkap, lambda, l_o
 
*       complex phase screen, real and imag are independent.
        real*4 pha_re(nsize,nsize),pha_im(nsize,nsize)
 
*       input control parameter file and output file(s)
        character*30 ifile, ofile
 
*       I/O units
        integer*4 f_in
        parameter( f_in = 11)
 
*       variables local to main program
        character*80 comment,ocomments
        integer*4 seed
 
*       common variables are not compatible with parameter statement.
        pi = 3.141592653
        twopi = 2.*pi
        lcn2 = 5.0e-13
 
* JLP
        call jlp_begin
        call jlp_inquifmt
 
*       input file names
        write(*,10)
   10   format(1x,' input ascii file name : ',$)
        read(*,'(a)') ifile
 
*       read the input file
        open(unit=f_in, file=ifile, status='old')
 
   30   format(a80)
   35   format(1x,a80)
        do 40 j=1,100
 
          read(f_in,30,end=100) comment
          write(*,35) comment
 
          if( comment(1:1) .eq. '!') goto 40
C Keyword "Random", parameter: seed (random generator): 
          if( comment(1:1) .eq. 'R' .or. comment(1:1) .eq. 'r') then
            read(f_in,*) seed
          end if
C Keyword "Turbulence", parameters: outer and inner  scales
          if( comment(1:1) .eq. 'T' .or. comment(1:1) .eq. 't') then
            read(f_in,*) l_0, l_o
          end if
C Keyword "Wavelength", parameters: central wavel., fractional bandwidth, 
C sampling points
          if( comment(1:1) .eq. 'W' .or. comment(1:1) .eq. 'w') then
            read(f_in,*) lambda
          end if
 
   40   continue
  100   continue

            write(6,36) seed,l_0,l_o,lambda 
            write(ocomments,36) seed,l_0,l_o,lambda 
36	  format(' seed=',I6,' L0 = ',G10.4,' lo = ',G10.4,
     &    ' Lambda = ',G10.4)
* kk: wave number, kap: spatial frequency in aperture plane
 
        kk = twopi / lambda
 
        work=min(l_s,l_0)
        kapmin = twopi / work
 
        work=max(l_o,l_p)
        kapmax = twopi / work 
 
* dkap: one pixel in spatial frequency in aperture plane
        dkap = twopi / l_s
*
* create a complex phase screen
        call jlp_random_init(seed)
        call screen(pha_re,pha_im,nsize,kk,kapmax,kapmin,lcn2,dkap)

        ofile=' '
        print *,' output of phase screen (real part)'
        call jlp_writeimag(pha_re,nsize,nsize,nsize,ofile,ocomments)
        ofile=' '
        print *,' output of phase screen (imaginary part)'
        call jlp_writeimag(pha_im,nsize,nsize,nsize,ofile,ocomments)
*
        call jlp_end
        stop
        end
 
 
***************************************************************************
* screen
* creates complex phase screen
***************************************************************************
        subroutine screen(pha_re,pha_im,nsize,kk,kapmax,kapmin,lcn2,dkap)
        real*4 pha_re(nsize,nsize),pha_im(nsize,nsize)
        real*4 fs, fso, pi, twopi
        real*4 kk, kapmax, kapmin, lcn2, dkap, u, wrand, work
 
*       radial spatial frequency in aperture plane
        real*4 kap
 
*       dimension array for fft
*       origin of the phase screen
        integer*4 Nc

        pi = 3.141592653
        twopi = 2.*pi
 
*------ generation of random phase in spatial domain -------------
*       Complex Gaussian random numbers are generated by
*       the inverse transform method.
 
        Nc = nsize/2+1
        fs0 = 0.033*pi*lcn2*kk*kk*2.0*dkap*dkap
 
        do 100 jy=1,nsize
          do 110 jx=1,nsize
 
            kap = dkap*sqrt( float((jx-Nc)*(jx-Nc)+(jy-Nc)*(jy-Nc)) )
* kap**(-11/3) is undefined at kap=0.
            kap = max(kap,dkap)
 
* introducing smooth cut-offs below kapmin and beyond kapmax.
            atten = 1.0
            if(kap .gt. kapmax) then
                atten = exp(-(kap-kapmax)/(5.0*dkap))
            elseif(kap .lt. kapmin) then 
                atten = exp(-(kapmin-kap)/(5.0*dkap))
            endif
 

C Generation of random phase screen:
10          call jlp_random(u) 
            if ( u.le.0.0 )goto 10

            radi = sqrt(-2.0*log(u))
            fs = fs0*kap**(-11.0/3.0)
            work = atten*sqrt(fs)*radi

            call jlp_random(wrand)
C psc(jx,jy) = atten*sqrt(fs)*radi*exp(i*2.*pi*wrand)

            wrand = wrand * twopi
            pha_re(jx,jy) = work*cos(wrand)
            pha_im(jx,jy) = work*sin(wrand)
 
110         continue
100       continue
 
*---------- Fourier transform ------------------------------
        write(*,555) nsize,nsize 
555     format('+', 'start ',i3,'x',i3,' sized fft              ')
        kod = 1
        call fft_2d(pha_re,pha_im,nsize,nsize,nsize,kod)
        write(*,666)
666     format('+', 'end fft for screen                         ')

C We are only interested in the phase term ,
C so we take arbitrary exp(i * phi) with phi = pha_re
C (We could have taken also pha_im)
        do jy=1,nsize
          do jx=1,nsize
            work = pha_re(jx,jy)
            pha_re(jx,jy) = cos(work)
            pha_im(jx,jy) = sin(work) 
          end do
        end do

        return
        end
********************************************************************** 
        include 'hrsa:fft_jlp.for'
