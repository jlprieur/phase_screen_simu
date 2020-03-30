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
*       For instance, the refractive effect can be seen by setting lo = 0.1,
*       while the diffractive effect can be investigated by setting L0 = 0.4.
*******************************************************************************
*------ declarations ---------*
*       constants
        real*4    pi, rapdeg, rapmin, rapsec
        complex iu
 
*       record size
        integer*4  RECS
        parameter( RECS = 512 )
 
*       integrated structure constant of the phase fluctuation
        real*4 LCn2
 
*       size of phase screen and resolution (m/pix)
        real*4 Ls, Lp
        parameter( Ls = 10.24 )
        parameter( Lp = 0.02  )
 
*       size of phase screen, size of power spectrum
        integer*4 Ns, Np
        parameter( Ns = 512  )
        parameter( Np = Ns/2 )
 
*       wave number, spatial frequencies and scale length of tubulence
        real k, kapmax, kapmin, L0, dkap, lambda, lo
 
*       complex phase screen, real and imag are independent.
        complex psc(Ns,Ns)
 
*       power spectrum array
        real*4 pwr(Np,Np)
 
*       input control parameter file and output power spectrum file
        character*30 ifile, ofile
 
*       I/O units
        integer*4 f_in, f_out
        parameter( f_in = 11, f_out = 12 )
 
*       # of frames, # of screens, # of frames per screen
        integer*4  nframe, nscreen, frm_per_scr
 
*       maximum # of components in source, maximum number of apertures
        integer*4 Maxcomp, Maxap
        parameter( Maxcomp = 50, Maxap = 15 )
 
*       source information (intensity and position in arcsec)
        real*4 s_int(Maxcomp), s_pos(Maxcomp,2)
 
*       aperture information (coordinate, radius, type)
        real*4 c_pos(Maxap,2), c_rad(Maxap)
        integer*4 c_typ(Maxap)
*
        common pi,rapdeg,rapmin,rapsec, iu
        common k, kapmax, kapmin, L0, LCn2, dkap, lambda
        common nframe, nscreen, frm_per_scr
        common x_min, x_max, xstep, y_min, y_max, ystep
        common ncomp, ncircle, nsamp, fraction
 
*       work area for subroutine power
        real*4  image1(Np,Np), image2(Np,Np),  ary1(Np,Np), ary2(Np,Np)
        complex pdata1(Np,Np), pdata2(Np,Np), cary1(Np,Np), cary2(Np,Np)
 
*       variables local to main program
        character*80 comment
        integer*4 seed, no_dist
 
*       common variables are not compatible with parameter statement.
        pi = 3.141592653
        rapdeg = pi/180.0
        rapmin = rapdeg/60.0
        rapsec = rapdeg/3600.0
        iu = (0.0,1.0)
        LCn2 = 5.0e-13
 
* JLP
        call jlp_begin
        call jlp_inquifmt
 
*       input file names
        write(*,10)
   10   format(1x,' input ascii file name : ',$)
        read(*,'(a)') ifile
        write(*,20)
   20   format(1x,'output power file name : ',$)
        read(*,'(a)') ofile
        write(*,25)
   25   format(1x)
 
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
            backspace f_in
            read(f_in,30,end=100) comment
            write(*,35) comment
          end if
C Keyword "Turbulence", parameters: outer and inner  scales
          if( comment(1:1) .eq. 'T' .or. comment(1:1) .eq. 't') then
            read(f_in,*) L0, lo
            backspace f_in
            read(f_in,30,end=100) comment
            write(*,35) comment
          end if
C Keyword "Wavelength", parameters: central wavel., fractional bandwidth, 
C sampling points
          if( comment(1:1) .eq. 'W' .or. comment(1:1) .eq. 'w') then
            read(f_in,*) lambda, fraction, nsamp
            backspace f_in
            read(f_in,30,end=100) comment
            write(*,35) comment
          end if
C Keyword "Source", parameters: number of points, followed by
C Intensity1, X1, Y1 
C Intensity2, X2, Y2... 
          if( comment(1:1) .eq. 'S' .or. comment(1:1) .eq. 's') then
            read(f_in,*) ncomp
            backspace f_in
            read(f_in,30,end=100) comment
            write(*,35) comment
            do n=1,ncomp
              read(f_in,*) s_int(n),s_pos(n,1),s_pos(n,2)
              backspace f_in
              read(f_in,30,end=100) comment
              write(*,35) comment
            end do
          end if
C Keyword "Aperture", parameters: number of apertures, followed by
C Xcenter, Ycenter, radius, type (0 or 1)
          if( comment(1:1) .eq. 'A' .or. comment(1:1) .eq. 'a') then
            read(f_in,*) ncircle
              backspace f_in
              read(f_in,30,end=100) comment
              write(*,35) comment
            do n=1,ncircle
              read(f_in,*) c_pos(n,1),c_pos(n,2),c_rad(n),c_typ(n)
              backspace f_in
              read(f_in,30,end=100) comment
              write(*,35) comment
            end do
          end if
C Keyword "Phase screen", parameters: xmin, xmax, step 
          if( comment(1:1) .eq. 'P' .or. comment(1:1) .eq. 'p') then
            read(f_in,*) x_min, x_max, xstep
            backspace f_in
            read(f_in,30,end=100) comment
            write(*,35) comment
            read(f_in,*) y_min, y_max, ystep
            backspace f_in
            read(f_in,30,end=100) comment
            write(*,35) comment
          end if
C Keyword "Frame", parameters: nber of frames 
          if( comment(1:1) .eq. 'F' .or. comment(1:1) .eq. 'f') then
            read(f_in,*) nframe
            backspace f_in
            read(f_in,30,end=100) comment
            write(*,35) comment
          end if
C Keyword "No disturbance" (no atmosphere)
          if( comment(1:1) .eq. 'N' .or. comment(1:1) .eq. 'n') then
            no_dist = 1
          end if
 
   40   continue
  100   continue
 
*       k: wave number, kap: spatial frequency in aperture plane
 
        k = 2.0 * pi / lambda
 
        if (Ls .lt. L0) then
          kapmin = 2.0 * pi / Ls
        else
          kapmin = 2.0 * pi / L0
        end if
 
        if (lo .lt. Lp) then
           kapmax = 2.0 * pi / Lp
        else
           kapmax = 2.0 * pi / lo
        end if
 
*       dkap: one pixel in spatial frequency in aperture plane
        dkap = 2.0 * pi / Ls
*
        frm_per_scr =
     1  (int((x_max-x_min)/xstep)+1)*(int((y_max-y_min)/ystep)+1)
        print *, 'frm_per_scr ',frm_per_scr
*
        if( mod(nframe,frm_per_scr) .eq. 0 ) then
          nscreen = nframe / frm_per_scr
        else
          nscreen = nframe / frm_per_scr + 1
        end if
 
*       Two real phase screens per complex random number generation in Fourier
*       domain.
 
        if( mod(nscreen,2) .eq. 0 ) then
          ncall = nscreen/2
        else
          ncall = nscreen/2 + 1
        end if
 
        do 300 i = 1, ncall
*------------------------------------
        if ( no_dist .eq. 1 ) then
*       no atmospheric disturbance
          do iz=1, Ns
            do jz=1, Ns
              psc(iz,jz) = (0.0,0.0)
            end do
          end do
        else
*       create a complex phase screen
        write(*,310) i
  310   format('+', 'generating ',i3,'th complex phase screen')
        call  screen(psc, Ns, seed)
        end if
 
*       sample phase screen, form images, calculate power
 
        call  power(psc, Ns, pwr, Np, iframe, Maxcomp, Maxap,
     1     s_int, s_pos, c_pos, c_rad, c_typ, pdata1, pdata2,
     2             image1, image2, ary1, ary2, cary1, cary2)
 
          if(iframe .ge. nframe) goto 400
*
*------------------------------------
  300   continue
*
 
*       normalize integrated power spectrum
  400   zero = pwr(Np/2+1,Np/2+1)
        do lx=1, Np
          do ly=1, Np
            pwr(lx,ly) = pwr(lx,ly)/zero
          end do
        end do
 
*       create a figaro file
 
        No = Np*Np
        call jlp_writeimag(pwr,Np,Np,Np,ofile,ocomments)
*
        call jlp_end
        stop
        end
 
 
        subroutine screen(psc, Ns, seed)
***************************************************************************
*       creates complex phase screen
***************************************************************************
        complex psc(Ns,Ns)
        real*4 k, kapmax, kapmin, L0, LCn2, dkap, u, work
        integer*4 frm_per_scr
        complex iu
        common pi,rapdeg,rapmin,rapsec, iu
        common k, kapmax, kapmin, L0, LCn2, dkap, lambda
        common nframe, nscreen, frm_per_scr
        common x_min, x_max, xstep, y_min, y_max, ystep
        common ncomp, ncircle, nsamp, fraction
 
*       radial spatial frequency in aperture plane
        real*4 kap
 
*       seed for the random number generator, dimension array for fft
*       origin of the phase screen
        integer*4 seed, N(2), Nc
 
*------ generation of random phase in spatial domain -------------
*       Complex Gaussian random numbers are generated by
*       the inverse transform method.
 
        N(1) = Ns
        N(2) = Ns
        Nc = Ns/2+1
 
        do 100 jx=1, Ns
          do 110 jy=1, Ns
 
            kap = dkap*sqrt( float( (jx-Nc)**2+(jy-Nc)**2 ) )
 
*           introducing smooth cut-offs below kapmin and beyond kapmax.
 
            if(kap .gt. kapmax) atten = exp(-(kap-kapmax)/(5.0*dkap))
            if(kap .lt. kapmin) atten = exp(-(kapmin-kap)/(5.0*dkap))
            if(kap.le.kapmax .and. kap.ge.kapmin) atten = 1.0
 
*           kap**(-11/3) is undefined at kap=0.
            if(kap .eq. 0.0) kap = dkap
 
            call jlp_random_init(seed)
10          call jlp_random(u) 
            if ( u.gt.0.0 ) then
              radi = sqrt(-2.0*log(u))
            else
              goto 10
            end if
            Fs = 0.033*pi*LCn2*k*k*2.0*kap**(-11.0/3.0)*dkap*dkap
            call jlp_random(work)
            psc(jx,jy) = atten*sqrt(Fs)*radi*cexp(iu*2.0*pi*work)
 
  110     continue
  100   continue
 
*---------- Fourier transform ------------------------------
*       shift origin of Fourier space to (Nc,Nc)
        do i=1,Ns,2
          do j=1,Ns,2
            psc(i+1,j) = -psc(i+1,j)
            psc(i,j+1) = -psc(i,j+1)
          end do
        end do
 
        write(*,555) N(1),N(2)
  555   format('+', 'start ',i3,'x',i3,' sized fft              ')
        call fourn1(psc,N,2,-1,+1)
        write(*,666)
  666   format('+', 'end fft for screen                         ')
 
*       now psc is a complex random phase screen
*       to avoid discontinuity shift origin of phase screen to (Nc,Nc)
        do i=1,Ns,2
          do j=1,Ns,2
            psc(i+1,j) = -psc(i+1,j)
            psc(i,j+1) = -psc(i,j+1)
          end do
        end do
 
        return
        end
 
 
        subroutine power(psc, Ns, pwr, Np, iframe, Maxcomp, Maxap,
     1          s_int, s_pos, c_pos, c_rad, c_typ, pdata1, pdata2,
     2                   image1, image2, ary1, ary2, cary1, cary2)
*******************************************************************************
*       produces fringe/speckle patterns of multiple point sources.
*       Finite bandwidth effect is taken into account.
*******************************************************************************
*       source information (intensity and position in arcsec)
        real*4 s_int(Maxcomp), s_pos(Maxcomp,2)
*       aperture information (coordinate, radius, type)
        real*4 c_pos(Maxap,2), c_rad(Maxap)
        integer*4 c_typ(Maxap)
        real*4  lpix
        complex psc(Ns,Ns)
        real*4 k, kapmax, kapmin, Lo, LCn2, dkap
        integer*4 frm_per_scr
        complex iu
 
        common pi,rapdeg,rapmin,rapsec,iu
        common k, kapmax, kapmin, L0, LCn2, dkap, lambda
        common nframe, nscreen, frm_per_scr
        common x_min, x_max, xstep, y_min, y_max, ystep
        common ncomp, ncircle, nsamp, fraction
*
        integer*4 N(2)
 
*       power spectrum
        real*4 pwr(Np,Np)
 
*       temporary arrays
        real*4  image1(Np,Np), image2(Np,Np),  ary1(Np,Np), ary2(Np,Np)
        complex pdata1(Np,Np), pdata2(Np,Np), cary1(Np,Np), cary2(Np,Np)
 
 
        N(1) = Np
        N(2) = Np
        lpix = 10.24/Ns
        Wc = Np/2+1
        Scale = Ls/lambda*Np/Ns*rapsec
 
*       setting apertures on Np sized window in Ns sized complex phase screen
 
        do 10 x=x_min,x_max,xstep
        do 20 y=y_min,y_max,ystep
**************************** major loop ***************************************
 
        write(*,111) x,y
  111   format('+',' creating a pair of speckle/fringe patterns',1x,
     1         'xy',f10.3,1x,f10.3)
 
        write(*,333) ncircle, nsamp, c_typ(1), c_typ(2)
  333   format(1x,'ncircle nsamp typ(1) typ(2) ',1x,4(i2,1x) )
 
 
        do 100 loop1=1,nsamp
**************** intermediate loop *****************
*        nsamp data points within the bandpass
         factor =  1.0 + fraction*(loop1 - nsamp/2)
 
         npx = 0
 
          do 30 i=1,ncircle
*-----------------------------------------
*       filling pdata with aperture phases
 
            do 40 xa = x+c_pos(i,1)-c_rad(i), x+c_pos(i,1)+c_rad(i), lpix
            do 50 ya = y+c_pos(i,2)-c_rad(i), y+c_pos(i,2)+c_rad(i), lpix
              ox = xa-x-c_pos(i,1)
              oy = ya-y-c_pos(i,2)
              if( ox*ox+oy*oy .lt. c_rad(i)*c_rad(i) ) then
                ix = nint( (xa-x)/lpix )+Wc
                iy = nint( (ya-y)/lpix )+Wc
 
*       shift origin of the window to (Wc,Wc)
                lx = nint(xa/lpix)
                ly = nint(ya/lpix)
 
                if( c_typ(i) .eq. 1) npx = npx+1
                if( c_typ(i) .eq. 0) npx = npx-1
 
                if( mod(npx,100) .eq. 0 ) then
                  write(*,777) npx
                end if
  777   format('+', '# of  positive pixels',1x,i5)
 
                if( mod(ix+iy,2) .eq. 0 ) then
                  pdata1(ix,iy) = cexp( iu*factor*real(psc(lx,ly)) )*c_typ(i)
                  pdata2(ix,iy) = cexp(iu*factor*aimag(psc(lx,ly)) )*c_typ(i)
                else
                  pdata1(ix,iy) =-cexp( iu*factor*real(psc(lx,ly)) )*c_typ(i)
                  pdata2(ix,iy) =-cexp(iu*factor*aimag(psc(lx,ly)) )*c_typ(i)
                end if
 
              end if
   50       continue
   40       continue
*-----------------------------------------
   30     continue
 
*       simulate propagation from aperture plane to detector plane
 
        call fourn1(pdata1,N,2,-1,+1)
        call fourn1(pdata2,N,2,-1,+1)
 
*       form speckle/fringe patterns
        do ix=1, Np
          do iy=1, Np
            image1(ix,iy) = image1(ix,iy) +
     1                      real(pdata1(ix,iy)*conjg(pdata1(ix,iy)))
            image2(ix,iy) = image2(ix,iy) +
     2                      real(pdata2(ix,iy)*conjg(pdata2(ix,iy)))
          end do
        end do
 
*       reset pdata1 and pdata2
        do ix=1, Np
          do iy=1, Np
            pdata1(ix,iy) = (0.0,0.0)
            pdata2(ix,iy) = (0.0,0.0)
          end do
        end do
 
******************* end of loop1 ********************
  100   continue
 
        do 60 loop2=1,ncomp
*--------------------------------------------
*       ncomp components in source structure
          dx = int(Scale * s_pos(loop2,1))
          dy = int(Scale * s_pos(loop2,2))
 
        if( dx .ge. 0 .and. dy. ge. 0) then
          do ix=1,Np-dx
          do iy=1,Np-dy
            ary1(ix+dx,iy+dy) = ary1(ix,iy)+s_int(loop2)*image1(ix,iy)
            ary2(ix+dx,iy+dy) = ary2(ix,iy)+s_int(loop2)*image2(ix,iy)
          end do
          end do
        end if
 
        if( dx .ge. 0 .and. dy. lt. 0) then
          do ix=1,Np-dx
          do iy=-dy+1,Np
            ary1(ix+dx,iy+dy) = ary1(ix,iy)+s_int(loop2)*image1(ix,iy)
            ary2(ix+dx,iy+dy) = ary2(ix,iy)+s_int(loop2)*image2(ix,iy)
          end do
          end do
        end if
 
        if( dx.lt.0 .and. dy.ge.0) then
          do ix=-dx+1,Np
          do iy=1,Np-dy
            ary1(ix+dx,iy+dy) = ary1(ix,iy)+s_int(loop2)*image1(ix,iy)
            ary2(ix+dx,iy+dy) = ary2(ix,iy)+s_int(loop2)*image2(ix,iy)
          end do
          end do
        end if
 
        if( dx.lt.0 .and. dy.lt.0) then
          do ix=-dx+1,Np
          do iy=-dy+1,Np
            ary1(ix+dx,iy+dy) = ary1(ix,iy)+s_int(loop2)*image1(ix,iy)
            ary2(ix+dx,iy+dy) = ary2(ix,iy)+s_int(loop2)*image2(ix,iy)
          end do
          end do
        end if
*-------------------------------------------
   60   continue
 
 
        write(*,320) iframe+1
 
*       shift origin of the image to (Wc,Wc) before fourier transform
        do ix=1,Np,2
        do iy=1,Np,2
          cary1(ix,iy)     = ary1(ix,iy)
          cary1(ix+1,iy+1) = ary1(ix+1,iy+1)
          cary1(ix+1,iy)   =-ary1(ix+1,iy)
          cary1(ix,iy+1)   =-ary1(ix,iy+1)
        end do
        end do
 
        call fourn1(cary1,N,2,-1,+1)
 
        do jx=1,Np
        do jy=1,Np
          pwr(jx,jy) = pwr(jx,jy) + real(cary1(jx,jy)*conjg(cary1(jx,jy)))
        end do
        end do
 
        iframe = iframe + 1
        if( iframe .ge. nframe) return
 
        write(*,320) iframe+1
  320   format('+', 'processing ',i4,' th frame                    ')
 
*       shift origin of the image to (Wc,Wc) before fourier transform
        do ix=1,Np,2
        do iy=1,Np,2
          cary2(ix,iy)     = ary2(ix,iy)
          cary2(ix+1,iy+1) = ary2(ix+1,iy+1)
          cary2(ix+1,iy)   =-ary2(ix+1,iy)
          cary2(ix,iy+1)   =-ary2(ix,iy+1)
        end do
        end do
 
        call fourn1(cary2,N,2,-1,+1)
 
        do jx=1,Np
        do jy=1,Np
          pwr(jx,jy) = pwr(jx,jy) + real(cary2(jx,jy)*conjg(cary2(jx,jy)))
        end do
        end do
 
        iframe = iframe + 1
        if( iframe .ge. nframe) return
 
*       reset temporary arrays
        do jx=1,Np
          do jy=1,Np
            image1(jx,jy) = 0.0
            image2(jx,jy) = 0.0
              ary1(jx,jy) = 0.0
              ary2(jx,jy) = 0.0
             cary1(jx,jy) = (0.0,0.0)
             cary2(jx,jy) = (0.0,0.0)
          end do
        end do
 
************************* end of major loop ***********************************
   20   continue
   10   continue
 
        return
        end
********************************************************************** 
        include 'hrsa:fft_jlp.for'
