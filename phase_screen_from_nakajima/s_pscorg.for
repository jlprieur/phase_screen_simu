        program s_pscorg
*******************************************************************************
*
* s_pscorg is a phase screen simulator which produces artificial fringe
* or speckle patterns and obtain integrated power spectrum or(/and)
* optionally specified bispectral components.
*
* The atmospheric disturbance is given by the Kolmogorov spectrum.
* Seeing condition, aperture configuration and the object ( at this
* point, only multiple point sources ) and the usage of the phase screen
* must be specified by the user in an ascii input file.
*
* Following is a sample input file. (free format, ! indicates a
* comment line.)
*
* ! Field of view: side of image in arcseconds
* Field of view
* 12
* ! seed for random number generator (large odd integer)
* Random
* 7597551
* ! outer scale, inner scale of the turbulence(m) - recommended values
* Turbulence scale
* 10.0  0.01
* ! Fried radius (r_zero) in m 
* Fried 
* 0.1
* ! central wavelength(m), fractional bandwidth, sampling points
* Wavelength
* 0.65e-6 0.016 3
* ! source: # of point sources
* !         intensity(arbitrary unit), position offset(x,y) in arcsecond
* Source
* 2
* 1.0  -0.7 -0.6
* 0.5   0.8 0.4
* ! # of circular apertures (including voids)
* ! configuraton: center of aperture(x,y), radius(m) and type
* ! 1 is an aperture,0 is a void. Voids should follow ordinary apertures.
* ! Example: annular mask
* Aperture
* 2
* 0.0 0.0 2.0  1
* 0.0 0.0 1.9  0
* ! Multiple subsets are taken out of one phase screen.
* ! Equally spaced sampling is specified by
* ! x_min, x_max, step size in x(m)
* ! y_min, y_max, step size in y(m)
* ! JLP/WARNING: x_min, y_min should be larger than the aperture radii !!!
* Phase_screen
* 1.25 8.75  0.14
* 1.25 8.75  2.5
* ! total number of frames to be processed
* Frames
* 100
* ------ as an option `no atmosphere' can be specified by ------
* ! obtaining beam pattern
* No_disturbance
* --------------------------------------------------------------
*
* This is a re-organized version of the pscsn.for developed on Tacos.
* The program is planned to be as users oriented as possible for the
* future development of Optical/IR interferometry.
*
* version 1.0    Dec 25, 1988 (Merry Christmas!)     Tadashi Nakajima
* Today's resolution : `I will have real holidays next year and live
*                       a human life.'
* version 2.0    Dec 30, 1988   Now sharp cut-offs of fluctuation power
* spectrum at kapmin and kapmax are introduced to investigatate the
* relative significance of the diffractive effect and refractive effect.
* For instance, the refractive effect can be seen by setting lo = 0.1,
* while the diffractive effect can be investigated by setting L0 = 0.4.
*
* JLP: 
*   - I introduce fourn1 and jlp_random
*   - I remove the common blocks, split the program into subroutines
*      and complete the documentation 
*   - AND I CORRECT DEEPLY THIS PROGRAM SINCE IT RELIED ON AN UNKNOWN
*      STRUCTURE OF FOURIER COMPLEX ARRAYS AND FOURIER TRANSFORM!
*   - I add r_zero and field of view into the parameter file
*   - I change the screen generator (jlp_phase_screen) so that the 
*   simulation is compatible with the r_zero value (the old version was not)
*   - I replace the pseudo-convolution of the objet with the PSF by a true 
*   convolution
*   - I implement a good formula for the image smearing due to finite bandwidth
*   (the original was wrong!) 
*
* JLP
* Version 26/05/2008
*******************************************************************************
*------ declarations ---------*
        implicit none
* Constants
        real*4 pi
 
* Integrated structure constant of the phase fluctuation
* LCn2 = \int Cn2(h) dh
* And also: 
* r_0 = (0.423 kk^2 LCn2)^{-3/5}
        real*4 LCn2, r_zero
 
* Size of phase screen and resolution (m/pix)
        real*4 Ls, Lp
 
* Size of phase screen, size of power spectrum
        integer*4 Ns, Np
        parameter(Ns = 512)
        parameter(Np = Ns/2)
 
* Wave number, scale length of turbulence
* fov: field of view (for the simulated images)
        real kk, L0, lambda, lo, fov
 
* psc: complex phase screen, real and imag are independent.
* FT_object(Np,Np): Fourier Transform of the object
        complex psc(Ns,Ns), FT_object(Np,Np)
 
* Power spectrum array
        real*4 pwr(Np,Np)
 
* Input control parameter file and output power spectrum file
        character ifile*60, ofile*60, ocomments*80
 
* I/O units
        integer*4 f_out
        parameter(f_out = 12)
 
* nframes: number of frames
* nscreens: number of phase screens 
* frm_per_scr: number of frames per phase screen
        integer*4  nframes, nscreens, frm_per_scr
 
* Maxcomp: maximum number of components in source 
* Maxap: maximum number of apertures
        integer*4 Maxcomp, Maxap
        parameter(Maxcomp = 50, Maxap = 15)
 
* Source information (intensity and position in arcsec)
        real*4 s_int(Maxcomp), s_pos(Maxcomp,2), object(Np,Np)
 
* Aperture information (coordinate, radius, type)
        real*4 c_pos(Maxap,2), c_rad(Maxap)
        integer*4 c_typ(Maxap)
* Miscellaneous:
        real*4 zero, x_min, x_max, xstep, y_min, y_max, ystep
        real*4 rapsec, fractional_bw
        integer*4 ncalls, ncircle, nsamp, lx, ly, i, iz, jz, iframe
        integer*4 ncomp, seed, no_dist
 
* Initialization of variables:
        pi = 3.141592653

* rapsec: radians per arcsecond
        rapsec = pi/(180. * 3600.0)
 
* JLP
* Program can then be run as "runs s_pscorg s_screen.input power"
        call jlp_begin
        call jlp_inquifmt
 
* Input file names (from terminal or "jlp_lu5.tmp" file when run with "runs")
        write(6,10)
   10   format(1x,'Input ascii file name : ',$)
        read(5,'(a)') ifile
        write(6,20)
   20   format(1x,'Output power FITS file name : ',$)
        read(5,'(a)') ofile
        write(6,25)
   25   format(1x)
 
C Read input parameter file:
        call read_input_file(ifile, seed, L0, lo, lambda, fractional_bw,
     1                nsamp, ncomp, s_int, s_pos, ncircle, c_pos, 
     1                c_rad, c_typ, x_min, x_max, xstep, y_min, 
     1                y_max, ystep, fov, r_zero, nframes, no_dist, 
     1                Maxcomp, Maxap)

C Initialization of random generator:
        call jlp_random_init(seed)

*  kk: wave number 
        kk = 2.0 * pi / lambda
 
* Original version with LCn2 = 5.0e-13
* r_0 = (0.423 kk^2 LCn2)^{-3/5}
        LCn2 = (r_zero**(3./5.)) / (0.423 * kk * kk) 
        write(6,287) lambda, r_zero, LCn2
287     format('lambda=',e12.5,' r_zero=',f8.3,' LCn2=',e12.5)

* Scale of image in pixels/arcsecond: Scale = (Ls/lambda)*(Np/Ns)*rapsec
* Field of view = Np/Scale
* Size of phase screen and resolution (m/pix)
        Ls = Ns / ((fov / lambda) * rapsec)
        write(6,*) 'Size of pupil: Ls = ',Ls,' m'
        Lp = Ls / real(Ns)
        write(6,*) 'Resolution of pupil: Lp = ',Lp,' m'

C Create an image of the source from ncomp, s_int, s_pos
         call create_object(object, FT_object, Np, Ns, ncomp, 
     1              s_int, s_pos, Ls, lambda, Maxcomp)

* Number of frames per phase screen 
        frm_per_scr =
     1  (int((x_max-x_min)/xstep)+1)*(int((y_max-y_min)/ystep)+1)
*
        if( mod(nframes,frm_per_scr) .eq. 0 ) then
          nscreens = nframes / frm_per_scr
        else
          nscreens = nframes / frm_per_scr + 1
        end if
 
        do iz=1,Np
        do jz=1,Np
          pwr(jz,jz) = 0. 
        end do
        end do
* Two real phase screens per complex random number generation in Fourier
* domain.
 
        if( mod(nscreens,2) .eq. 0 ) then
          ncalls = nscreens/2
        else
          ncalls = nscreens/2 + 1
        end if
        write(6,*)' ******** Will now simulate ', nframes,
     1            ' elementary frames ************'
        write(6,*)' Number of frames per phase screen: ', frm_per_scr 
        write(6,*)' Number of phase screen needed: ', nscreens
        write(6,*)' Number of calls of phase_screen routine needed: ', 
     1            ncalls
 
        iframe = 0
        do 300 i = 1, ncalls
*------------------------------------
        if (no_dist.eq.1) then
* no atmospheric disturbance
        write(6,312) i, iframe+1
  312   format('+', 'generating ',i3,
     1         'th zero phase screen (iframe=',i3,')')
          do iz=1, Ns
            do jz=1, Ns
              psc(iz,jz) = (0.0,0.0)
            end do
          end do
        else
* Create a complex phase screen
        write(6,310) i, iframe+1
  310   format('+', 'generating ',i3,
     1         'th complex phase screen (iframe=',i3,')')
C New version:
        call jlp_phase_screen(psc, Ns, seed, L0, r_zero) 
C Old version:
C        call  phase_screen(psc, Ns, seed, Ls, L0, lo, 
C     1                     kk, LCn2)
        end if
 
* sample phase screen, form images, calculate power
        call  power(psc, Ns, FT_object, pwr, Np, iframe, nframes, Ls,
     1                  lambda, fractional_bw, nsamp, ncomp,
     1                  s_int, s_pos, ncircle, c_pos, c_rad, c_typ,
     1                  x_min, x_max, xstep, y_min, y_max, ystep,
     1                  Maxcomp, Maxap)
 
          if(iframe .ge. nframes) goto 400
*
*------------------------------------
* End of loop
  300   continue
*

* Exit from loop
  400   continue

        call jlp_recent_real(pwr, Np)
* Normalize integrated power spectrum
        zero = pwr(Np/2+1,Np/2+1)
        do lx=1, Np
          do ly=1, Np
            pwr(lx,ly) = pwr(lx,ly)/zero
          end do
        end do
 
* Create a FITS file
        write(6,*)' Writing output file:',ofile
        write(ocomments,88)ifile,char(0)
88      format('From ',A,A)
        call jlp_writeimag(pwr,Np,Np,Np,ofile,ocomments)
*
        call jlp_end
        stop
        end
***************************************************************************
* Creates complex phase screen (JLP's version)
*
* INPUT:
* Ns: size of phase screen (512)
* seed: seed for the random number generator, 
*
* OUPUT:
*
***************************************************************************
        subroutine jlp_phase_screen(psc, Ns, seed, L0, r_zero) 
        implicit none
        integer*4 Ns, seed
        complex psc(Ns,Ns),iiu
        real*4 L0, r_zero, pi, work(Ns,Ns), mean, sigma
        real*4 rad2
        real*8 sum, sumsq
        logical debug
        character out_file*60, out_comments*80
        real*4 urad, radi, Fs, Fs0, utheta 
 
* N(): dimension array for fft
* Ns_cent: origin of the phase screen
        integer*4 N(2), Ns_cent, jx, jy, i, j
 
        debug=.true.
* Constants:
* iiu = sqrt(-1) or "i"
        iiu = (0.0,1.0)
        pi = 3.141592653

* r_0 = (0.423 kk^2 LCn2)^{-3/5}
* Hence Fs0 is proportional to kk^2
C JLP2008 / I put Ns^4 instead of Ns^2 ...
        Fs0 = 0.023 * Ns * Ns * Ns * Ns * (L0 / r_zero)**(5./3.)

*------ generation of random phase in spatial domain -------------
* Complex Gaussian random numbers are generated by
* the inverse transform method.
 
        N(1) = Ns
        N(2) = Ns
        Ns_cent = Ns/2+1
 
        do 100 jx=1, Ns
          do 110 jy=1, Ns
 
            rad2=(jx-Ns_cent)*(jx-Ns_cent)+(jy-Ns_cent)*(Jy-Ns_cent)
 
* rad**(-11/3) is undefined at rad2=0.
            if(rad2 .eq. 0.0) rad2 = 0.5 
 
10          call jlp_random(urad) 
            if ( urad.gt.0.0 ) then
              radi = sqrt(-2.0*log(urad))
            else
              goto 10
            end if
            Fs = Fs0 * rad2**(-11.0/6.0)
            call jlp_random(utheta)
* iiu = sqrt(-1) or "i"
* Hence psc is proportional to kk = 2 pi / lambda
            psc(jx,jy) = sqrt(Fs)*radi*cexp(iiu*2.0*pi*utheta)
  110     continue
  100   continue
 
*---------- Fourier transform ------------------------------
* Shift origin of Fourier space to (0,0) before inverse Fourier transform
* Old recentering version from original program was wrong!
        call jlp_recent_cmplx(psc, Ns)

        write(6,555) N(1),N(2)
  555   format('+', 'start ',i3,'x',i3,' sized fft              ')
C        call fourn1(psc,N,2,-1,+1)
C JLP2008 (by comparing with s_psc.for):
C fourn1(data,nn,ndim,isign)
C SHOULD USE MY fourn1 routine (in fft_jlp.for)
C NOT s_fourn1.for !!!
        call fourn1(psc,N,2,-1)
        write(6,666)
  666   format('+', 'end fft for screen                         ')
 
* Now psc is a complex random phase screen
* to avoid discontinuity shift origin of phase screen to (Ns/2,Ns/2)
* Old recentering version from original program was wrong!
        call jlp_recent_cmplx(psc, Ns)

C Output phase screen:
        if(debug)then
          sum = 0.
          sumsq = 0.
          do i=1, Ns
            do j=1, Ns
              work(i, j) = real(psc(i,j))
              sum = sum + work(i,j)
              sumsq = sumsq + work(i,j) * work(i,j)
            end do
          end do
          mean = sum / real(Ns * Ns)
          sigma = sqrt(sumsq / real(Ns * Ns) - mean * mean)
          write(6,*)'Phase screen: mean=', mean,' sigma=',sigma
          out_file='screen_real_psc'
          out_comments='Real part of phase screen'
          call jlp_writeimag(work,Ns,Ns,Ns,out_file,out_comments)
        endif

        return
        end
***************************************************************************
* Creates complex phase screen (original version)
*
* INPUT:
* Ns: size of phase screen (512)
* seed: seed for the random number generator, 
* kk: wave number
* LCn2 integrated CN2 over all the altitudes (originally: LCn2=5.0e-13)
*
* OUPUT:
*
***************************************************************************
        subroutine phase_screen(psc, Ns, seed, Ls, L0, lo, 
     1                          kk, LCn2)
        implicit none
        integer*4 Ns, seed
        complex psc(Ns,Ns),iiu
        real*4 dkap, kapmin, kapmax, kk, LCn2, pi, Ls, L0, lo, Lp
        real*4 work(Ns,Ns), mean, sigma, ww
        real*8 sum,sumsq
        logical debug
        character out_file*60, out_comments*80
 
* kap: radial spatial frequency in aperture plane ( * 2 pi / L)
        real*4 kap, atten, urad, radi, Fs, Fs0, utheta 
 
* N(): dimension array for fft
* Ns_cent: origin of the phase screen
        integer*4 N(2), Ns_cent, jx, jy, i, j
 
        debug=.true.
* Constants:
* iiu = sqrt(-1) or "i"
        iiu = (0.0,1.0)
        pi = 3.141592653

* Lp: pixel scale of phase screen
        Lp = Ls/Ns

* kapmin, kapmax: min and max spatial frequencies in aperture plane
* Ls: size of phase screen (in m)
* Lp: resolution, i.e., size of one pixel in pupil (in m/pix)
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
 
* dkap: one pixel in spatial frequency in aperture plane (dkap = 2 pi / Ls)
        dkap = 2.0 * pi / Ls

* dkap: one pixel in spatial frequency in aperture plane (dkap = 2 pi / Ls)
* r_0 = (0.423 kk^2 LCn2)^{-3/5}
        Fs0 = 0.033*pi*LCn2*kk*kk*2.0*dkap*dkap

*------ generation of random phase in spatial domain -------------
* Complex Gaussian random numbers are generated by
* the inverse transform method.
 
        N(1) = Ns
        N(2) = Ns
        Ns_cent = Ns/2+1
 
        do 100 jx=1, Ns
          do 110 jy=1, Ns
 
            kap = dkap*sqrt( float( (jx-Ns_cent)**2+(jy-Ns_cent)**2 ) )
 
* Introducing smooth cut-offs below kapmin and beyond kapmax.
 
            if(kap .gt. kapmax) atten = exp(-(kap-kapmax)/(5.0*dkap))
            if(kap .lt. kapmin) atten = exp(-(kapmin-kap)/(5.0*dkap))
            if(kap.le.kapmax .and. kap.ge.kapmin) atten = 1.0
 
* kap**(-11/3) is undefined at kap=0.
            if(kap .eq. 0.0) kap = dkap
 
10          call jlp_random(urad) 
            if ( urad.gt.0.0 ) then
              radi = sqrt(-2.0*log(urad))
            else
              goto 10
            end if
* Fs0 = 0.033*pi*LCn2*k*k*2.0*dkap*dkap
* Fs = 0.033*pi*LCn2*kk*kk*2.0*kap**(-11.0/3.0)*dkap*dkap
            Fs = Fs0 * kap**(-11.0/3.0)
            call jlp_random(utheta)
* iiu = sqrt(-1) or "i"
            psc(jx,jy) = atten*sqrt(Fs)*radi*cexp(iiu*2.0*pi*utheta)
  110     continue
  100   continue
 
*---------- Fourier transform ------------------------------
* Shift origin of Fourier space to (0,0)
* Old recentering version from original program was wrong!
        call jlp_recent_cmplx(psc, Ns)

        write(6,555) N(1),N(2)
  555   format('+', 'start ',i3,'x',i3,' sized fft              ')
C        call fourn1(psc,N,2,-1,+1)
C JLP2008 (by comparing with s_psc.for):
C fourn1(data,nn,ndim,isign)
C SHOULD USE MY fourn1 routine (in fft_jlp.for)
C NOT s_fourn1.for !!!
        call fourn1(psc,N,2,-1)
        write(6,666)
  666   format('+', 'end fft for screen                         ')
 
C JLP2008: I multiply with Ns*Ns, since otherwise
C the values of psc are much too small!
        ww = real(Ns * Ns)
          do i=1, Ns
            do j=1, Ns
              psc(i,j) = psc(i,j) * ww 
            end do
          end do

* Now psc is a complex random phase screen
* to avoid discontinuity shift origin of phase screen to (Ns/2,Ns/2)
* Old recentering version from original program was wrong!
        call jlp_recent_cmplx(psc, Ns)

C Output phase screen:
        if(debug)then
          sum = 0.
          sumsq = 0.
          do i=1, Ns
            do j=1, Ns
              work(i, j) = real(psc(i,j))
              sum = sum + work(i,j)
              sumsq = sumsq + work(i,j) * work(i,j)
            end do
          end do
          mean = sum / real(Ns * Ns)
          sigma = sqrt(sumsq / real(Ns * Ns) - mean * mean)
          write(6,*)'Phase screen: mean=', mean,' sigma=',sigma
          out_file='screen_real_psc'
          out_comments='Real part of phase screen'
          call jlp_writeimag(work,Ns,Ns,Ns,out_file,out_comments)
        endif

        return
        end
*******************************************************************************
* power
* Produces fringe/speckle patterns of multiple point sources.
* Finite bandwidth effect is taken into account.
*
* INPUT:
* psc(Ns,Ns): phase screen
* Ns: size of phase screen (=512)
* Np: size of power spectrum (=Ns/2)
* iframe: index of current frame to be processed
* nframes: number of frames to be simulated
* ncomp: number of components for the source
* s_int, s_pos: source information (intensity and position in arcsec)
* ncirc: number of circles for the aperture
* c_pos, c_rad, c_typ: aperture information (coordinate, radius, type)
* Ls: size of phase screen (in m)
* lambda, fractional_bw, nsamp: Wavelength parameters: central wavel., 
*                          fractional bandwidth and sampling points
* x_min, x_max, xstep : sampling parameters to extract a phase screen 
*                       from a larger one 
* y_min, y_max, ystep : sampling parameters to extract a phase screen 
*                       from a larger one 
* Maxcomp: maximum number of components in source 
* Maxap: maximum number of apertures
*
* OUTPUT:
* pwr(Np,Np): power spectrum
*******************************************************************************
        subroutine power(psc, Ns, FT_object, pwr, Np, iframe, 
     1                  nframes, Ls, 
     1                  lambda, fractional_bw, nsamp, ncomp,
     1                  s_int, s_pos, ncircle, c_pos, c_rad, c_typ,
     1                  x_min, x_max, xstep, y_min, y_max, ystep,
     1                  Maxcomp, Maxap)
        implicit none
        integer Ns, Np, Maxcomp, Maxap
        complex psc(Ns,Ns), FT_object(Np,Np)
        real*4 pwr(Np,Np), Ls, lambda, fractional_bw
        real*4 x_min, x_max, xstep, y_min, y_max, ystep
        integer*4 iframe, nframes, nsamp, ncomp, ncircle

* Aperture information (coordinate, radius, type)
        real*4 c_pos(Maxap,2), c_rad(Maxap)
        integer*4 c_typ(Maxap)

* Source information (intensity and position in arcsec)
        real*4 s_int(Maxcomp), s_pos(Maxcomp,2)
 
* Temporary arrays
        real*4  PSF1(Np,Np), PSF2(Np,Np)

C Miscellaneous:
        integer*4 N(2), Np_cent
        integer*4 output_screen
        real*4  Lp, Scale, x, y

* Constants
        real*4 pi, rapdeg, rapsec
        complex iiu
* iiu = sqrt(-1) or "i"
        iiu = (0.0,1.0)
        pi = 3.141592653

* rapsec: radians per arcsecond
        rapdeg = pi/180.0
        rapsec = rapdeg/3600.0
 
        output_screen = 0
        N(1) = Np
        N(2) = Np
* Ls: size of phase screen (in m)
* Lp: pixel scale of phase screen
        Lp = Ls/Ns
        Np_cent = Np/2+1
* lambda: wavelength (in meters)
* Scale: in pixels/arcsecond
        Scale = (Ls/lambda)*(real(Np)/real(Ns))*rapsec

* Setting apertures on Np sized window in Ns sized complex phase screen
 
        do 110 x=x_min,x_max,xstep
        do 120 y=y_min,y_max,ystep
**************************** major loop ***************************************
 
        write(6,111) x,y
  111   format('+',' creating a pair of speckle/fringe patterns at',
     1         'x=',f10.3,' y=',f10.3)
 
        write(6,333) ncircle, nsamp, c_typ(1), c_typ(2)
  333   format(1x,'ncircle nsamp typ(1) typ(2) ',1x,4(i2,1x) )
 
C Compute two elementary PSFs (PSF1 and PSF2)
        call compute_elementary_PSFs(psc, x, y, Ls,
     1                  lambda, fractional_bw, nsamp, 
     1                  ncircle, c_pos, c_rad, c_typ,
     1                  PSF1, PSF2, Ns, Np, Maxap, output_screen)

C Convolution of the object with the elementary PSFs 
        
        write(6,320) iframe+1
  320   format('+', 'processing ',i5,' th frame                    ')
       call process_elementary_image(pwr, PSF1, FT_object,
     1                      Np, output_screen)
        iframe = iframe + 1
        if( iframe .ge. nframes) return

        write(6,320) iframe+1
       call process_elementary_image(pwr, PSF2, FT_object,
     1                      Np, output_screen)
        iframe = iframe + 1
        if( iframe .ge. nframes) return
 
************************* end of major loop ***********************************
  120   continue
  110   continue

        return
        end
***************************************************************************
* read_input_file: read input parameter file
*
* INPUT:
* ifile: name of input parameter file 
*
* OUTPUT:
* seed: seed for random number generator (large odd integer)
* L0, lo: external and internal scale lengths of turbulence
* lambda, fractional_bw, nsamp: Wavelength parameters: central wavel., 
*                          fractional bandwidth and sampling points
* x_min, x_max, xstep : sampling parameters to extract a phase screen 
*                       from a larger one 
* y_min, y_max, ystep : sampling parameters to extract a phase screen 
*                       from a larger one 
* ncomp: number of components for the source
* s_int, s_pos: source information (intensity and position in arcsec)
* ncirc: number of circles for the aperture
* c_pos, c_rad, c_typ: aperture information (coordinate, radius, type)
* fov: field of view (in arcseconds)
* nframes: number of frames to be simulated
* no_dist: flag set to one if no atmospheric disturbance
* Maxcomp: maximum number of components in source 
* Maxap: maximum number of apertures
*
***************************************************************************
        subroutine read_input_file(ifile, seed, L0, lo, 
     1                  lambda, fractional_bw, nsamp, ncomp,
     1                  s_int, s_pos, ncircle, c_pos, c_rad, c_typ,
     1                  x_min, x_max, xstep, y_min, y_max, ystep,
     1                  fov, r_zero, nframes, no_dist, Maxcomp, Maxap)
        implicit none
        integer*4 Maxcomp, Maxap
        character ifile*60, comment*80
        real*4 L0, lo, lambda, fractional_bw
        real*4 s_int(Maxcomp), s_pos(Maxcomp,2)
        real*4 c_pos(Maxap,2), c_rad(Maxap)
        real*4 x_min, x_max, xstep, y_min, y_max, ystep, fov, r_zero
        integer*4 seed, nsamp, ncomp, ncircle, nframes, no_dist
        integer*4 c_typ(Maxap), j, n
* I/O units
        integer*4 f_in
        parameter(f_in = 11)
 
* Initialize no_dist to default value:
        no_dist = 0

* Read the input file
        open(unit=f_in, file=ifile, status='old')
 
   30   format(a80)
        do 40 j=1,100
 
          read(f_in,30,end=100) comment
          if( comment(1:1) .eq. '!') goto 40

          write(6,30) comment
 
C Keyword "Field of view", parameter: fov 
          if(comment(1:5) .eq. 'Field')then
            read(f_in,*) fov 
            write(6,*) 'OK: fov=',fov,' arcseconds'
C Keyword "Random", parameter: seed (random generator): 
          else if(comment(1:6) .eq. 'Random')then
            read(f_in,*) seed
            write(6,*) 'OK: seed=',seed
C Keyword "Turbulence scale", parameters: outer and inner  scales
          else if(comment(1:10) .eq. 'Turbulence')then
            read(f_in,*) L0, lo
            write(6,*) 'OK: L0=',L0,' lo=',lo
C Keyword "Fried radius", parameters: outer and inner  scales
          else if(comment(1:5) .eq. 'Fried')then
            read(f_in,*) r_zero
            write(6,*) 'OK: r0=',r_zero
C Keyword "Wavelength", parameters: central wavel., fractional bandwidth, 
C sampling points
C lambda: wavelength (in meters)
C Fractional bandwidth = 2 (Fh - Fl)/(Fh + Fl)
C where Fh is the highest frequency limit with signal 10dB below peak emission
C where Fl is the lowest frequency limit with signal 10dB below peak emission
          else if(comment(1:10) .eq. 'Wavelength')then
            read(f_in,*) lambda, fractional_bw, nsamp
            write(6,*) 'OK: lambda=',lambda,
     1           'm  fractional_bw=',fractional_bw,' nsamp=', nsamp,
     1           '(hence bandwidth =',fractional_bw * lambda,' m)'
C Keyword "Source", parameters: number of points, followed by
C Intensity1, X1, Y1 
C Intensity2, X2, Y2... 
          else if(comment(1:6) .eq. 'Source')then
            read(f_in,*) ncomp
            write(6,*) 'OK: Source with ',ncomp,' components'
            do n=1,ncomp
              read(f_in,*) s_int(n),s_pos(n,1),s_pos(n,2)
              write(6,*)' int, x, y: ',s_int(n),s_pos(n,1),s_pos(n,2)
            end do
C Keyword "Aperture", parameters: number of apertures, followed by
C Xcenter, Ycenter, radius, type (0 or 1)
          else if(comment(1:8) .eq. 'Aperture')then
            read(f_in,*) ncircle
            write(6,*) 'OK: Aperture with ',ncircle,' circles'
            do n=1,ncircle
              read(f_in,*) c_pos(n,1),c_pos(n,2),c_rad(n),c_typ(n)
              write(6,*)' x, y, rad, type: ',
     1                    c_pos(n,1),c_pos(n,2),c_rad(n),c_typ(n)
            end do
C Keyword "Phase screen", parameters: x_min, x_max, step 
          else if(comment(1:5) .eq. 'Phase')then
            read(f_in,*) x_min, x_max, xstep
            read(f_in,*) y_min, y_max, ystep
            write(6,*) 'OK: xmin, xmax, xstep:', x_min, x_max, xstep
            write(6,*) 'OK: ymin, ymax, ystep:', y_min, y_max, ystep
C Keyword "Frame", parameters: nber of frames 
          else if(comment(1:6) .eq. 'Frames')then
            read(f_in,*) nframes
            write(6,*) 'OK: nframes=',nframes
C Keyword "No disturbance" (no atmosphere)
          else if( comment(1:2) .eq. 'No') then
            no_dist = 1
            write(6,*) 'OK: no atmospheric disturbance'
          end if
C EOF loop 
   40   continue

C Switch to label=100 in case of end of file:
  100   close(unit=f_in)
        return
        end
*****************************************************************
* Recenter Fourier transform
* and shift origin from (0,0) to (Ns/2+1, Ns/2+1) or inversely
*
* Old recentering version from original program was not compatible
* with "fourn1" (and relative to an unknown data structure) 
*
* Quadrants 
*     3 4        2 1 
*     1 2    to  4 3 
*****************************************************************
       subroutine jlp_recent_cmplx(psc, Ns)
       implicit none
       integer Ns
       complex psc(Ns,Ns), tmp
       integer i, j
C Quadrants 
C     3 4        2 1 
C     1 2    to  4 3 
        do i=1,Ns/2
          do j=1,Ns/2
C Flip quadrants 2 and 3:
            tmp = psc(i+Ns/2,j)
            psc(i+Ns/2,j) = psc(i,j+Ns/2)
            psc(i,j+Ns/2) = tmp
C Flip quadrants 1 and 4:
            tmp = psc(i+Ns/2,j+Ns/2)
            psc(i+Ns/2,j+Ns/2) = psc(i,j)
            psc(i,j) = tmp
          end do
        end do

        return
        end
*****************************************************************
* Recenter Fourier transform
* and shift origin from (0,0) to (Ns/2+1, Ns/2+1) or inversely
*
* Old recentering version from original program was not compatible
* with "fourn1" (and relative to an unknown data structure) 
*
* Quadrants 
*     3 4        2 1 
*     1 2    to  4 3 
*****************************************************************
       subroutine jlp_recent_real(pwr, Np)
       implicit none
       integer Np
       real pwr(Np,Np), tmp
       integer i, j
C Quadrants 
C     3 4        2 1 
C     1 2    to  4 3 
        do i=1,Np/2
          do j=1,Np/2
C Flip quadrants 2 and 3:
            tmp = pwr(i+Np/2,j)
            pwr(i+Np/2,j) = pwr(i,j+Np/2)
            pwr(i,j+Np/2) = tmp
C Flip quadrants 1 and 4:
            tmp = pwr(i+Np/2,j+Np/2)
            pwr(i+Np/2,j+Np/2) = pwr(i,j)
            pwr(i,j) = tmp
          end do
        end do

        return
        end
***********************************************************************
* To create an image corresponding to the source
*
* INPUT:
* Np: size of object array (=Ns/2)
* Ns: size of phase screen (512)
* ncomp: number of components for the source
* s_int, s_pos: source information (intensity and position in arcsec)
* Ls: size of phase screen (in m)
* lambda: wavelength (in meters)
* Maxcomp: maximum number of components in source 
*
* OUTPUT:
* object(Np,Np): array with ones at the location of the components
* FT_object(Np,Np): Fourier Transform of the object
***********************************************************************
        subroutine create_object(object, FT_object, Np, Ns, ncomp, 
     1                          s_int, s_pos, Ls, lambda, Maxcomp)
        implicit none
        integer Np, Ns, ncomp, Maxcomp
        real*4 s_int(Maxcomp), s_pos(Maxcomp,2)
        real*4 object(Np,Np), Ls, lambda, Scale, pi, rapsec
        complex FT_object(Np,Np)
        integer i, j, ix, iy, Ncent, N(2)
        character out_file*60, out_comments*80

* Initialization of variables:
        pi = 3.141592653

* rapsec: radians per arcsecond
        rapsec = pi/(180. * 3600.0)
 
C Initializing the array:
        do i=1,Np
          do j=1,Np
            object(i,j) = 0.
          end do
        end do

* Scale: in pixels/arcsecond
        Scale = (Ls/lambda)*(real(Np)/real(Ns))*rapsec
        Ncent = Np/2 + 1
        write(6,*)' Scale:', Scale,' pixels/arcsec' 
 
C Loop on all components:
        do i=1,ncomp
* ncomp components in source structure
          ix = Ncent + int(Scale * s_pos(i,1))
          iy = Ncent + int(Scale * s_pos(i,2))
          ix = min(Np, max(1, ix))
          iy = min(Np, max(1, iy))
          write(6,*)' Component #',i,' ix=',ix,' iy=',iy
          object(ix,iy) = object(ix,iy) + s_int(i) 
        end do

C Transfer object to FT_object for Fourier transform
        do i=1,Np
          do j=1,Np
            FT_object(i,j)  = complex(object(i,j), 0.)
          end do
        end do

        N(1) = Np
        N(2) = Np
        call fourn1(FT_object,N,2,1)
 
C Output object for a diagnostic:
          out_file='s_object'
          out_comments='Source'
          call jlp_writeimag(object,Np,Np,Np,out_file,out_comments)
        return
        end
**************************************************************************
* Compute two elementary PSFs PSF1 and PSF2
* from a subwindow of the phase screen starting at x,y
* Take into account wavelength smearing due to finite bandpass
*
* INPUT:
* psc: phase screen to be used for computing the elementary PSFs 
* lambda: wavelength (in meters)
* fractional_bw: fractional bandwidth 
* nsamp: number of sampling points to take the finite bandwidth into account
* x, y : origin of the subwindow used to extract a phase screen from psc 
* ncirc: number of circles for the aperture
* c_pos, c_rad, c_typ: aperture information (coordinate, radius, type)
*  
* OUTPUT:
* PSF1, PSF2: elementary PSFs corresponding to real and imaginary
*                  parts of the phase screen
**************************************************************************
        subroutine compute_elementary_PSFs(psc, x, y, Ls,
     1                  lambda, fractional_bw, nsamp, 
     1                  ncircle, c_pos, c_rad, c_typ,
     1                  PSF1, PSF2, Ns, Np, Maxap, output_screen)
        implicit none
        real*4 x, y, Ls, lambda, fractional_bw
        integer Ns, Np, Maxap
        complex psc(Ns,Ns)
        real*4  PSF1(Np,Np), PSF2(Np,Np)
        integer nsamp, ncircle, output_screen
* Aperture information (coordinate, radius, type)
        real*4 c_pos(Maxap,2), c_rad(Maxap)
        integer*4 c_typ(Maxap)

* Temporary arrays
        complex PupilFunc1(Np,Np), PupilFunc2(Np,Np)

* Miscellaneous:
        real work(Np,Np), ww
        real factor, xa, ya, ox, oy, Lp
        integer loop1, jx, jy, ix, iy, lx, ly, i, Np_cent
        integer iz, jz, N(2)
        real Delta_lambda, lambda_min, l_step
        character out_file*60, out_comments*80
        complex iiu
* iiu = sqrt(-1) or "i"
        iiu = (0.0,1.0)

        Lp = Ls/Ns
        Np_cent = Np/2+1

C Reset image arrays
        do jx=1,Np
          do jy=1,Np
            PSF1(jx,jy) = 0.0
            PSF2(jx,jy) = 0.0
          end do
        end do
 
* Maximum variation of lambda is Delta_lambda
        Delta_lambda = fractional_bw * lambda
        lambda_min = lambda * ( 1. - fractional_bw / 2.)
        l_step = Delta_lambda / nsamp

******************* BEGIN OF LOOP1 ********************
********* on all the (nsamp) wavelengths belonging to the bandwidth *******
        do 100 loop1=1,nsamp
**************** intermediate loop *****************
* nsamp data points within the bandpass
* Old version: wrong!
*         factor =  1.0 + fractional_bw*(loop1 - INT(nsamp/2))
* JLP's version since psc is proportional to kk = 2 pi / lambda:
         factor =  lambda / (lambda_min + loop1 * l_step)
         write(6,*)'loop1=',loop1,' factor=',factor

C Reset temporary arrays
        do jx=1,Np
          do jy=1,Np
             PupilFunc1(jx,jy) = (0.0,0.0)
             PupilFunc2(jx,jy) = (0.0,0.0)
          end do
        end do
 
C Loop on all the circular apertures
          do 30 i=1,ncircle
C FOR DEBUG:
C            write(6,*)'Circle #',i,' centered on x+cposx, y+cposy'
C            write(6,*) '  x=',x,' y=',y,' cposx=',c_pos(i,1)
C            write(6,*)' c_posy=',c_pos(i,2),' c_rad(i)=',c_rad(i)
C 
*-----------------------------------------
*       filling pdata with aperture phases
 
           do 40 xa = x+c_pos(i,1)-c_rad(i), x+c_pos(i,1)+c_rad(i), Lp
           do 50 ya = y+c_pos(i,2)-c_rad(i), y+c_pos(i,2)+c_rad(i), Lp
             ox = xa-x-c_pos(i,1)
             oy = ya-y-c_pos(i,2)
             if( ox*ox+oy*oy .lt. c_rad(i)*c_rad(i) ) then
* Np_cent = Np/2+1
               ix = nint( (xa-x)/Lp )+Np_cent
               iy = nint( (ya-y)/Lp )+Np_cent
 
* Shift origin of the window to (Np_cent,Np_cent)
* Lp: pixel scale of phase screen (Lp = Ls/Ns)
               lx = nint(xa/Lp)
               ly = nint(ya/Lp)
               if(lx.lt.0.or.ly.lt.0) then
                write(6,*)'Fatal error: lx=',lx,' ly=',ly
                write(6,*)'Check if xmin > radius in parameter file...'
                stop
               endif
 
* iiu = sqrt(-1) or "i"
* NB: psc is proportional to kk = 2 pi / lambda
* correction by factor = lambda / (lambda + dlambda):
                PupilFunc1(ix,iy) = cexp(iiu*factor
     1                                *real(psc(lx,ly)) )*c_typ(i)
                PupilFunc2(ix,iy) = cexp(iiu*factor
     1                                *aimag(psc(lx,ly)) )*c_typ(i)
 
              end if
   50       continue
   40       continue
*-----------------------------------------
   30     continue

C Output real part of complex pupil function:
        if(output_screen.eq.0)then
          do iz=1, Np
            do jz=1, Np
              work(iz, jz) = real(PupilFunc1(iz,jz))
            end do
          end do
          out_file='s_PupilFunc1'
          out_comments='real part of exp (i phase screen)'
          call jlp_writeimag(work,Np,Np,Np,out_file,out_comments)
          output_screen = 1
        endif
 
* Simulate propagation from aperture plane to detector plane
 
C        call fourn1(PupilFunc1,N,2,-1,+1)
C        call fourn1(PupilFunc2,N,2,-1,+1)
C JLP2008 (by comparing with s_psc.for):
C fourn1(data,nn,ndim,isign)
        N(1) = Np
        N(2) = Np
        call fourn1(PupilFunc1,N,2,1)
        call fourn1(PupilFunc2,N,2,1)
 
* The intensity of the Fourier Transform of the Pupil function
* corresponds to an elementary image
*
* Recenter images (necessary for the subsequent convolution ...) 
        call jlp_recent_cmplx(PupilFunc1, Np)
        call jlp_recent_cmplx(PupilFunc2, Np)
 
* Form speckle/fringe patterns
* i.e., intensity of the complex field (PupilFunc1 or PupilFunc2)
        ww = 1. / real(Np * Np)
        do ix=1, Np
          do iy=1, Np
            PSF1(ix,iy) = PSF1(ix,iy) +
     1         ww*real(PupilFunc1(ix,iy)*conjg(PupilFunc1(ix,iy)))
            PSF2(ix,iy) = PSF2(ix,iy) +
     2         ww*real(PupilFunc2(ix,iy)*conjg(PupilFunc2(ix,iy)))
          end do
        end do
 
  100   continue
******************* END OF LOOP1 (on nsamp) ********************

C Output PSF1:
        if(output_screen.eq.1)then
          out_file='s_psf'
          out_comments='Intensity of complex field in image plane'
          call jlp_writeimag(PSF1,Np,Np,Np,out_file,out_comments)
          output_screen = 2
        endif

        return
        end
********************************************************************
* Compute elementary images by convolving PSF1 with source
* and add the corresponding power spectrum to pwr
*
* INPUT:
* PSF1: elementary PSF 
* FT_object: Fourier Transform of the source
* Np: size of power spectrum (=Ns/2)
*
* INPUT/OUTPUT:
* pwr(Np,Np): power spectrum
********************************************************************
       subroutine process_elementary_image(pwr, PSF1, 
     1                      FT_object, Np, output_screen)
       implicit none
       integer Np, output_screen
       real*4  pwr(Np,Np), PSF1(Np,Np)
       complex FT_object(Np,Np)

* Temporary arrays
       real*4 ary1(Np,Np)
       complex cary1(Np,Np), TransferFunc(Np,Np)

* Miscellaneous:
       real*4 work(Np,Np)
       integer N(2), jx, jy, iz, jz
       character out_file*60, out_comments*80

        N(1) = Np
        N(2) = Np

C JLP2008
        do jx=1,Np
          do jy=1,Np
             TransferFunc(jx,jy) = complex(PSF1(jx,jy),0.0)
          end do
        end do
        call fourn1(TransferFunc,N,2,1)

C Convolution is done in Fourier domain:
        do jx=1,Np
          do jy=1,Np
C Product is a complex product here!
             cary1(jx,jy) = TransferFunc(jx,jy) * FT_object(jx,jy) 
          end do
        end do

        do jx=1,Np
        do jy=1,Np
          pwr(jx,jy) = pwr(jx,jy) 
     1                 + real(cary1(jx,jy)*conjg(cary1(jx,jy)))
        end do
        end do
 
C Output of an elementary frame:
        if(output_screen.eq.2)then

        call fourn1(cary1,N,2,-1)
        do jx=1,Np
          do jy=1,Np
            ary1(jx,jy) = real(cary1(jx,jy)) 
          end do
        end do
        call jlp_recent_real(ary1, Np)

          do iz=1, Np
            do jz=1, Np
              work(iz, jz) = ary1(iz,jz)
            end do
          end do
          out_file='s_elem_frame'
          out_comments='PSF1 convolved by source'
          call jlp_writeimag(work,Np,Np,Np,out_file,out_comments)
          output_screen = 3
        endif
 
        return
        end
