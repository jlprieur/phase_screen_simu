        program s_simu
*******************************************************************************
*
* s_simu.for (from s_pscorg.for) is a phase screen simulator which produces 
* artificial fringe
* or speckle patterns and obtain integrated power spectrum or(/and)
* optionally specified bispectral components.
*
* The atmospheric disturbance is given by the Kolmogorov spectrum.
* Seeing condition, aperture configuration and the object ( at this
* point, only multiple point sources ) and the usage of the phase screen
* must be specified by the user in an ascii input file.
*
* Example: 
*     runs s_simu simu_triple.input triple phot_gv7_128.fits
*
* Following is a sample input file. (free format, ! indicates a
* comment line.)
*
* ! Size of images
* Size
* 128
* ! mean number of photons per frame
* Photons
* 200
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
* ! central wavelength(nm), bandwidth (nm), sampling points
* Wavelength
* 650 70 3
* ! source: # of point sources
* !         intensity(arbitrary unit), position offset(x,y) in arcsecond
* Source
* 2
* 1.0  -0.7 -0.6
* 0.5   0.8 0.4
* ! # of circular apertures (including voids)
* ! configuraton: center of aperture(x,y), diameter(m) and type
* ! 1 is an aperture,0 is a void. Voids should follow ordinary apertures.
* ! Example: annular mask
* Aperture
* 2
* 0.0 0.0 2.0  1
* 0.0 0.0 1.9  0
* ! total number of frames to be processed
* Frames
* 100
* ------ as an option `no atmosphere' can be specified by ------
* ! obtaining beam pattern
* No_disturbance
* --------------------------------------------------------------
*
* From s_pscorg.for (but nearly nothing left from original version...)
* I still leave the old comments...
*
* This is a re-organized version of the pscsn.for developed on Tacos.
* The program is planned to be as users oriented as possible for the
* future development of Optical/IR interferometry.
*
* version 1.0    Dec 25, 1988 (Merry Christmas!)     Tadashi Nakajima
* Today's resolution : `I will have real holidays next year and live
*                       a human life.'
* version 2.0    Dec 30, 1988   Now sharp cut-offs of fluctuation power
* spectrum at kapmin and kapmax are introduced to investigate the
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
*   I replace the fractional bandwidth = 2 (Fh - Fl)/(Fh + Fl)
*   fractional bandwidth = 2 (Fh - Fl)/(Fh + Fl)
*   where Fh is the highest frequency limit with signal 10dB below peak emission
*   where Fl is the lowest frequency limit with signal 10dB below peak emission
*   with the bandwidth \approx (Fh - Fl) in the parameter file
*   - I change the screen generator (jlp_phase_screen) so that the 
*   simulation is compatible with the r_zero value (the old version was not)
*   - I replace the pseudo-convolution of the objet with the PSF by a true 
*   convolution
*   - I implement a good formula for the image smearing due to finite bandwidth
*   (the original was wrong!) 
*   - I also put the good formula in jlp_phase_screen to obtain simulated
*   frames that are consistent with the value of r0!
*
* JLP
* Version 21/07/2008
*******************************************************************************
*------ declarations ---------*
        implicit none
* Constants
        real pi
 
* Integrated structure constant of the phase fluctuation
* LCn2 = \int Cn2(h) dh
* And also: 
* r_0 = (0.423 kk^2 LCn2)^{-3/5}
        real LCn2, r_zero
 
* Size of phase screen and resolution (m/pix)
        real Ls, Scale_s
 
* Size of phase screen, size of power spectrum
        integer Ns, Np, idim2, idim, NGMAX
* For IR=30: NGMAX=388400 
        parameter(idim2 = 256, NGMAX = 388400)
C        parameter(idim2 = 512, NGMAX = 388400)
        parameter(idim = idim2/2)
 
* Wave number corresponding to lambda_cent, scale length of turbulence
* fov: field of view (for the simulated images)
        real kk_cent, L0, lambda_cent, lo, fov
 
* psc: complex phase screen, real and imag are independent.
* FT_object(idim,idim): Fourier Transform of the object
        complex psc(idim2,idim2), FT_object(idim,idim)
* Detector transfer function:
        complex DetectorTransf(idim,idim)
 
* pwr: mean power spectrum of analog elementary frames 
* modsq: mean power spectrum with photon clipping and detector photon response
* long_int: long integration with photon clipping and detector photon response 
        real*8 pwr(idim,idim), modsq(idim,idim), long_int(idim,idim)
        real*8 snrm(idim,idim), bisp1(4*NGMAX)
 
* Input control parameter file and output power spectrum file
        character input_file*60, out_prefix*40, photon_modsq*60
 
* Maxcomp: maximum number of components in source 
* Maxap: maximum number of apertures
        integer Maxcomp, Maxap
        parameter(Maxcomp = 50, Maxap = 15)

* nframes: number of frames
* frm_per_screen: number of frames per phase screen in X or in Y
* nph: mean number of photons per frame
        integer  nframes, frm_per_screen, nph
 
* Source information (intensity and position in arcsec)
        real s_int(Maxcomp), s_pos(Maxcomp,2), object(idim,idim)
 
* Aperture information (coordinate, radius, type)
        real c_pos(Maxap,2), c_rad(Maxap)
        integer c_typ(Maxap)
* Miscellaneous:
        real x_min, x_max, xstep, y_min, y_max, ystep
        real rapsec, bandwidth, image_scale, fwhm_pix, r0_pix 
        integer ncircle, nsamp, i, iframe
        integer ncomp, seed, no_dist, output_demo
        integer ir, max_nclosure, nbeta, ngamma
 
* Initialization of variables:
        pi = 3.141592653

* rapsec: radians per arcsecond
        rapsec = pi/(180. * 3600.0)
 
* JLP
* Program can then be run as "runs s_simu s_screen.input simu1"
        call jlp_begin
        call jlp_inquifmt
 
* Input file names (from terminal or "jlp_lu5.tmp" file when run with "runs")
        write(6,20)
   20   format(1x,'Input ascii file name: ',$)
        read(5,'(a)') input_file
        write(6,21)
   21   format(1x,'Prefix for output file names: ',$)
        read(5,'(a)') out_prefix
        write(6,22)
   22   format(1x,'Power spectrum of detector photon response',
     1         ' (0 if none):',$)
        read(5,'(a)') photon_modsq 
        write(6,25)
   25   format(1x)

        output_demo = 0
 
C First read the parameter file (in order to obtain the value of Np):
        call read_param_file(input_file, Np, seed, nph, L0, lo, 
     1                lambda_cent, bandwidth, nsamp, 
     1                ncomp, s_int, s_pos, 
     1                ncircle, c_pos, c_rad, c_typ, fov, r_zero, 
     1                nframes, no_dist, Maxcomp, Maxap, out_prefix)

C Too many problems with fourn1 if idim!=Np, so...
        if(Np.ne.idim)then
          write(6,*)'Fatal error: Np=',Np,' not equal to idim=',idim
          write(6,*)'Recompile the program with idim=',Np
          stop
        endif
        Ns = 2 * Np

C Then read the photon response of the detector and check that nx=Np:
        call read_photon_response(photon_modsq, Np, 
     1               DetectorTransf, idim)

C Computing the uv coverage: 
        write(6,*) ' Radius of uv-coverage (IR) in pixels:',ir
        ir = 30 
        max_nclosure = 200
        call covera(ir, max_nclosure, nbeta, ngamma)
        write(6,*) ' nbeta=',nbeta,' ngamma=',ngamma
        if(ngamma.gt.NGMAX)then
          write(6,*)' Fatal error, ngamma maxi =',NGMAX
          stop
        endif

C Initialization of random generator:
        call jlp_random_init(seed)

*  kk_cent: wave number corresponding to lambda_cent 
        kk_cent = 2.0 * pi / lambda_cent
 
* Original version with LCn2 = 5.0e-13
* r_0 = (0.423 kk^2 LCn2)^{-3/5}
        LCn2 = (r_zero**(3./5.)) / (0.423 * kk_cent * kk_cent) 
        write(6,287) lambda_cent, r_zero, LCn2
287     format('lambda_cent=',e12.5,' r_zero=',f8.3,' LCn2=',e12.5)

* rapsec = pi/(180. * 3600.0)
* Scale of image in pixels/arcsecond: Scale = (Ls/lambda_cent)*(Np/Ns)*rapsec
        image_scale = Np/fov
* Field of view = Np/Scale
* Size of phase screen and resolution (m/pix)
        Ls = Ns / ((fov / lambda_cent) * rapsec)
        write(6,*) 'Size of screen: Ls = ',Ls,' m'
        Scale_s = Ls / real(Ns)
        write(6,*) 'Resolution of screen: Scale_s = ',Scale_s,' m/pixel'
        r0_pix = r_zero/Scale_s 
        write(6,*) 'r0 corresponds to: ',r0_pix,' pixels'
        write(6,*) 'Image scale:  ',image_scale,' pix/"'
        write(6,*) ' Check:',
     1             (Ls/lambda_cent)*(real(Np)/real(Ns))*rapsec,' pix/"'
* rapsec = pi/(180. * 3600.0)
* and image_scale in pixel/arcsec
        fwhm_pix = ((lambda_cent / r_zero) / rapsec) * image_scale
        write(6,*) 'FWHM:  ',fwhm_pix,' pixels'

C Create an image of the source from ncomp, s_int, s_pos
         call create_object(object, FT_object, Np, Ns, ncomp, 
     1             s_int, s_pos, Ls, lambda_cent, Maxcomp, 
     1             out_prefix, idim)

* frm_per_screen: number of frames per phase screen in X or in Y
        frm_per_screen = 4
C Compute sampling parameters corresponding to frm_per_screen: 
        call compute_sampling(frm_per_screen, Ls, Ns, 
     1                      x_min, x_max, xstep, y_min, y_max, ystep,
     1                      ncircle, c_pos, c_rad, Maxap)
*
* Initialize arrays:
        call init_arrays(pwr, modsq, snrm, long_int, bisp1, 
     1                   Np, ngamma, idim)

* Two real phase screens per complex random number generation in Fourier
* domain.
 
        write(6,*)' ******** Will now simulate ', nframes,
     1            ' elementary frames ************'
        write(6,*)' Number of frames per phase screen: ',
     1              frm_per_screen*frm_per_screen 
        write(6,*)' Number of calls of phase_screen routine needed: ', 
     1            max(1,(nframes/(frm_per_screen * frm_per_screen))/2)
 
 
        iframe = 0
* Main loop (maximum number of iterations will be nframes)
        do 300 i = 1, nframes 
*------------------------------------
* Create a complex phase screen
        write(6,310) i, iframe+1
  310   format('+', 'generating ',i6,
     1         'th complex phase screen (iframe=',i6,')')
C New version:
        call jlp_phase_screen(psc, Ns, idim2, seed, Ls, r_zero, 
     1                     r0_pix, no_dist, out_prefix, output_demo) 
C Old version:
C        call  phase_screen(psc, Ns, idim2, seed, Ls, L0, lo, 
C     1                     kk_cent, LCn2, out_prefix)
 
* Sample phase screen, form images, calculate power spectrum
        call  process_frames(psc, Ns, FT_object, DetectorTransf, 
     1              pwr, modsq, long_int, snrm, bisp1, Np, iframe, 
     1              nframes, Ls, seed, nph, ir, nbeta, ngamma,
     1              x_min, x_max, xstep, y_min, y_max, ystep,
     1              lambda_cent, bandwidth, nsamp, ncomp,
     1              s_int, s_pos, ncircle, c_pos, c_rad, c_typ,
     1              Maxcomp, Maxap, out_prefix,
     1              output_demo, idim, idim2)
 
          if(iframe .ge. nframes) goto 400
*
*------------------------------------
* End of loop
  300   continue
*

* Exit from loop
  400   continue


C Output FITS files:
        call output_results(pwr, modsq, long_int, snrm, 
     1                      bisp1, Np, ngamma, nframes, 
     1                      input_file, out_prefix, idim)

*
        call jlp_end
        stop
        end
***************************************************************************
* Creates complex phase screen (JLP's version)
*
* INPUT:
* Ns: size of phase screen 
* seed: seed for the random number generator, 
* output_demo = 0 changed to 1 in exit
* r_zero: Fried's parameter (in meters)
* Scale_s: Screen scale (in meter/pixel)
* Ls: size of screen in meters (corresponding to Ns pixels)
* no_dist: flag set to one if no atmospheric disturbance
* Np: size of image
* fov: field of view
*
* OUPUT:
*
***************************************************************************
        subroutine jlp_phase_screen(psc, Ns, idim2, seed, Ls, r_zero, 
     1                      r0_pix, no_dist, out_prefix, output_demo)
        implicit none
        integer Ns, idim2, seed, output_demo, no_dist
        integer fourn_format
        complex psc(idim2,idim2),iiu
        real work(idim2,idim2), structure_funct(idim2,idim2)
        real Ls, r_zero, pi, rad2, mean, sigma, r0_pix 
        real*8 sum, sumsq
        logical debug
        character out_file*60, out_comments*80
        character out_prefix*40
        real urad, radi, Fs, Fs0, utheta 
 
* N(): dimension array for fft
* Ns_cent: origin of the phase screen
        integer N(2), Ns_cent, jx, jy, i, j
 
* No atmospheric disturbance
        if (no_dist.eq.1) then
          do jx=1, Ns
            do jy=1, Ns
              psc(jx,jy) = (0.0,0.0)
            end do
          end do
          return
        endif

        debug=.true.
* Constants:
* iiu = sqrt(-1) or "i"
        iiu = (0.0,1.0)
        pi = 3.141592653

* r_0 = (0.423 kk^2 LCn2)^{-3/5}
* Hence Fs0 is proportional to kk^2
C JLP2008 / old version (June 2008):
C        Fs0 = 0.023 * Ns * Ns * Ns * Ns * (L0 / r_zero)**(5./3.)
C JLP2008: new version (July 2008)
* by this value:
*        Fs0 = 0.023 * Ns * Ns / r_zero**(5./3.)
* Formula derived by JLP:
        Fs0 = sqrt(0.023 * Ns * Ns * (Ls/ r_zero)**(5./3.))
* Formula from Porro (2000, Applied Optics 39, 10, p 1644)
*        Fs0 = sqrt(0.023) * (Ls / r_zero)**(5./6.)
* JLP's ad'hoc correction:
*        Fs0 = Fs0 * Ns

*------ generation of random phase in spatial domain -------------
* Complex Gaussian random numbers are generated by
* the inverse transform method.
 
        N(1) = Ns
        N(2) = Ns
        Ns_cent = Ns/2+1
 
        do 100 jy=1, Ns
          do 110 jx=1, Ns
            rad2=(jx-Ns_cent)*(jx-Ns_cent)+(jy-Ns_cent)*(Jy-Ns_cent)
 
* rad**(-11/3) is undefined at rad2=0.
            if(rad2 .eq. 0.0) rad2 = 0.5 
 
10            call jlp_random(urad) 
              if ( urad.gt.0.0 ) then
                radi = sqrt(-2.0*log(urad))
              else
                goto 10
              end if
C 11/6=1.833333
              Fs = rad2**(-11.0/6.0)
              call jlp_random(utheta)

* iiu = sqrt(-1) or "i"
* WARNING: next statement stores a complex number to psc 
* since cexp is a complex exponential!
              psc(jx,jy) = Fs0 * sqrt(Fs)*radi*cexp(iiu*2.0*pi*utheta)
C For DEBUG (calibration)
C              psc(jx,jy) = radi*cexp(iiu*2.0*pi*utheta)
  110     continue
  100   continue
 
*---------- Fourier transform ------------------------------
* Shift origin of Fourier space to (0,0) before inverse Fourier transform
        call jlp_recent_cmplx(psc, Ns, idim2)

C JLP2008 (by comparing with s_psc.for):
C fourn1(data,nn,ndim,isign)
C SHOULD USE MY fourn1 routine (in fft_jlp.for)
C NOT s_fourn1.for !!!
        call fourn1(psc,N,2,-1)
 
C Correction since division by Ns with inverse FFT:
C (Calibrated in July 2008 with random Gauss, and obtained sigma=1 in Fourier
C domain from sigma=1 in direct domain)
        do j=1, Ns
          do i=1, Ns
           psc(i,j) = psc(i,j) * Ns
          end do
        end do

* Now psc is a complex random phase screen
* Center array to avoid a discontinuity of phase screen at (Ns/2,Ns/2)
        call jlp_recent_cmplx(psc, Ns, idim2)

        write(6,*)'OK6: psc(10,20) =',real(psc(10,20)),
     1              imag(psc(10,20))
C Output phase screen:
        if(output_demo.eq.0)then
          sum = 0.
          sumsq = 0.
          do j=1, Ns
            do i=1, Ns
              work(i, j) = real(psc(i,j))
              sum = sum + work(i,j)
              sumsq = sumsq + work(i,j) * work(i,j)
            end do
          end do
          mean = sum / real(Ns * Ns)
          sigma = sqrt(sumsq / real(Ns * Ns) - mean * mean)
          write(6,*)'Phase screen: mean=', mean,' sigma=',sigma
          write(out_file, 87) out_prefix(1:index(out_prefix,' ')-1)
87        format(A,'_phase_screen')
          out_comments='Real part of phase screen'
          call jlp_writeimag(work,Ns,Ns,idim2,out_file,out_comments)

C Structure function of the real part of the complex phase array:
          fourn_format = 1
          call fstructure_real_with_ft(psc,Ns,Ns,idim2,
     1          structure_funct,fourn_format)
          write(out_file, 88) out_prefix(1:index(out_prefix,' ')-1)
88        format(A,'_phase_struct')
          out_comments='Struct. func. of real part of cplx phase screen'
          write(6,*)' Output of structure function to: ',out_file
          call jlp_writeimag(structure_funct,Ns,Ns,idim2,
     1                       out_file,out_comments)
          output_demo = 1
        endif

        return
        end
***************************************************************************
* Creates complex phase screen (original version: problems of
* compatibility of LCn2 and the r_zero obtained in long integrations...))
*
* INPUT:
* Ns: size of phase screen
* seed: seed for the random number generator, 
* kk_cent: wave number corresponding to lambda_cent
* LCn2 integrated CN2 over all the altitudes (originally: LCn2=5.0e-13)
*
* OUPUT:
*
***************************************************************************
        subroutine phase_screen(psc, Ns, idim2, seed, Ls, L0, lo, 
     1                          kk_cent, LCn2, out_prefix)
        implicit none
        integer Ns, idim2, seed
        complex psc(idim2,idim2),iiu
        real dkap, kapmin, kapmax, kk_cent, LCn2, pi, Ls, L0, lo, Scale_s
        real work(idim2,idim2), mean, sigma, ww
        real*8 sum,sumsq
        logical debug
        character out_file*60, out_comments*80
        character out_prefix*40
 
* kap: radial spatial frequency in aperture plane ( * 2 pi / L)
        real kap, atten, urad, radi, Fs, Fs0, utheta 
 
* N(): dimension array for fft
* Ns_cent: origin of the phase screen
        integer N(2), Ns_cent, jx, jy, i, j
 
        debug=.true.
* Constants:
* iiu = sqrt(-1) or "i"
        iiu = (0.0,1.0)
        pi = 3.141592653

* Scale_s: pixel scale of phase screen
        Scale_s = Ls/Ns

* kapmin, kapmax: min and max spatial frequencies in aperture plane
* Ls: size of phase screen (in m)
* Scale_s: screen resolution, size of one pixel (in m/pix)
        if (Ls .lt. L0) then
          kapmin = 2.0 * pi / Ls
        else
          kapmin = 2.0 * pi / L0
        end if
 
        if (lo .lt. Scale_s) then
           kapmax = 2.0 * pi / Scale_s
        else
           kapmax = 2.0 * pi / lo
        end if
 
* dkap: one pixel in spatial frequency in aperture plane (dkap = 2 pi / Ls)
        dkap = 2.0 * pi / Ls

* dkap: one pixel in spatial frequency in aperture plane (dkap = 2 pi / Ls)
* r_0 = (0.423 kk^2 LCn2)^{-3/5}
        Fs0 = 0.033*pi*LCn2*kk_cent*kk_cent*2.0*dkap*dkap

*------ generation of random phase in spatial domain -------------
* Complex Gaussian random numbers are generated by
* the inverse transform method.
 
        N(1) = Ns
        N(2) = Ns
        Ns_cent = Ns/2+1
 
        do 100 jy=1, Ns
          do 110 jx=1, Ns
 
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
* Fs = 0.033*pi*LCn2*kk_cent*kk_cent*2.0*kap**(-11.0/3.0)*dkap*dkap
            Fs = Fs0 * kap**(-11.0/3.0)
            call jlp_random(utheta)
* iiu = sqrt(-1) or "i"
            psc(jx,jy) = atten*sqrt(Fs)*radi*cexp(iiu*2.0*pi*utheta)
  110     continue
  100   continue
 
*---------- Fourier transform ------------------------------
* Shift origin of Fourier space to (0,0)
* Old recentering version from original program was wrong!
        call jlp_recent_cmplx(psc, Ns, idim2)

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
        call jlp_recent_cmplx(psc, Ns, idim2)

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
          write(out_file, 87) out_prefix(1:index(out_prefix,' ')-1)
87        format(A,'_phase_screen')
          out_comments='Real part of phase screen'
          call jlp_writeimag(work,Ns,Ns,idim2,out_file,out_comments)
        endif

        return
        end
*******************************************************************************
* process_frames
* 
* Produces fringe/speckle patterns of multiple point sources.
* Finite bandwidth effect is taken into account.
*
* INPUT:
* psc(idim2,idim2): phase screen
* Ns: size of phase screen
* FT_object: Fourier Transform of the source
* DetectorTransf: detector transfer function (i.e., modulus of Fourier 
*                 Tranform of PSF)
* Np: size of elementary frames (=Ns/2)
* iframe: index of current frame to be processed
* nframes: number of frames to be simulated
* ncomp: number of components for the source
* nph: mean number of photons per frame
* s_int, s_pos: source information (intensity and position in arcsec)
* ncirc: number of circles for the aperture
* c_pos, c_rad, c_typ: aperture information (coordinate, radius, type)
* Ls: size of phase screen (in m)
* lambda_cent, bandwidth, nsamp: wavelength parameters: central wavelength (m),
*                          bandwidth (m) and sampling points
* x_min, x_max, xstep : sampling parameters to extract a phase screen 
*                       from a larger one 
* y_min, y_max, ystep : sampling parameters to extract a phase screen 
*                       from a larger one 
* Maxcomp: maximum number of components in source 
* Maxap: maximum number of apertures
* output_demo: 1 at first call, and it is incremented in exit
*
* OUTPUT:
* pwr: mean power spectrum of analog elementary frames 
* modsq: mean power spectrum with photon clipping and detector photon response
* long_int: long integration with photon clipping and detector photon response 
*******************************************************************************
        subroutine process_frames(psc, Ns, FT_object, DetectorTransf,
     1                  pwr, modsq, long_int, snrm, bisp1, Np, iframe, 
     1                  nframes, Ls, seed, nph, ir, nbeta, ngamma, 
     1                  x_min, x_max, xstep, y_min, y_max, ystep,
     1                  lambda_cent, bandwidth, nsamp, ncomp,
     1                  s_int, s_pos, ncircle, c_pos, c_rad, c_typ,
     1                  Maxcomp, Maxap, out_prefix,
     1                  output_demo, idim, idim2)
        implicit none
        integer Ns, Np, Maxcomp, Maxap, seed, nph
        integer idim, idim2, ir, nbeta, ngamma
        complex psc(idim2,idim2), FT_object(idim,idim)
        complex DetectorTransf(idim,idim)
        real*8 pwr(idim,idim), modsq(idim,idim), long_int(idim,idim)
        real*8 snrm(idim,idim), bisp1(4*ngamma)
        real Ls, lambda_cent, bandwidth 
        real x_min, x_max, xstep, y_min, y_max, ystep
        integer iframe, nframes, nsamp, ncomp, ncircle, output_demo
        character out_prefix*40

* Aperture information (coordinate, radius, type)
        real c_pos(Maxap,2), c_rad(Maxap)
        integer c_typ(Maxap)

* Source information (intensity and position in arcsec)
        real s_int(Maxcomp), s_pos(Maxcomp,2)
 
* Temporary arrays
        real  PSF1(idim,idim), PSF2(idim,idim)

C Miscellaneous:
        integer N(2), Np_cent
        real  Scale_s, Scale, x, y
        logical debug

* Constants
        real pi, rapdeg, rapsec
        complex iiu
* iiu = sqrt(-1) or "i"
        iiu = (0.0,1.0)
        pi = 3.141592653

* rapsec: radians per arcsecond
        rapdeg = pi/180.0
        rapsec = rapdeg/3600.0
        debug = .false.
 
        N(1) = Np
        N(2) = Np
* Ls: size of phase screen (in m)
* Scale_s: pixel scale of phase screen
        Scale_s = Ls/Ns
        Np_cent = Np/2+1
* lambda_cent: wavelength (in meters)
* Scale: in pixels/arcsecond
        Scale = (Ls/lambda_cent)*(real(Np)/real(Ns))*rapsec

* Setting apertures on Np sized window in Ns sized complex phase screen
 
        do 110 x=x_min,x_max,xstep
        do 120 y=y_min,y_max,ystep
**************************** major loop ***************************************
 
        if(debug)then
        write(6,111) x,y
  111   format('+',' creating a pair of speckle/fringe patterns at',
     1         'x=',f10.3,' y=',f10.3)
        endif
 
C DEBUG:
C        write(6,333) ncircle, nsamp, c_typ(1), c_typ(2)
C  333   format(1x,'ncircle nsamp typ(1) typ(2) ',1x,4(i2,1x) )
 
C Compute two elementary PSFs (PSF1 and PSF2)
        call compute_elementary_PSFs(psc, x, y, Ls, lambda_cent, 
     1                  bandwidth, nsamp, ncircle, c_pos, c_rad, c_typ,
     1                  PSF1, PSF2, Ns, Np, Maxap, output_demo, 
     1                  out_prefix, idim, idim2)

C Convolution of the object with the elementary PSFs 
        
C Processing a first frame with real part of phase screen:
        if((iframe.lt.10).or.(mod(iframe,20).eq.0))then
          write(6,320) iframe+1
        endif
  320   format('+', 'processing ',i5,' th frame                    ')
       call process_elementary_frame(pwr, modsq, long_int, snrm,
     1                        bisp1, PSF1, FT_object, DetectorTransf,
     1                        Np, seed, nph, ir, nbeta, ngamma,
     1                        output_demo, out_prefix, idim)
        iframe = iframe + 1
        if( iframe .ge. nframes) return

C Processing a second frame with imag part of phase screen:
       call process_elementary_frame(pwr, modsq, long_int, snrm,
     1                        bisp1, PSF2, FT_object, DetectorTransf,
     1                        Np, seed, nph, ir, nbeta, ngamma,
     1                        output_demo, out_prefix, idim)
        iframe = iframe + 1
        if( iframe .ge. nframes) return
 
************************* end of major loop ***********************************
  120   continue
  110   continue

        return
        end
***************************************************************************
* read_param_file: read input parameter file
*
* INPUT:
* input_file: name of input parameter file 
*
* OUTPUT:
* Np: size of elementary frames (=Ns/2)
* seed: seed for random number generator (large odd integer)
* L0, lo: external and internal scale lengths of turbulence
* lambda_cent, bandwidth, nsamp: wavelength parameters: central wavelength (m),
*                          bandwidth (m) and sampling points
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
        subroutine read_param_file(input_file, Np, seed, nph, L0, lo, 
     1                  lambda_cent, bandwidth, nsamp, ncomp,
     1                  s_int, s_pos, ncircle, c_pos, c_rad, c_typ,
     1                  fov, r_zero, nframes, no_dist, Maxcomp, Maxap,
     1                  out_prefix)
        implicit none
        integer Np, Maxcomp, Maxap
        character out_prefix*40
        character input_file*60, comment*80
        real L0, lo, lambda_cent, bandwidth 
        real s_int(Maxcomp), s_pos(Maxcomp,2)
        real c_pos(Maxap,2), c_rad(Maxap), fov, r_zero
        integer seed, nph, nsamp, ncomp, ncircle, nframes, no_dist
        integer c_typ(Maxap), j, n
* I/O units
        integer f_in
        parameter(f_in = 11)
 
* Initialize no_dist to default value:
        no_dist = 0
        nph = 200
        Np = 128

* Read the input file
        open(unit=f_in, file=input_file, status='old')
 
   30   format(a80)
        do 40 j=1,100
 
          read(f_in,30,end=100) comment
          if( comment(1:1) .eq. '!') goto 40

          write(6,30) comment
 
C Keyword "Size of images", parameter: Np 
          if(comment(1:4) .eq. 'Size')then
            read(f_in,*) Np 
            write(6,*) 'OK: Np=',Np
C Keyword "Photons per frame"
          else if(comment(1:7) .eq. 'Photons')then
            read(f_in,*)  nph 
            write(6,*) 'OK: nph=',nph,' photons/frame'
C Keyword "Field of view", parameter: fov 
          else if(comment(1:5) .eq. 'Field')then
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
C Keyword "Fried radius", parameter: r_zero 
          else if(comment(1:5) .eq. 'Fried')then
            read(f_in,*) r_zero
            write(6,*) 'OK: r0=',r_zero
C Keyword "Wavelength", parameters: 
* lambda_cent, bandwidth, nsamp: central wavelength (nm), bandwidth (nm) 
* and sampling points
          else if(comment(1:10) .eq. 'Wavelength')then
            read(f_in,*) lambda_cent, bandwidth, nsamp
C Conversion from nm to m:
            lambda_cent = lambda_cent * 1.e-9
            bandwidth = bandwidth * 1.e-9
            write(6,*) 'OK: lambda_cent=',lambda_cent,
     1           'm  bandwidth=',bandwidth,'m  nsamp=', nsamp
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
C Xcenter, Ycenter, diameter, type (0 or 1)
          else if(comment(1:8) .eq. 'Aperture')then
            read(f_in,*) ncircle
            write(6,*) 'OK: Aperture with ',ncircle,' circles'
            do n=1,ncircle
              read(f_in,*) c_pos(n,1),c_pos(n,2),c_rad(n),c_typ(n)
C Conversion from diameter to radius:
              c_rad(n) = c_rad(n)/2.
              write(6,*)' x, y, rad, type: ',
     1                    c_pos(n,1),c_pos(n,2),c_rad(n),c_typ(n)
            end do
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
       subroutine jlp_recent_cmplx(psc, Ns, idim2)
       implicit none
       integer Ns, idim2
       complex psc(idim2,idim2), tmp
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
       subroutine jlp_recent_real(ft_array, Np, idim)
       implicit none
       integer Np, idim
       real ft_array(idim,idim), tmp
       integer i, j
C Quadrants 
C     3 4        2 1 
C     1 2    to  4 3 
        do i=1,Np/2
          do j=1,Np/2
C Flip quadrants 2 and 3:
            tmp = ft_array(i+Np/2,j)
            ft_array(i+Np/2,j) = ft_array(i,j+Np/2)
            ft_array(i,j+Np/2) = tmp
C Flip quadrants 1 and 4:
            tmp = ft_array(i+Np/2,j+Np/2)
            ft_array(i+Np/2,j+Np/2) = ft_array(i,j)
            ft_array(i,j) = tmp
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
       subroutine jlp_recent_real8(ft_array, Np, idim)
       implicit none
       integer Np, idim
       real*8 ft_array(idim,idim), tmp
       integer i, j
C Quadrants 
C     3 4        2 1 
C     1 2    to  4 3 
        do i=1,Np/2
          do j=1,Np/2
C Flip quadrants 2 and 3:
            tmp = ft_array(i+Np/2,j)
            ft_array(i+Np/2,j) = ft_array(i,j+Np/2)
            ft_array(i,j+Np/2) = tmp
C Flip quadrants 1 and 4:
            tmp = ft_array(i+Np/2,j+Np/2)
            ft_array(i+Np/2,j+Np/2) = ft_array(i,j)
            ft_array(i,j) = tmp
          end do
        end do

        return
        end
***********************************************************************
* To create an image corresponding to the source
*
* INPUT:
* Np: size of elementary frames (=Ns/2)
* Ns: size of phase screen 
* ncomp: number of components for the source
* s_int, s_pos: source information (intensity and position in arcsec)
* Ls: size of phase screen (in m)
* lambda_cent: wavelength (in meters)
* Maxcomp: maximum number of components in source 
*
* OUTPUT:
* object(idim,idim): array with ones at the location of the components
* FT_object(idim,idim): Fourier Transform of the object
***********************************************************************
        subroutine create_object(object, FT_object, Np, Ns, ncomp, 
     1                          s_int, s_pos, Ls, lambda_cent, 
     1                          Maxcomp, out_prefix, idim)
        implicit none
        integer Np, Ns, ncomp, Maxcomp, idim
        real s_int(Maxcomp), s_pos(Maxcomp,2)
        real object(idim,idim), Ls, lambda_cent, Scale, pi, rapsec
        complex FT_object(idim,idim)
        integer i, j, ix, iy, Ncent, N(2)
        character out_prefix*40
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
        Scale = (Ls/lambda_cent)*(real(Np)/real(Ns))*rapsec
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
        write(out_file, 87) out_prefix(1:index(out_prefix,' ')-1)
87      format(A,'_object')
        out_comments='Source'
        call jlp_writeimag(object,Np,Np,idim,out_file,out_comments)

        return
        end
**************************************************************************
* Compute two elementary PSFs PSF1 and PSF2
* from a subwindow of the phase screen starting at x,y
* Take into account wavelength smearing due to finite bandpass
*
* INPUT:
* psc: phase screen to be used for computing the elementary PSFs 
* lambda_cent: central wavelength (in meters)
* bandwidth: half width of filter used (in meters) 
* nsamp: number of sampling points to take the finite bandwidth into account
* x, y : offset for the subwindow used to extract a phase screen from psc 
* ncirc: number of circles for the aperture
* c_pos, c_rad, c_typ: aperture information (coordinate, radius, type)
*  
* OUTPUT:
* PSF1, PSF2: elementary PSFs corresponding to real and imaginary
*                  parts of the phase screen
**************************************************************************
        subroutine compute_elementary_PSFs(psc, x, y, Ls,
     1                  lambda_cent, bandwidth, nsamp, 
     1                  ncircle, c_pos, c_rad, c_typ,
     1                  PSF1, PSF2, Ns, Np, Maxap, output_demo,
     1                  out_prefix, idim, idim2)
        implicit none
        real x, y, Ls, lambda_cent, bandwidth 
        integer Ns, Np, Maxap, idim, idim2
        complex psc(idim2,idim2)
        real  PSF1(idim,idim), PSF2(idim,idim)
        integer nsamp, ncircle, output_demo
* Aperture information (coordinate, radius, type)
        real c_pos(Maxap,2), c_rad(Maxap)
        integer c_typ(Maxap)
        character out_prefix*40

* Temporary arrays
        complex PupilFunc1(idim,idim), PupilFunc2(idim,idim)

* Miscellaneous:
        real work(idim,idim), ww
        real factor, xa, ya, ox, oy, Scale_s
        integer isamp, jx, jy, ix, iy, lx, ly, icirc, Np_cent
        integer iz, jz, N(2)
        real lambda_min, l_step, xcent, ycent
        character out_file*60, out_comments*80
        complex iiu
* iiu = sqrt(-1) or "i"
        iiu = (0.0,1.0)

        Scale_s = Ls/Ns
        Np_cent = Np/2+1

C Reset image arrays
        do jx=1,Np
          do jy=1,Np
            PSF1(jx,jy) = 0.0
            PSF2(jx,jy) = 0.0
          end do
        end do
 
* Maximum variation of lambda is bandwidth 
        lambda_min = lambda_cent - bandwidth / 2.
        l_step = bandwidth / real(nsamp)

******************* BEGIN OF LOOP1 ON ISAMP ********************
********* on all the (nsamp) wavelengths belonging to the bandwidth *******
        do 100 isamp=1,nsamp
**************** intermediate loop *****************
* nsamp data points within the bandpass
* Old version: wrong!
*         factor =  1.0 + fractional_bw*(isamp - INT(nsamp/2))
* JLP's version since psc is proportional to kk_cent = 2 pi / lambda_cent:
         factor =  lambda_cent / (lambda_min + isamp * l_step)

C DEBUG
C         write(6,*)'isamp=',isamp,' factor=',factor

C Reset temporary arrays
        do jx=1,Np
          do jy=1,Np
             PupilFunc1(jx,jy) = (0.0,0.0)
             PupilFunc2(jx,jy) = (0.0,0.0)
          end do
        end do
 
C Loop on all the circular apertures
          do 30 icirc=1,ncircle
              xcent = x+c_pos(icirc,1)
              ycent = y+c_pos(icirc,2)
C FOR DEBUG:
            write(6,*)'Circle #',icirc,' centered on x+cposx, y+cposy'
            write(6,*) '  x=',x,' y=',y
            write(6,*) '  xcent=',xcent,' ycent=',ycent
C 
*-----------------------------------------
*       filling pdata with aperture phases
 
           do 40 xa = xcent-c_rad(icirc), xcent+c_rad(icirc), Scale_s
           do 50 ya = ycent-c_rad(icirc), ycent+c_rad(icirc), Scale_s
             ox = xa - xcent
             oy = ya - ycent
             if( ox*ox+oy*oy .lt. c_rad(icirc)*c_rad(icirc) ) then
* Np_cent = Np/2+1
               ix = nint( (xa - x)/Scale_s )+Np_cent
               iy = nint( (ya - y)/Scale_s )+Np_cent
               if(ix.lt.0.or.iy.lt.0.or.ix.gt.Np.or.iy.gt.Np) then
                write(6,*)'Fatal error: ix=',ix,' iy=',iy
                write(6,*)'xa-x/Scale_s=',(xa - x )/ Scale_s
                write(6,*)'x=',x,' y=',y,' xcent=',xcent,' ycent=',ycent
                write(6,*)'xa=',xa,' ya=',ya,' cpos1',c_pos(icirc,1),
     1               ' cpos2',c_pos(icirc,2),' crad=',c_rad(icirc)
                write(6,*)' (NB: Np=',Np,', Scale_s=',Scale_s,
     1                    ' Np_cent=',Np_cent,')'
                write(6,*)"Check pupil parameters in parameter file..."
                stop
               endif
 
* Shift origin of the window to (Np_cent,Np_cent)
* Scale_s: pixel scale of phase screen (Scale_s = Ls/Ns)
               lx = nint(xa/Scale_s)
               ly = nint(ya/Scale_s)
               if(lx.lt.0.or.ly.lt.0.or.lx.gt.Ns.or.ly.gt.Ns) then
                write(6,*)'Fatal error: lx=',lx,' ly=',ly
                write(6,*)' (NB: Ns=',Ns,')'
                write(6,*)'Check if xmin > radius in parameter file...'
                stop
               endif
 
* iiu = sqrt(-1) or "i"
* NB: psc is proportional to kk = 2 pi / lambda
* correction by factor = lambda_cent / (lambda_cent + dlambda):
                if(ix.eq.50.and.iy.eq.80)then
                  write(6,*)' JLPPP: lx, ly=',lx,ly,' ix, iy=',ix,iy
                  write(6,*)' factor=',factor,' psc:',psc(lx,ly)
                endif
                PupilFunc1(ix,iy) = cexp(iiu*factor
     1                                *real(psc(lx,ly)) )*c_typ(icirc)
                PupilFunc2(ix,iy) = cexp(iiu*factor
     1                                *aimag(psc(lx,ly)) )*c_typ(icirc)
 
              end if
   50       continue
   40       continue
*-----------------------------------------
   30     continue

C Output real part of complex pupil function:
        if(output_demo.eq.1)then
          do iz=1, Np
            do jz=1, Np
              work(iz, jz) = real(PupilFunc1(iz,jz))
            end do
          end do
          write(out_file, 87) out_prefix(1:index(out_prefix,' ')-1)
87        format(A,'_PupilFunc')
          out_comments='real part of exp (i phase screen)'
          call jlp_writeimag(work,Np,Np,idim,out_file,out_comments)
          output_demo = 2
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
        call jlp_recent_cmplx(PupilFunc1, Np, idim)
        call jlp_recent_cmplx(PupilFunc2, Np, idim)
 
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
******************* END OF LOOP1 ON ISAMP ********************

C Output PSF1:
        if(output_demo.eq.2)then
          write(out_file, 88) out_prefix(1:index(out_prefix,' ')-1)
88        format(A,'_psf')
          out_comments='Intensity of complex field in image plane'
          call jlp_writeimag(PSF1,Np,Np,idim,out_file,out_comments)
        endif

        return
        end
********************************************************************
* Compute an elementary frame by convolving PSF1 with source
* and pile up the corresponding power spectrum and long integration 
*
* INPUT:
* PSF1: elementary PSF 
* FT_object: Fourier Transform of the source
* DetectorTransf: detector transfer function (i.e., modulus of Fourier 
*                 Tranform of PSF)
* nph: number of photons per frame
* Np: size of elementary frames (=Ns/2)
*
* INPUT/OUTPUT:
* pwr: mean power spectrum of analog elementary frames 
* modsq: mean power spectrum with photon clipping and detector photon response
* long_int: long integration with photon clipping and detector photon response 
* output_demo: value of a dummy variable used to know whether
*                an image has to be output
********************************************************************
       subroutine process_elementary_frame(pwr, modsq, long_int, snrm, 
     1                        bisp1, PSF1, FT_object, DetectorTransf, 
     1                        Np, seed, nph, ir, nbeta, ngamma,
     1                        output_demo, out_prefix, idim)
       implicit none
       integer Np, output_demo, seed, nph, idim
       integer ir, nbeta, ngamma
       real*8 pwr(idim,idim), modsq(idim,idim), long_int(idim,idim)
       real*8 snrm(idim,idim), bisp1(4*ngamma)
       real PSF1(idim,idim)
       complex FT_object(idim,idim), DetectorTransf(idim,idim)
       character out_prefix*40

* Temporary arrays
       real ary1(idim,idim), ary1_ph(idim,idim)
       real*8 re(idim,idim),im(idim,idim)
       complex cary1(idim,idim), cary1_ph(idim,idim)
       complex TransferFunc(idim,idim)

* Miscellaneous:
       integer N(2), jx, jy, nphr
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

C Pile up the analog power spectrum:
        do jx=1,Np
        do jy=1,Np
          pwr(jx,jy) = pwr(jx,jy) 
     1                 + real(cary1(jx,jy)*conjg(cary1(jx,jy)))
        end do
        end do
 
C Inverse Fourier Transform to create an analog image:
        call fourn1(cary1,N,2,-1)
        do jx=1,Np
          do jy=1,Np
            ary1(jx,jy) = real(cary1(jx,jy)) 
          end do
        end do
        call jlp_recent_real(ary1, Np, idim)

C Output one of those elementary frames:
        if(output_demo.eq.2)then
          write(out_file, 87) out_prefix(1:index(out_prefix,' ')-1)
87        format(A,'_elem_frame_an')
          out_comments='PSF1 convolved by source'
          call jlp_writeimag(ary1,Np,Np,idim,out_file,out_comments)
        endif

C Photon clipping:
        call photon_clipping(ary1, ary1_ph, Np, nph, nphr, seed, idim)
 
C Output one of those photon clipped images:
       if(output_demo.eq.2)then
         write(out_file, 88) out_prefix(1:index(out_prefix,' ')-1)
88       format(A,'_elem_frame_ph')
         out_comments='PSF1 convolved by source and photon clipped'
         call jlp_writeimag(ary1_ph,Np,Np,idim,out_file,out_comments)
       endif

C Fourier transform to compute mosq:
        do jx=1,Np
          do jy=1,Np
             cary1_ph(jx,jy) = cmplx(ary1_ph(jx,jy),0.) 
          end do
        end do

        call fourn1(cary1_ph,N,2,1)

C Convolution with detector response (i.e product in Fourier domain):
        do jx=1,Np
          do jy=1,Np
C Product is a complex product here!
             cary1(jx,jy) = cary1_ph(jx,jy) * DetectorTransf(jx,jy) 
          end do
        end do

C Loading Fourier transform to re,im arrays:
        do jx=1,Np
          do jy=1,Np
             re(jx,jy) = real(cary1(jx,jy))
             im(jx,jy) = imag(cary1(jx,jy))
          end do
        end do


* Processing the spectrum of this image now:
* int bispec3_single(float *re, float *im, float *modsq, float *snrm,
*            INT4 *nx, INT4 *ny, float *bispp, INT4 *ir,
*            INT4 *nbeta, INT4 *ngamma)
        call BISPEC3(re,im,modsq,snrm,Np,Np,bisp1,ir,
     1               nbeta,ngamma)
 
C Back to direct domain:
        call fourn1(cary1,N,2,-1)

C Pile up the long integration 
        do jx=1,Np
        do jy=1,Np
          ary1(jx,jy) = real(cary1(jx,jy)) 
          long_int(jx,jy) = long_int(jx,jy) + ary1(jx,jy) 
        end do
        end do

C Output one of those elementary frames (photon clipped 
C and convolved with detector response) 
       if(output_demo.eq.2)then
         write(out_file, 89) out_prefix(1:index(out_prefix,' ')-1)
89       format(A,'_elem_frame')
         out_comments='Simulation of output from detector'
         call jlp_writeimag(ary1,Np,Np,idim,out_file,out_comments)
         output_demo = 3
       endif


        return
        end
**************************************************************************
* Before calling poidev, the sum of in_array has to be normalised to nph
*
* INPUT:
* in_array: analog elementary frame according to Poisson statistics
* Np: size of elementary frames (=Ns/2)
* nph: number of photons/frame wanted
* seed: seed used for Poisson random generator
*
* OUTPUT:
* out_array: elementary frame that was digitized according to Poisson statistics
* nphr: number of photons of the output frame
**************************************************************************
       subroutine photon_clipping(in_array, out_array, Np, nph, nphr, 
     1                            seed, idim)
       implicit none
       integer Np, nph, nphr, seed, idim
       real in_array(idim,idim), out_array(idim,idim)
       real*8 sum
       real w_scale
       integer iz, jz
       logical debug
       debug = .false.

* First normalise the sum of image (ary1) to nph photons (but still analogic) 
       sum = 0.
        do iz=1, Np
          do jz=1, Np
            sum = sum + in_array(iz,jz)
          end do
        end do
    
       w_scale = real(nph) / sum
        do iz=1, Np
          do jz=1, Np
            in_array(iz,jz) = in_array(iz,jz) * w_scale
          end do
        end do

* Before calling poidev, the sum of prf has to be normalised to nphr
       nphr = 0

* Poisson generator
        do iz=1, Np
          do jz=1, Np
            call jlp_poisson(in_array(iz,jz), out_array(iz,jz), seed)
* Computing total number of photons: 
            nphr = nphr + out_array(iz,jz)
          end do
        end do

       if(debug)then
         write(6,*) ' # of photons after clipping: ', nphr 
       endif

       return
       end
**************************************************************************
* INPUT:
* frm_per_screen: number of frames per phase screen in X or in Y
* Ls: width of phase screen (in m)
* Ns: width of phase screen (in pixels)
* ncirc: number of circles for the aperture
* c_pos, c_rad, c_typ: aperture information (coordinate, radius, type)
* Maxap: maximum number of apertures
*
* OUTPUT:
* x_min, x_max, xstep : sampling parameters to extract a phase screen 
*                       from a larger one 
* y_min, y_max, ystep : sampling parameters to extract a phase screen 
*                       from a larger one 
**************************************************************************
       subroutine compute_sampling(frm_per_screen, Ls, Ns, 
     1                   x_min, x_max, xstep, y_min, y_max, ystep, 
     1                   ncircle, c_pos, c_rad, Maxap)
       implicit none
       integer  frm_per_screen, Ns, ncircle, Maxap
       real x_min, x_max, xstep, y_min, y_max, ystep
       real c_pos(Maxap,2), c_rad(Maxap)
       real Ls, aperture_diameter
       integer i

       aperture_diameter = 0.
       do i=1, ncircle
        aperture_diameter = max(aperture_diameter, 2. * c_rad(i))
       end do
C Np is Ns/2 hence screen size is Ls/2 
       if(aperture_diameter.gt.Ls/2.)then
         write(6,*)'Fatal error, aperture diameter=',aperture_diameter
         write(6,*)'larger than size of phase screens: ',Ls/2.
         stop
       endif

C x_min,ymin: lower left coordinates of the subwindow to be used
C to extract a small phase screen from the large phase screen
C First compute x_min, and y_min:
       x_min = Ls 
       y_min = Ls
       do i=1, ncircle
        x_min = min(x_min, c_pos(i, 1) - c_rad(i))
        y_min = min(y_min, c_pos(i, 2) - c_rad(i))
       end do
C Tolerance is 0.1 m to avoid zero at the edge...
       if(x_min.le.0.)then
          x_min = 0.1 - x_min
       else
          x_min = 0.1
       endif
       if(y_min.le.0.)then
          y_min = 0.1 - y_min
       else
          y_min = 0.1
       endif

C Then x_max, y_max:
       x_max = 0.
       y_max = 0.
       do i=1, ncircle
        x_max = max(x_max, c_pos(i, 1) + c_rad(i))
        y_max = max(y_max, c_pos(i, 2) + c_rad(i))
       end do
       if(x_max.ge.0.)then
         x_max = Ls - 0.1 - x_max
       else
         x_max = Ls - 0.1
       endif
       if(y_max.ge.0.)then
         y_max = Ls - 0.1 - y_max
       else
         y_max = Ls - 0.1
       endif

C Then xstep, ystep:
       if(frm_per_screen.gt.1)then
         xstep = (x_max - x_min) / (frm_per_screen - 1)
         ystep = (y_max - y_min) / (frm_per_screen - 1)
       else
         x_max = x_min
         y_max = y_min
         xstep = 1000. 
         ystep = 1000. 
       endif

C Error analysis:
       if(xstep.le.0..or.ystep.le.0.
     1    .or.x_min.ge.Ls.or.y_min.ge.Ls
     1    .or.x_max.le.0..or.y_max.le.0.)then
         write(6,*)'Fatal error computing sampling parameters'
         write(6,*)'x_min, x_max, y_min, y_max:',x_min, x_max, 
     1              y_min, y_max
         write(6,*)'xstep, ystep:',xstep, ystep
         stop
       endif

         write(6,*)'x_min, x_max, y_min, y_max:',x_min, x_max, 
     1              y_min, y_max
         write(6,*)'xstep, ystep:',xstep, ystep
       return
       end
***********************************************************************
* To read the FITS file containing the photon response of the detector 
*
* INPUT:
*  photon_modsq: name of the FITS file containing the detector response
*
* OUTPUT:
*  DetectorTransf: detector transfer function 
***********************************************************************
        subroutine read_photon_response(photon_modsq, Np, 
     1               DetectorTransf, idim)
        implicit none
        integer Np, idim
        complex DetectorTransf(idim,idim)
        character photon_modsq*60
* Miscellaneous:
        real detector_modsq(idim,idim), w_scale
        character comments*80
        integer iz, jz, nx, ny

        if(photon_modsq(1:1).eq.'0'.or.photon_modsq(1:1).eq.' ')then
          do iz=1,Np
          do jz=1,Np
             DetectorTransf(iz,jz) = cmplx(1.,0.)
          end do
          end do
        else
          call jlp_readimag(detector_modsq,nx,ny,idim,
     1                      photon_modsq,comments)
          if(nx.ne.Np.or.ny.ne.Np)then
             write(6,*)'Fatal error inconsistent size: nx,ny=',nx,ny
             write(6,*)'whereas Np=',Np
             stop
          endif

C Shift center to (1,1)
          call jlp_recent_real(detector_modsq, Np, idim)

C Normalization to unity:
          w_scale = 1./sqrt(detector_modsq(1,1))
C Transfer function corresponds to sqrt(modsq): 
          do iz=1,Np
          do jz=1,Np
             DetectorTransf(iz,jz) = 
     1               cmplx(w_scale * sqrt(detector_modsq(iz,jz)),0.)
          end do
          end do
        endif
        return
        end
********************************************************************
********************************************************************
        subroutine output_results(pwr, modsq, long_int, snrm, 
     1                            bisp1, Np, ngamma, nframes, 
     1                            input_file, out_prefix, idim)
        implicit none
        integer Np, ngamma, nframes, idim
        real*8 pwr(idim,idim), modsq(idim,idim), long_int(idim,idim)
        real*8 snrm(idim,idim), bisp1(4*ngamma)
        character out_prefix*40, input_file*60

* Miscellaneous:
        real*8 w_scale, w1, w2, w3, bisp2(ngamma,3)
        integer lx, ly, i1, i2, i3, i4
        character out_file*80, comments*80

C Recenter Fourier transforms:
        call jlp_recent_real8(pwr, Np, idim)
        call jlp_recent_real8(modsq, Np, idim)
        call jlp_recent_real8(snrm, Np, idim)

* Compute the mean power spectrum and long integration: 
* To normalize power spectrum to unity, set:  w_scale = modsq(Np/2+1,Np/2+1)
* Otherwise:
        w_scale = 1./real(nframes)
        do lx=1, Np
          do ly=1, Np
            pwr(lx,ly) = pwr(lx,ly) * w_scale
            modsq(lx,ly) = modsq(lx,ly) * w_scale
            snrm(lx,ly) = snrm(lx,ly) * w_scale
            long_int(lx,ly) = long_int(lx,ly) * w_scale
          end do
        end do
* Computing the mean:
        do lx=1, 4 * ngamma 
            bisp1(lx) = bisp1(lx) * w_scale
         end do

* The variance:
        do lx=1, ngamma 
* sum(re)
            i1 = 4*lx - 3
* sum(im)
            i2 = 4*lx - 2
* sumsq(re)
            i3 = 4*lx - 1
* sumsq(im)
            i4 = 4*lx 
* variance of re:
            w1 = bisp1(i3) - bisp1(i1)*bisp1(i1)
* variance of im:
            w2 = bisp1(i4) - bisp1(i2)*bisp1(i2)
* Full variance (variance of real + variance of imag):
            w1 = w1 + w2
            w1 = max(w1,1.e-10)
* Standard deviation:
            w1 = sqrt(w1)
* Computing the phase factor on the first two lines:
            w3 = bisp1(i1)*bisp1(i1) + bisp1(i2)*bisp1(i2)
            w3 = max(w3,1.e-10)
            w3 = sqrt(w3)
C            bisp2(lx,1) = bisp1(i1)/w3
C            bisp2(lx,2) = bisp1(i2)/w3
C As in vcrb: I save the bispectrum instead of the phase factor:
C to be able to display the bispectrun with "bisp_to_image.c"
C and to be compatible with other programs:
            bisp2(lx,1) = bisp1(i1)
            bisp2(lx,2) = bisp1(i2)
* SNR of bispectrum on the third line:
            bisp2(lx,3) = w3 / w1
        end do
 
* Output the mean long integration to a FITS file:
        write(out_file, 81) out_prefix(1:index(out_prefix,' ')-1)
81      format(A,'_l')
        write(comments,82) input_file,char(0)
82      format('Mean long integration with: ',A,A)
        call jlp_d_writeimag(long_int,Np,Np,idim,out_file,comments)

* Output the analog power spectrum to a FITS file:
        write(out_file, 83) out_prefix(1:index(out_prefix,' ')-1)
83      format(A,'_power')
        write(comments,84) input_file,char(0)
84      format('Analog power spectrum with: ',A,A)
        call jlp_d_writeimag(pwr,Np,Np,idim,out_file,comments)

* Output the SNR of the power spectrum to a FITS file:
        write(out_file, 85) out_prefix(1:index(out_prefix,' ')-1)
85      format(A,'_snrm')
        write(comments,86) input_file,char(0)
86      format('SNR of modsq with: ',A,A)
        call jlp_d_writeimag(snrm,Np,Np,idim,out_file,comments)

* Output the fully simulated power spectrum to a FITS file:
        write(out_file, 87) out_prefix(1:index(out_prefix,' ')-1)
87      format(A,'_m')
        write(comments,88) input_file,char(0)
88      format('Simulated power spectrum with: ',A,A)
        call jlp_d_writeimag(modsq,Np,Np,idim,out_file,comments)

* Output the bispectrum to a FITS file:
        write(out_file, 89) out_prefix(1:index(out_prefix,' ')-1)
89      format(A,'_b')
        write(comments,90) input_file,char(0)
90      format('Bispectrum with: ',A,A)
        call jlp_d_writeimag(bisp2,ngamma,3,ngamma,out_file,comments)

        return
        end
*************************************************************************
        subroutine init_arrays(pwr, modsq, snrm, long_int, bisp1, 
     1                         Np, ngamma, idim)
        implicit none
        integer iz, jz, Np, ngamma, idim
        real*8 pwr(idim,idim), modsq(idim,idim), snrm(idim,idim)
        real*8 long_int(idim,idim), bisp1(4*ngamma)

        do iz=1,Np
         do jz=1,Np
           pwr(jz,jz) = 0. 
           modsq(jz,jz) = 0. 
           snrm(jz,jz) = 0. 
           long_int(jz,jz) = 0. 
         end do
        end do
        do iz=1,4*ngamma
          bisp1(iz) = 0.
        end do
        return
        end
