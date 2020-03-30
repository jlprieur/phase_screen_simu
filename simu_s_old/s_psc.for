C Head Simulation of disturbed phase screen
        program s_psc
*******************************************************************************
*
*       is a phase screen simulator which produces artifitial fringe
*       or speckle patterns and obtain integrated power spectrum or(/and)
*       optionally specified bispectral components.
*       The atmospheric disturbance is given by the Kolmogorov spectrum.
*       Seeing condition, aperture configuration and the object ( at this
*       point, only multiple point sources ) and the usage of the phase screen
*       must be specified by the user.
*
*       Written in standard FORTRAN and expected to run both on VMS and UNIX
*       simply by replacing FFT routines and random number generator.
*       On UNIX all tempory arrays are to be dynamically allocated by 'calloc()'
*       Source code should be called pscreen.for on VMS, pscreen.f on UNIX.
*       On a vector machine, effeciency of loops should be reconsidered.
*
*CMV No need to 'calloc()' arrays on Alliant
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
*
*       Changed to work with SCASIS/GHRIL interface, options to write
*       phase screen and detector plane.  Marco de Vos, 10-12-90
*
*       Optimised many nested loops by changing their order. Removed bug
*       in simulation of double source. Marco de Vos, may/jun 91
*
*******************************************************************************
*------ declarations ---------*
* Constants
        real      pi, rapdeg, rapmin, rapsec
        complex iu
 
* Integrated structure constant of the phase fluctuation
        real   LCn2
 
* Size of phase screen (m) and resolution (m/pix)
        real   Ls, Lp
 
* Size of phase screen (pix), size of power spectrum (pix)
        integer   Ns, Np
        parameter( Ns = 512  )
        parameter( Np = Ns/2 )
 
* Wave number, spatial frequencies and scale length of tubulence
        real k, kapmax, kapmin, L0, dkap, lo
 
* Wavelength, fractional bandwidth, # of samples
        real lambda, fraction
        integer nsamp
 
* Complex phase screen, real and imag are independent.
        complex psc(Ns,Ns)
 
* Power spectrum array
        real   pwr(Np,Np)
 
* Maximum # of components in source, maximum number of apertures
        integer   Maxcomp, Maxap
        parameter( Maxcomp = 50, Maxap = 15 )
 
* Source information (intensity and position in arcsec)
        real   s_int(Maxcomp), s_pos(Maxcomp,2)
 
* Aperture information (coordinate, radius, type)
        real   c_pos(Maxap,2), c_rad(Maxap)
        integer   c_typ(Maxap)
 
* Definition of subframes
        real   x_min, x_max, xstep, y_min, y_max, ystep
 
* Number of sources, apertures, frames
        integer ncomp, ncircle, nframe
*
        common/psc/ pi,rapdeg,rapmin,rapsec, iu,
     &              k, kapmax, kapmin, L0, LCn2, Ls, Lp,
     &              dkap, lambda, nframe,
     &              x_min, x_max, xstep, y_min, y_max, ystep,
     &              ncomp, ncircle, nsamp, fraction
 
* Work area for subroutine power
        real    image1(Np,Np), image2(Np,Np),  ary1(Np,Np), ary2(Np,Np)
        complex pdata1(Np,Np), pdata2(Np,Np), cary1(Np,Np), cary2(Np,Np)
 
* Variables local to main program
        real      zero
        integer   seed, n, iscreen, iframe, iz, jz, lx, ly
        logical   no_dist,w_screen,w_image,w_power
        character comment*80
 
C Head Initialisation, user input
 
* Common variables are not compatible with parameter statement.
        pi = 3.141592653
        rapdeg = pi/180.0
        rapmin = rapdeg/60.0
        rapsec = rapdeg/3600.0
        iu = (0.0,1.0)
        LCn2 = 5.0e-13
        Ls = 10.24
        Lp = 0.02
 
* Open log-file etc, get input parameters
        call init()
 
        nframe=50
        x_min=2.5
        y_min=2.5
        x_max=8.5
        y_max=8.5
        xstep=.15
        ystep=8.0
        call userint(nframe,1,1,'NFRAME=','Enter the number of frames')
        call userrpair(x_min,y_min,1,'SUB_MIN=',
     &        'Enter minimum X,Y of subframe (m)')
        call userrpair(x_max,y_max,1,'SUB_MAX=',
     &        'Enter maximum X,Y of subframe (m)')
        call userrpair(xstep,ystep,1,'SUB_STEP=',
     &        'Enter stepsize in X,Y (m)')
 
        L0=0.60
        lo=0.02
        call userrpair(L0,lo,2,'SCALE=',
     &                'Enter outer,inner scale of turbulence (m)')
 
        lambda=0.63e-6
        fraction=0.00476
        nsamp=3
        call userreal(lambda,1,2,'LAMBDA=','Enter wavelength (m)')
        call userreal(fraction,1,2,'FBAND=','Enter fractional bandwidth')
        call userint(nsamp,1,2,'NSAMP=','Enter number of sampling points')
 
        ncircle=1
        c_rad(1)=1.
        c_pos(1,1)=0.0
        c_pos(1,2)=0.0
        c_typ(1)=1
        call userint(ncircle,1,1,'NCIRCLE=','Enter number of apertures')
        do n=1,ncircle
           call userrpair(c_pos(n,1),c_pos(n,2),1,'CPOS=',
     &                   'Enter position in X,Y (m)')
           call userreal(c_rad(n),1,1,'CRAD=','Enter radius (m)')
           call userint(c_typ(n),1,1,'CTYPE=','Enter type (1=hole,0=void)')
           c_rad(n+1)=c_rad(n)
           c_pos(n+1,1)=c_pos(n,1)
           c_pos(n+1,2)=c_pos(n,1)
           c_typ(n+1)=c_typ(n)
           call cancel('CPOS=')
           call cancel('CRAD=')
           call cancel('CTYPE=')
        end do
 
        ncomp=1
        call userint(ncomp,1,1,'NCOMP=','Enter number of sources')
        do n=1,ncomp
           s_int(n)=1.
           s_pos(n,1)=0.0
           s_pos(n,2)=0.0
           call userrpair(s_pos(n,1),s_pos(n,2),1,'SPOS=',
     &                   'Enter position in X,Y (m)')
           call userreal(s_int(n),1,1,'SINT=','Enter intensity')
           call cancel('SPOS=')
           call cancel('SINT=')
        end do
 
        no_dist=.FALSE.
        w_screen=.FALSE.
        w_image=.FALSE.
        w_power=.TRUE.
        call userlog(no_dist,1,2,'NDIST=','No disturbance?')
        call userlog(w_screen,1,2,'WSCREEN=','Write screen data?')
        call userlog(w_image, 1,2,'WIMAGE=', 'Write image plane?')
        call userlog(w_power, 1,2,'WPOWER=', 'Write power spectrum?')
 
        seed=7597551
        call userint(seed,1,2,'SEED=','Enter seed for random generator')
 
* k: wave number, kap: spatial frequency in aperture plane
 
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
 
* dkap: one pixel in spatial frequency in aperture plane
        dkap = 2.0 * pi / Ls
 
C Head Loop over frames
 
*
        write(1,*) 'frm_per_scr ',
     1     (int((x_max-x_min)/xstep)+1)*(int((y_max-y_min)/ystep)+1)
C JLP 91        call easyou(0)
 
        iframe=1
        iscreen=1
*------------------------------------
287     if ( no_dist ) then
*            no atmospheric disturbance
             do iz=1, Ns
               do jz=1, Ns
                 psc(iz,jz) = (0.0,0.0)
               end do
             end do
           else
* Create a complex phase screen
             write(comment,310) iscreen
  310        format(' Generating ',i3,'th complex phase screen')
             call  screen(psc, Ns, seed)
             iscreen=iscreen+1
             if (w_screen)
     &           call psc_wscreen(psc,Ns,Ns,pdata1,pdata2,Np,Np)
           end if
 
* Sample phase screen, form images, calculate power
           call  power(psc, Ns, pwr, Np, iframe, Maxcomp, Maxap,
     1            s_int, s_pos, c_pos, c_rad, c_typ, pdata1, pdata2,
     2            image1, image2, ary1, ary2, cary1, cary2)
           if (w_image) call psc_wimage(ary1,ary2,Np,Np)
*
        if(iframe.le.nframe)goto 287
*-----------------------------------------
*
 
* Write normalized integrated power spectrum
 
        zero = pwr(Np/2+1,Np/2+1)
        if (zero.eq.0.0) zero=1.0
        write(1,*) 'Scaling: divide by ',zero
C JLP 91        call easyou(0)
* CMV REV
        do ly=1, Np
          do lx=1, Np
            pwr(lx,ly) = pwr(lx,ly)/zero
          end do
        end do
        print *,'Scaled'
        if (w_power) call psc_wpower(pwr,Np,Np)
*
        n=0
        call data_close(1,n)
        n=0
        call data_close(2,n)
        n=0
        call data_close(3,n)
        call finis()
 
        end
 
C Head Compute phase screen
 
        subroutine screen(psc, Ns, seed)
***************************************************************************
* Creates complex phase screen
***************************************************************************
        implicit  none
        integer Ns,seed
        complex psc(Ns,Ns)
 
        real    pi,rapdeg,rapmin,rapsec,Ls,Lp
        complex iu
        real    k, kapmax, kapmin, L0, LCn2, dkap, lambda, fraction
        real    x_min,x_max,xstep, y_min,y_max,ystep
        integer ncomp,ncircle,nsamp,nframe
*
        common/psc/ pi,rapdeg,rapmin,rapsec, iu,
     &              k, kapmax, kapmin, L0, LCn2, Ls, Lp,
     &              dkap, lambda, nframe,
     &              x_min, x_max, xstep, y_min, y_max, ystep,
     &              ncomp, ncircle, nsamp, fraction
 
* Radial spatial frequency in aperture plane
        real   kap
 
* Dimension array for fft, origin of the phase screen
        integer   N(2), Nc
 
        integer   i,j,jx,jy
        real      Fs,Fs0,atten,radi,u,ran
        character comment*80
 
*------ generation of random phase in spatial domain -------------
* Complex Gaussian random numbers are generated by
* the inverse transform method.
 
        N(1) = Ns
        N(2) = Ns
        call four2(psc,Ns,Ns,0)
        Nc = Ns/2+1
        Fs0 = 0.033*pi*LCn2*k*k*2.0*dkap*dkap
* CMV REV
        do 100 jy=1, Ns
 
          if (mod(jy,100).eq.0) then
             write(comment,*) 'Making line ',jy
             call status(comment)
          end if
 
          do 110 jx=1, Ns
 
            kap = dkap*sqrt( float( (jx-Nc)**2+(jy-Nc)**2 ) )
 
* Introducing smooth cut-offs below kapmin and beyond kapmax.
 
            if(kap .gt. kapmax) atten = exp(-(kap-kapmax)/(5.0*dkap))
            if(kap .lt. kapmin) atten = exp(-(kapmin-kap)/(5.0*dkap))
            if(kap.le.kapmax .and. kap.ge.kapmin) atten = 1.0
 
* kap**(-11/3) is undefined at kap=0.
            if(kap .eq. 0.0) kap = dkap
 
            u=0.0
288           u=ran(seed)
              if(u.le.0.0) goto 288

            radi = sqrt(-2.0*log(u))
            Fs = Fs0*kap**(-11.0/3.0)
            psc(jx,jy) = atten*sqrt(Fs)*radi*cexp(iu*2.0*pi*ran(seed))
 
  110     continue
  100   continue
 
C Page
 
*---------- Fourier transform ------------------------------
* Shift origin of Fourier space to (Nc,Nc)
* CMV REV
        do j=1,Ns,2
          do i=1,Ns,2
            psc(i+1,j) = -psc(i+1,j)
            psc(i,j+1) = -psc(i,j+1)
          end do
        end do
 
        write(comment,555) N(1),N(2)
  555   format(' ', 'start ',i3,'x',i3,' sized fft              ')
        call four2(psc,Ns,Ns,-1)
        write(comment,666)
  666   format(' ', 'end fft for screen                         ')
 
* Now psc is a complex random phase screen
* to avoid discontinuity shift origin of phase screen to (Nc,Nc)
* CMV REV
        do j=1,Ns,2
          do i=1,Ns,2
            psc(i+1,j) = -psc(i+1,j)
            psc(i,j+1) = -psc(i,j+1)
          end do
        end do
 
        return
        end
 
C Head Integrate power spectrum
 
        subroutine power(psc, Ns, pwr, Np, iframe, Maxcomp, Maxap,
     1          s_int, s_pos, c_pos, c_rad, c_typ, pdata1, pdata2,
     2                   image1, image2, ary1, ary2, cary1, cary2)
*******************************************************************************
* Produces fringe/speckle patterns of multiple point sources.
* Finite bandwidth effect is taken into account.
*******************************************************************************
 
        implicit none
        integer Ns,Np,iframe,Maxcomp,Maxap
 
* Source information (intensity and position in arcsec)
        real   s_int(Maxcomp), s_pos(Maxcomp,2)
 
* Aperture information (coordinate, radius, type)
        real   c_pos(Maxap,2), c_rad(Maxap)
        integer   c_typ(Maxap)
 
* Phase screen
        complex psc(Ns,Ns)
 
* Power spectrum
        real   pwr(Np,Np)
 
* Temporary arrays
        real    image1(Np,Np), image2(Np,Np),  ary1(Np,Np), ary2(Np,Np)
        complex pdata1(Np,Np), pdata2(Np,Np), cary1(Np,Np), cary2(Np,Np)
 
* Common block
        real    pi,rapdeg,rapmin,rapsec, Ls, Lp
        complex iu
        real    k, kapmax, kapmin, L0, LCn2, dkap, lambda, fraction
        real    x_min,x_max,xstep, y_min,y_max,ystep
        integer ncomp,ncircle,nsamp,nframe
*
        common/psc/ pi,rapdeg,rapmin,rapsec, iu,
     &              k, kapmax, kapmin, L0, LCn2, Ls, Lp,
     &              dkap, lambda, nframe,
     &              x_min, x_max, xstep, y_min, y_max, ystep,
     &              ncomp, ncircle, nsamp, fraction
 
*
        integer N(2),Wc,Npx
        real    lpix,Scale,factor
        real    x,y,xa,ya,ox,oy
        integer ix,iy,jx,jy,lx,ly,dx,dy,iox,ioy
        integer loop1,loop2,i
        character comment*80
 
        N(1) = Np
        N(2) = Np
        call four2(pdata1,Np,Np,0)
        lpix = Ls/Ns
        Wc = Np/2+1
        Scale = Ls/lambda*Np/Ns*rapsec
 
* Setting apertures on Np sized window in Ns sized complex phase screen
 
        do 10 x=x_min,x_max,xstep
        do 20 y=y_min,y_max,ystep
**************************** major loop ***************************************
 
           write(1,111) x,y
C JLP 91        call easyou(0)
  111      format(' ',' creating a pair of speckle/fringe patterns',1x,
     1            'xy',f10.3,1x,f10.3)
 
* Reset temporary arrays
* CMV REV
           do jy=1,Np
             do jx=1,Np
               image1(jx,jy) = 0.0
               image2(jx,jy) = 0.0
                 ary1(jx,jy) = 0.0
                 ary2(jx,jy) = 0.0
                cary1(jx,jy) = (0.0,0.0)
                cary2(jx,jy) = (0.0,0.0)
             end do
           end do
 
C Head Loop over waveband
 
           do 100 loop1=1,nsamp
 
* Reset pdata1 and pdata2
* CMV REV
             do iy=1, Np
               do ix=1, Np
                 pdata1(ix,iy) = (0.0,0.0)
                 pdata2(ix,iy) = (0.0,0.0)
               end do
             end do
 
**************** intermediate loop *****************
* nsamp data points within the bandpass
             factor =  1.0 + fraction*float(loop1 - nsamp/2)
 
             npx = 0
 
             do 30 i=1,ncircle
*-----------------------------------------
* Filling pdata with aperture phases
 
            do 40 xa = x+c_pos(i,1)-c_rad(i), x+c_pos(i,1)+c_rad(i), lpix
            do 50 ya = y+c_pos(i,2)-c_rad(i), y+c_pos(i,2)+c_rad(i), lpix
              ox = xa-x-c_pos(i,1)
              oy = ya-y-c_pos(i,2)
              if( ox*ox+oy*oy .lt. c_rad(i)*c_rad(i) ) then
                ix = nint( (xa-x)/lpix )+Wc
                iy = nint( (ya-y)/lpix )+Wc
 
* Shift origin of the window to (Wc,Wc)
                lx = nint(xa/lpix)
                ly = nint(ya/lpix)
 
                if( c_typ(i) .eq. 1) npx = npx+1
                if( c_typ(i) .eq. 0) npx = npx-1
 
                if( npx.eq.1 .or. mod(npx,100) .eq. 0 ) then
                  write(comment,777) npx
                  call status(comment)
                end if
  777   format('# of  positive pixels',1x,i5)
 
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
 
* Simulate propagation from aperture plane to detector plane
 
        call four2(pdata1,Np,Np,-1)
        call four2(pdata2,Np,Np,-1)
 
* Form speckle/fringe patterns
* CMV REV
        do iy=1, Np
          do ix=1, Np
            image1(ix,iy) = image1(ix,iy) +
     1                      real( pdata1(ix,iy)*conjg(pdata1(ix,iy)) )
            image2(ix,iy) = image2(ix,iy) +
     2                      real( pdata2(ix,iy)*conjg(pdata2(ix,iy)) )
          end do
        end do
 
******************* end of loop1 ********************
  100   continue
 
C Head Put in sources and make power spectrum
 
        do 60 loop2=1,ncomp
*--------------------------------------------
*       ncomp components in source structure
          dx = int(Scale * s_pos(loop2,1))
          dy = int(Scale * s_pos(loop2,2))
 
          iox=0
          ioy=0
          if (dx.lt.0) iox=dx
          if (dy.lt.0) ioy=dy
* CMV REV
          do iy=1-ioy,Np-dy+ioy
          do ix=1-iox,Np-dx+iox
            ary1(ix+dx,iy+dy) = ary1(ix,iy)+s_int(loop2)*image1(ix,iy)
            ary2(ix+dx,iy+dy) = ary2(ix,iy)+s_int(loop2)*image2(ix,iy)
          end do
          end do
*-------------------------------------------
   60   continue
 
        write(comment,320) iframe
 
* Shift origin of the image to (Wc,Wc) before fourier transform
* CMV REV
        do iy=1,Np,2
        do ix=1,Np,2
          cary1(ix,iy)     = ary1(ix,iy)
          cary1(ix+1,iy+1) = ary1(ix+1,iy+1)
          cary1(ix+1,iy)   =-ary1(ix+1,iy)
          cary1(ix,iy+1)   =-ary1(ix,iy+1)
        end do
        end do
 
        call four2(cary1,Np,Np,-1)
 
* CMV REV
        do jy=1,Np
        do jx=1,Np
          pwr(jx,jy) = pwr(jx,jy) + real(cary1(jx,jy)*conjg(cary1(jx,jy)))
        end do
        end do
 
        iframe = iframe + 1
        if( iframe .gt. nframe ) return
 
        write(comment,320) iframe
  320   format(' ', 'processing ',i4,' th frame                    ')
 
* Shift origin of the image to (Wc,Wc) before fourier transform
* CMV REV
        do iy=1,Np,2
        do ix=1,Np,2
          cary2(ix,iy)     = ary2(ix,iy)
          cary2(ix+1,iy+1) = ary2(ix+1,iy+1)
          cary2(ix+1,iy)   =-ary2(ix+1,iy)
          cary2(ix,iy+1)   =-ary2(ix,iy+1)
        end do
        end do
 
        call four2(cary2,Np,Np,-1)
*CMV REV
        do jy=1,Np
        do jx=1,Np
          pwr(jx,jy) = pwr(jx,jy) + real(cary2(jx,jy)*conjg(cary2(jx,jy)))
        end do
        end do
 
        iframe = iframe + 1
        if( iframe .gt. nframe) return
 
************************* end of major loop ***********************************
   20   continue
   10   continue
 
        end
 
CHead Subroutine to write phase screen
 
        subroutine psc_wscreen(array,dx,dy,rbuf,ibuf,nx,ny)
c
        implicit  none
        integer   dx,dy,nx,ny
        complex   array(dx,dy)
        real      rbuf(nx,ny),ibuf(nx,ny)
        character ofile*25,otype*10,oform*5
        integer   err,ix,iy,x0,y0
        logical   first
        data      first/.TRUE./
        save      first
c
        err=0
        if (first) then
           first=.false.
           ofile=' '
           otype='GHRIL'
           oform='-32'
           call usertext(ofile,1,'SFILE=',
     &                  'Enter filename for phase-screen')
           call usertext(otype,2,'STYPE=','Enter file-type')
           call usertext(oform,2,'SFORM=','Enter format for GHRIL-file')
           call data_open(1,ofile,otype,'WRITE',err)
           call data_wkey(1,'OBJECT','Screen',' ',err)
           call date_header(1,'DATE','TIME')
        end if
*
*       Write each screen in eight parts: re,im ...
*
        do x0=0,dx-1,nx
          do y0=0,dy-1,nx
*CMV REV
             do iy=1,ny
               do ix=1,nx
                  rbuf(ix,iy)= real(array(x0+ix,y0+iy))
                  ibuf(ix,iy)=aimag(array(x0+ix,y0+iy))
               end do
             end do
             call data_write(1,0,0,1,nx,ny,rbuf,0,0,nx,oform,err)
             call data_write(1,0,0,1,nx,ny,ibuf,0,0,nx,oform,err)
          end do
        end do
        call data_error(err)
        return
        end
 
CHead Subroutine to write image plane
 
        subroutine psc_wimage(array1,array2,dx,dy)
c
        implicit  none
        real      array1(*),array2(*)
        integer   dx,dy
        character ofile*25,otype*10,oform*5
        integer   err
        logical   first
        data      first/.TRUE./
        save      first
c
        err=0
        if (first) then
           first=.false.
           ofile=' '
           otype='GHRIL'
           oform='-32'
           call usertext(ofile,1,'IFILE=',
     &                  'Enter filename for image-plane')
           call usertext(otype,2,'ITYPE=','Enter file-type')
           call usertext(oform,2,'SFORM=','Enter format for GHRIL-file')
           call data_open(2,ofile,otype,'WRITE',err)
           call data_wkey(2,'OBJECT','Image',' ',err)
           call date_header(2,'DATE','TIME')
        end if
        call data_write(2,0,0,1,dx,dy,array1,0,0,dx,oform,err)
        call data_write(2,0,0,1,dx,dy,array2,0,0,dx,oform,err)
        call data_error(err)
        return
        end
 
CHead Subroutine to write phase screen
 
        subroutine psc_wpower(array,dx,dy)
c
        implicit  none
        real      array(*)
        integer   dx,dy
        character ofile*25,otype*10,oform*5
        integer   err
        logical   first
        data      first/.TRUE./
        save      first
c
        err=0
        if (first) then
           first=.false.
           ofile=' '
           otype='GHRIL'
           oform='-32'
           call usertext(ofile,1,'PFILE=',
     &                  'Enter filename for power spectrum')
           call usertext(otype,2,'PTYPE=','Enter file-type')
           call usertext(oform,2,'SFORM=','Enter format for GHRIL-file')
           call data_open(3,ofile,otype,'WRITE',err)
           call data_wkey(3,'OBJECT','Power',' ',err)
           call date_header(3,'DATE','TIME')
        end if
        call data_write(3,0,0,1,dx,dy,array,0,0,dx,oform,err)
        call data_error(err)
        return
        end
