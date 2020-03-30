      programs_psc
      realpi,rapdeg,rapmin,rapsec
      complexiu
      realLCn2
      realLs,Lp
      integerNs,Np
      parameter(Ns=512)
      parameter(Np=Ns/2)
      realk,kapmax,kapmin,L0,dkap,lo
      reallambda,fraction
      integernsamp
      complexpsc(Ns,Ns)
      realpwr(Np,Np)
      integerMaxcomp,Maxap
      parameter(Maxcomp=50,Maxap=15)
      reals_int(Maxcomp),s_pos(Maxcomp,2)
      realc_pos(Maxap,2),c_rad(Maxap)
      integerc_typ(Maxap)
      realx_min,x_max,xstep,y_min,y_max,ystep
      integerncomp,ncircle,nframe
      common/psc/pi,rapdeg,rapmin,rapsec,iu,k,kapmax,kapmin,L0,LCn2,Ls,L
     +p,dkap,lambda,nframe,x_min,x_max,xstep,y_min,y_max,ystep,ncomp,nci
     +rcle,nsamp,fraction
      realimage1(Np,Np),image2(Np,Np),ary1(Np,Np),ary2(Np,Np)
      complexpdata1(Np,Np),pdata2(Np,Np),cary1(Np,Np),cary2(Np,Np)
      realzero
      integerseed,n,iscreen,iframe,iz,jz,lx,ly
      logicalno_dist,w_screen,w_image,w_power
      charactercomment*80
      pi=3.141592653
      rapdeg=pi/180.0
      rapmin=rapdeg/60.0
      rapsec=rapdeg/3600.0
      iu=(0.0,1.0)
      LCn2=5.0e-13
      Ls=10.24
      Lp=0.02
      callinit()
      nframe=50
      x_min=2.5
      y_min=2.5
      x_max=8.5
      y_max=8.5
      xstep=.15
      ystep=8.0
      calluserint(nframe,1,1,'NFRAME=','Enter the number of frames')
      calluserrpair(x_min,y_min,1,'SUB_MIN=','Enter minimum X,Y of subfr
     +ame (m)')
      calluserrpair(x_max,y_max,1,'SUB_MAX=','Enter maximum X,Y of subfr
     +ame (m)')
      calluserrpair(xstep,ystep,1,'SUB_STEP=','Enter stepsize in X,Y (m)
     +')
      L0=0.60
      lo=0.02
      calluserrpair(L0,lo,2,'SCALE=','Enter outer,inner scale of turbule
     +nce (m)')
      lambda=0.63e-6
      fraction=0.00476
      nsamp=3
      calluserreal(lambda,1,2,'LAMBDA=','Enter wavelength (m)')
      calluserreal(fraction,1,2,'FBAND=','Enter fractional bandwidth')
      calluserint(nsamp,1,2,'NSAMP=','Enter number of sampling points')
      ncircle=1
      c_rad(1)=1.
      c_pos(1,1)=0.0
      c_pos(1,2)=0.0
      c_typ(1)=1
      calluserint(ncircle,1,1,'NCIRCLE=','Enter number of apertures')
      DO80000n=1,ncircle
      calluserrpair(c_pos(n,1),c_pos(n,2),1,'CPOS=','Enter position in X
     +,Y (m)')
      calluserreal(c_rad(n),1,1,'CRAD=','Enter radius (m)')
      calluserint(c_typ(n),1,1,'CTYPE=','Enter type (1=hole,0=void)')
      c_rad(n+1)=c_rad(n)
      c_pos(n+1,1)=c_pos(n,1)
      c_pos(n+1,2)=c_pos(n,1)
      c_typ(n+1)=c_typ(n)
      callcancel('CPOS=')
      callcancel('CRAD=')
      callcancel('CTYPE=')
80000 CONTINUE
      ncomp=1
      calluserint(ncomp,1,1,'NCOMP=','Enter number of sources')
      DO80001n=1,ncomp
      s_int(n)=1.
      s_pos(n,1)=0.0
      s_pos(n,2)=0.0
      calluserrpair(s_pos(n,1),s_pos(n,2),1,'SPOS=','Enter position in X
     +,Y (m)')
      calluserreal(s_int(n),1,1,'SINT=','Enter intensity')
      callcancel('SPOS=')
      callcancel('SINT=')
80001 CONTINUE
      no_dist=.FALSE.
      w_screen=.FALSE.
      w_image=.FALSE.
      w_power=.TRUE.
      calluserlog(no_dist,1,2,'NDIST=','No disturbance?')
      calluserlog(w_screen,1,2,'WSCREEN=','Write screen data?')
      calluserlog(w_image,1,2,'WIMAGE=','Write image plane?')
      calluserlog(w_power,1,2,'WPOWER=','Write power spectrum?')
      seed=7597551
      calluserint(seed,1,2,'SEED=','Enter seed for random generator')
      k=2.0*pi/lambda
      if(Ls.lt.L0)then
      kapmin=2.0*pi/Ls
      else
      kapmin=2.0*pi/L0
      endif
      if(lo.lt.Lp)then
      kapmax=2.0*pi/Lp
      else
      kapmax=2.0*pi/lo
      endif
      dkap=2.0*pi/Ls
      write(1,*)'frm_per_scr ',(int((x_max-x_min)/xstep)+1)*(int((y_max-
     +y_min)/ystep)+1)
      iframe=1
      iscreen=1
  287 if(no_dist)then
      DO80002iz=1,Ns
      DO80003jz=1,Ns
      psc(iz,jz)=(0.0,0.0)
80003 CONTINUE
80002 CONTINUE
      else
      write(comment,310)iscreen
  310 format(' Generating ',i3,'th complex phase screen')
      callscreen(psc,Ns,seed)
      iscreen=iscreen+1
      if(w_screen)callpsc_wscreen(psc,Ns,Ns,pdata1,pdata2,Np,Np)
      endif
      callpower(psc,Ns,pwr,Np,iframe,Maxcomp,Maxap,s_int,s_pos,c_pos,c_r
     +ad,c_typ,pdata1,pdata2,image1,image2,ary1,ary2,cary1,cary2)
      if(w_image)callpsc_wimage(ary1,ary2,Np,Np)
      if(iframe.le.nframe)goto287
      zero=pwr(Np/2+1,Np/2+1)
      if(zero.eq.0.0)zero=1.0
      write(1,*)'Scaling: divide by ',zero
      DO80004ly=1,Np
      DO80005lx=1,Np
      pwr(lx,ly)=pwr(lx,ly)/zero
80005 CONTINUE
80004 CONTINUE
      print*,'Scaled'
      if(w_power)callpsc_wpower(pwr,Np,Np)
      n=0
      calldata_close(1,n)
      n=0
      calldata_close(2,n)
      n=0
      calldata_close(3,n)
      callfinis()
      end
      subroutinescreen(psc,Ns,seed)
      implicitnone
      integerNs,seed
      complexpsc(Ns,Ns)
      realpi,rapdeg,rapmin,rapsec,Ls,Lp
      complexiu
      realk,kapmax,kapmin,L0,LCn2,dkap,lambda,fraction
      realx_min,x_max,xstep,y_min,y_max,ystep
      integerncomp,ncircle,nsamp,nframe
      common/psc/pi,rapdeg,rapmin,rapsec,iu,k,kapmax,kapmin,L0,LCn2,Ls,L
     +p,dkap,lambda,nframe,x_min,x_max,xstep,y_min,y_max,ystep,ncomp,nci
     +rcle,nsamp,fraction
      realkap
      integerN(2),Nc
      integeri,j,jx,jy
      realFs,Fs0,atten,radi,u,ran
      charactercomment*80
      N(1)=Ns
      N(2)=Ns
      callfour2(psc,Ns,Ns,0)
      Nc=Ns/2+1
      Fs0=0.033*pi*LCn2*k*k*2.0*dkap*dkap
      do100jy=1,Ns
      if(mod(jy,100).eq.0)then
      write(comment,*)'Making line ',jy
      callstatus(comment)
      endif
      do110jx=1,Ns
      kap=dkap*sqrt(float((jx-Nc)**2+(jy-Nc)**2))
      if(kap.gt.kapmax)atten=exp(-(kap-kapmax)/(5.0*dkap))
      if(kap.lt.kapmin)atten=exp(-(kapmin-kap)/(5.0*dkap))
      if(kap.le.kapmax.and.kap.ge.kapmin)atten=1.0
      if(kap.eq.0.0)kap=dkap
      u=0.0
  288 u=ran(seed)
      if(u.le.0.0)goto288
      radi=sqrt(-2.0*log(u))
      Fs=Fs0*kap**(-11.0/3.0)
      psc(jx,jy)=atten*sqrt(Fs)*radi*cexp(iu*2.0*pi*ran(seed))
  110 continue
  100 continue
      DO80006j=1,Ns,2
      DO80007i=1,Ns,2
      psc(i+1,j)=-psc(i+1,j)
      psc(i,j+1)=-psc(i,j+1)
80007 CONTINUE
80006 CONTINUE
      write(comment,555)N(1),N(2)
  555 format(' ','start ',i3,'x',i3,' sized fft              ')
      callfour2(psc,Ns,Ns,-1)
      write(comment,666)
  666 format(' ','end fft for screen                         ')
      DO80008j=1,Ns,2
      DO80009i=1,Ns,2
      psc(i+1,j)=-psc(i+1,j)
      psc(i,j+1)=-psc(i,j+1)
80009 CONTINUE
80008 CONTINUE
      return
      end
      subroutinepower(psc,Ns,pwr,Np,iframe,Maxcomp,Maxap,s_int,s_pos,c_p
     +os,c_rad,c_typ,pdata1,pdata2,image1,image2,ary1,ary2,cary1,cary2)
      implicitnone
      integerNs,Np,iframe,Maxcomp,Maxap
      reals_int(Maxcomp),s_pos(Maxcomp,2)
      realc_pos(Maxap,2),c_rad(Maxap)
      integerc_typ(Maxap)
      complexpsc(Ns,Ns)
      realpwr(Np,Np)
      realimage1(Np,Np),image2(Np,Np),ary1(Np,Np),ary2(Np,Np)
      complexpdata1(Np,Np),pdata2(Np,Np),cary1(Np,Np),cary2(Np,Np)
      realpi,rapdeg,rapmin,rapsec,Ls,Lp
      complexiu
      realk,kapmax,kapmin,L0,LCn2,dkap,lambda,fraction
      realx_min,x_max,xstep,y_min,y_max,ystep
      integerncomp,ncircle,nsamp,nframe
      common/psc/pi,rapdeg,rapmin,rapsec,iu,k,kapmax,kapmin,L0,LCn2,Ls,L
     +p,dkap,lambda,nframe,x_min,x_max,xstep,y_min,y_max,ystep,ncomp,nci
     +rcle,nsamp,fraction
      integerN(2),Wc,Npx
      reallpix,Scale,factor
      realx,y,xa,ya,ox,oy
      integerix,iy,jx,jy,lx,ly,dx,dy,iox,ioy
      integerloop1,loop2,i
      charactercomment*80
      N(1)=Np
      N(2)=Np
      callfour2(pdata1,Np,Np,0)
      lpix=Ls/Ns
      Wc=Np/2+1
      Scale=Ls/lambda*Np/Ns*rapsec
      do10x=x_min,x_max,xstep
      do20y=y_min,y_max,ystep
      write(1,111)x,y
  111 format(' ',' creating a pair of speckle/fringe patterns',1x,'xy',f
     +10.3,1x,f10.3)
      DO80010jy=1,Np
      DO80011jx=1,Np
      image1(jx,jy)=0.0
      image2(jx,jy)=0.0
      ary1(jx,jy)=0.0
      ary2(jx,jy)=0.0
      cary1(jx,jy)=(0.0,0.0)
      cary2(jx,jy)=(0.0,0.0)
80011 CONTINUE
80010 CONTINUE
      do100loop1=1,nsamp
      DO80012iy=1,Np
      DO80013ix=1,Np
      pdata1(ix,iy)=(0.0,0.0)
      pdata2(ix,iy)=(0.0,0.0)
80013 CONTINUE
80012 CONTINUE
      factor=1.0+fraction*float(loop1-nsamp/2)
      npx=0
      do30i=1,ncircle
      do40xa=x+c_pos(i,1)-c_rad(i),x+c_pos(i,1)+c_rad(i),lpix
      do50ya=y+c_pos(i,2)-c_rad(i),y+c_pos(i,2)+c_rad(i),lpix
      ox=xa-x-c_pos(i,1)
      oy=ya-y-c_pos(i,2)
      if(ox*ox+oy*oy.lt.c_rad(i)*c_rad(i))then
      ix=nint((xa-x)/lpix)+Wc
      iy=nint((ya-y)/lpix)+Wc
      lx=nint(xa/lpix)
      ly=nint(ya/lpix)
      if(c_typ(i).eq.1)npx=npx+1
      if(c_typ(i).eq.0)npx=npx-1
      if(npx.eq.1.or.mod(npx,100).eq.0)then
      write(comment,777)npx
      callstatus(comment)
      endif
  777 format('# of  positive pixels',1x,i5)
      if(mod(ix+iy,2).eq.0)then
      pdata1(ix,iy)=cexp(iu*factor*real(psc(lx,ly)))*c_typ(i)
      pdata2(ix,iy)=cexp(iu*factor*aimag(psc(lx,ly)))*c_typ(i)
      else
      pdata1(ix,iy)=-cexp(iu*factor*real(psc(lx,ly)))*c_typ(i)
      pdata2(ix,iy)=-cexp(iu*factor*aimag(psc(lx,ly)))*c_typ(i)
      endif
      endif
   50 continue
   40 continue
   30 continue
      callfour2(pdata1,Np,Np,-1)
      callfour2(pdata2,Np,Np,-1)
      DO80014iy=1,Np
      DO80015ix=1,Np
      image1(ix,iy)=image1(ix,iy)+real(pdata1(ix,iy)*conjg(pdata1(ix,iy)
     +))
      image2(ix,iy)=image2(ix,iy)+real(pdata2(ix,iy)*conjg(pdata2(ix,iy)
     +))
80015 CONTINUE
80014 CONTINUE
  100 continue
      do60loop2=1,ncomp
      dx=int(Scale*s_pos(loop2,1))
      dy=int(Scale*s_pos(loop2,2))
      iox=0
      ioy=0
      if(dx.lt.0)iox=dx
      if(dy.lt.0)ioy=dy
      DO80016iy=1-ioy,Np-dy+ioy
      DO80017ix=1-iox,Np-dx+iox
      ary1(ix+dx,iy+dy)=ary1(ix,iy)+s_int(loop2)*image1(ix,iy)
      ary2(ix+dx,iy+dy)=ary2(ix,iy)+s_int(loop2)*image2(ix,iy)
80017 CONTINUE
80016 CONTINUE
   60 continue
      write(comment,320)iframe
      DO80018iy=1,Np,2
      DO80019ix=1,Np,2
      cary1(ix,iy)=ary1(ix,iy)
      cary1(ix+1,iy+1)=ary1(ix+1,iy+1)
      cary1(ix+1,iy)=-ary1(ix+1,iy)
      cary1(ix,iy+1)=-ary1(ix,iy+1)
80019 CONTINUE
80018 CONTINUE
      callfour2(cary1,Np,Np,-1)
      DO80020jy=1,Np
      DO80021jx=1,Np
      pwr(jx,jy)=pwr(jx,jy)+real(cary1(jx,jy)*conjg(cary1(jx,jy)))
80021 CONTINUE
80020 CONTINUE
      iframe=iframe+1
      if(iframe.gt.nframe)return
      write(comment,320)iframe
  320 format(' ','processing ',i4,' th frame                    ')
      DO80022iy=1,Np,2
      DO80023ix=1,Np,2
      cary2(ix,iy)=ary2(ix,iy)
      cary2(ix+1,iy+1)=ary2(ix+1,iy+1)
      cary2(ix+1,iy)=-ary2(ix+1,iy)
      cary2(ix,iy+1)=-ary2(ix,iy+1)
80023 CONTINUE
80022 CONTINUE
      callfour2(cary2,Np,Np,-1)
      DO80024jy=1,Np
      DO80025jx=1,Np
      pwr(jx,jy)=pwr(jx,jy)+real(cary2(jx,jy)*conjg(cary2(jx,jy)))
80025 CONTINUE
80024 CONTINUE
      iframe=iframe+1
      if(iframe.gt.nframe)return
   20 continue
   10 continue
      end
      subroutinepsc_wscreen(array,dx,dy,rbuf,ibuf,nx,ny)
      implicitnone
      integerdx,dy,nx,ny
      complexarray(dx,dy)
      realrbuf(nx,ny),ibuf(nx,ny)
      characterofile*25,otype*10,oform*5
      integererr,ix,iy,x0,y0
      logicalfirst
      datafirst/.TRUE./
      savefirst
      err=0
      if(first)then
      first=.false.
      ofile=' '
      otype='GHRIL'
      oform='-32'
      callusertext(ofile,1,'SFILE=','Enter filename for phase-screen')
      callusertext(otype,2,'STYPE=','Enter file-type')
      callusertext(oform,2,'SFORM=','Enter format for GHRIL-file')
      calldata_open(1,ofile,otype,'WRITE',err)
      calldata_wkey(1,'OBJECT','Screen',' ',err)
      calldate_header(1,'DATE','TIME')
      endif
      DO80026x0=0,dx-1,nx
      DO80027y0=0,dy-1,nx
      DO80028iy=1,ny
      DO80029ix=1,nx
      rbuf(ix,iy)=real(array(x0+ix,y0+iy))
      ibuf(ix,iy)=aimag(array(x0+ix,y0+iy))
80029 CONTINUE
80028 CONTINUE
      calldata_write(1,0,0,1,nx,ny,rbuf,0,0,nx,oform,err)
      calldata_write(1,0,0,1,nx,ny,ibuf,0,0,nx,oform,err)
80027 CONTINUE
80026 CONTINUE
      calldata_error(err)
      return
      end
      subroutinepsc_wimage(array1,array2,dx,dy)
      implicitnone
      realarray1(*),array2(*)
      integerdx,dy
      characterofile*25,otype*10,oform*5
      integererr
      logicalfirst
      datafirst/.TRUE./
      savefirst
      err=0
      if(first)then
      first=.false.
      ofile=' '
      otype='GHRIL'
      oform='-32'
      callusertext(ofile,1,'IFILE=','Enter filename for image-plane')
      callusertext(otype,2,'ITYPE=','Enter file-type')
      callusertext(oform,2,'SFORM=','Enter format for GHRIL-file')
      calldata_open(2,ofile,otype,'WRITE',err)
      calldata_wkey(2,'OBJECT','Image',' ',err)
      calldate_header(2,'DATE','TIME')
      endif
      calldata_write(2,0,0,1,dx,dy,array1,0,0,dx,oform,err)
      calldata_write(2,0,0,1,dx,dy,array2,0,0,dx,oform,err)
      calldata_error(err)
      return
      end
      subroutinepsc_wpower(array,dx,dy)
      implicitnone
      realarray(*)
      integerdx,dy
      characterofile*25,otype*10,oform*5
      integererr
      logicalfirst
      datafirst/.TRUE./
      savefirst
      err=0
      if(first)then
      first=.false.
      ofile=' '
      otype='GHRIL'
      oform='-32'
      callusertext(ofile,1,'PFILE=','Enter filename for power spectrum')
      callusertext(otype,2,'PTYPE=','Enter file-type')
      callusertext(oform,2,'SFORM=','Enter format for GHRIL-file')
      calldata_open(3,ofile,otype,'WRITE',err)
      calldata_wkey(3,'OBJECT','Power',' ',err)
      calldate_header(3,'DATE','TIME')
      endif
      calldata_write(3,0,0,1,dx,dy,array,0,0,dx,oform,err)
      calldata_error(err)
      return
      end
