      programs_pscorg
      real*4pi,twopi
      real*4lcn2
      real*4l_s,l_p
      parameter(l_s=10.24)
      parameter(l_p=0.02)
      integer*4nsise
      parameter(nsize=256)
      realkk,kapmax,kapmin,l_0,dkap,lambda,l_o
      real*4pha_re(nsize,nsize),pha_im(nsize,nsize)
      character*30ifile,ofile
      integer*4f_in
      parameter(f_in=11)
      character*80comment,ocomments
      integer*4seed
      pi=3.141592653
      twopi=2.*pi
      lcn2=5.0e-13
      calljlp_begin
      calljlp_inquifmt
      write(*,10)
   10 format(1x,' input ascii file name : ',$)
      read(*,'(a)')ifile
      open(unit=f_in,file=ifile,status='old')
   30 format(a80)
   35 format(1x,a80)
      do40j=1,100
      read(f_in,30,end=100)comment
      write(*,35)comment
      if(comment(1:1).eq.'!')goto40
      if(comment(1:1).eq.'R'.or.comment(1:1).eq.'r')then
      read(f_in,*)seed
      endif
      if(comment(1:1).eq.'T'.or.comment(1:1).eq.'t')then
      read(f_in,*)l_0,l_o
      endif
      if(comment(1:1).eq.'W'.or.comment(1:1).eq.'w')then
      read(f_in,*)lambda
      endif
   40 continue
  100 continue
      write(6,36)seed,l_0,l_o,lambda
      write(ocomments,36)seed,l_0,l_o,lambda
   36 format(' seed=',I6,' L0 = ',G10.4,' lo = ',G10.4,' Lambda = ',G10.
     +4)
      kk=twopi/lambda
      work=min(l_s,l_0)
      kapmin=twopi/work
      work=max(l_o,l_p)
      kapmax=twopi/work
      dkap=twopi/l_s
      calljlp_random_init(seed)
      callscreen(pha_re,pha_im,nsize,kk,kapmax,kapmin,lcn2,dkap)
      ofile=' '
      print*,' output of phase screen (real part)'
      calljlp_writeimag(pha_re,nsize,nsize,nsize,ofile,ocomments)
      ofile=' '
      print*,' output of phase screen (imaginary part)'
      calljlp_writeimag(pha_im,nsize,nsize,nsize,ofile,ocomments)
      calljlp_end
      stop
      end
      subroutinescreen(pha_re,pha_im,nsize,kk,kapmax,kapmin,lcn2,dkap)
      real*4pha_re(nsize,nsize),pha_im(nsize,nsize)
      real*4fs,fso,pi,twopi
      real*4kk,kapmax,kapmin,lcn2,dkap,u,wrand,work
      real*4kap
      integer*4Nc
      pi=3.141592653
      twopi=2.*pi
      Nc=nsize/2+1
      fs0=0.033*pi*lcn2*kk*kk*2.0*dkap*dkap
      do100jy=1,nsize
      do110jx=1,nsize
      kap=dkap*sqrt(float((jx-Nc)*(jx-Nc)+(jy-Nc)*(jy-Nc)))
      kap=max(kap,dkap)
      atten=1.0
      if(kap.gt.kapmax)then
      atten=exp(-(kap-kapmax)/(5.0*dkap))
      elseif(kap.lt.kapmin)then
      atten=exp(-(kapmin-kap)/(5.0*dkap))
      endif
   10 calljlp_random(u)
      if(u.le.0.0)goto10
      radi=sqrt(-2.0*log(u))
      fs=fs0*kap**(-11.0/3.0)
      work=atten*sqrt(fs)*radi
      calljlp_random(wrand)
      wrand=wrand*twopi
      pha_re(jx,jy)=work*cos(wrand)
      pha_im(jx,jy)=work*sin(wrand)
  110 continue
  100 continue
      write(*,555)nsize,nsize
  555 format('+','start ',i3,'x',i3,' sized fft              ')
      kod=1
      callfft_2d(pha_re,pha_im,nsize,nsize,nsize,kod)
      write(*,666)
  666 format('+','end fft for screen                         ')
      DO80000jy=1,nsize
      DO80001jx=1,nsize
      work=pha_re(jx,jy)
      pha_re(jx,jy)=cos(work)
      pha_im(jx,jy)=sin(work)
80001 CONTINUE
80000 CONTINUE
      return
      end
