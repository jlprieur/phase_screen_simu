#!/bin/csh
###################################################################
# simu1000 Command procedure for simulations
#
# JLP
# Version 24-05-2008
###################################################################
setenv JLP_FORMAT 8,8
unalias lpr
alias lpr lpr -h
unalias rm
alias rm rm -f
#
# BEGIN:
goto binary
#
#goto skip0
#goto profile 
#goto inv1 
#
# seed = s
# n = nframes
# p = photons
# s = seed for random generator
# r = Kolmogorov radius 
# -a for output of aperture functions and psf's for frames 1 to 4 
#
#double star:
####################################################################
#runs object_cerga
# Enter output size  nx,ny :128,128
# Enter output file name: dble
# Number of objects in the output image = 2
# Object # 1
# Analytical object (1) or pixel function (2) ? 2
# Coordonnees centrales x,y (centre fenetre = (0,0)) = -10,-10
# Brillance du pixel (float) = 2.
# Object # 2
# Analytical object (1) or pixel function (2) ? 2
# Coordonnees centrales x,y (centre fenetre = (0,0)) = 10,10
#   Brillance du pixel (float) = 1.
####################################################################
binary:
$EXEC/simu3.exe -n 1000 -p 130 -o dble -s 115 -r 0.1 -a
$EXEC/decode_simu2.exe simu2.data1 nphotons 12,100 4000
goto end
# Single star:
# $EXEC/simu3.exe -n 1000 -p 1000 -s 123 -r 0.2 -a
$EXEC/simu3.exe -n 10 -p 1000 -s 123 -r 0.2 -a
goto end 
runs pstcopy0 endfrm
lpr pstcopy.tmp
runs pstcopy0 aper1
lpr pstcopy.tmp
runs pstcopy0 psf1
lpr pstcopy.tmp
skip1:
# Generating b_real, b_modsq, b_bisp,...
runs bispect 12 obj1
goto end
inv1:
#bisp1_100: r0=0.1 100frames, 1000photons
#bisp1_200: r0=0.2 300frames, 1000photons
#bisp1_300: r0=0.2 1000frames, 1000photons
# 20,0.1 initially
#runs inv_bispec 2 12,20,0.008,220,0,5. modsq_100 b_real b_imag bisp1_100
# 2= simulations 
# 12= Radius of uv-coverage
# 0= weights set to unity
# 0.008= test to exit from main loop
# 20= number of bispectral terms
# 0= lower cut to take into account spectral terms
rm err_mod.dat
rm sigma.dat
rm error1.dat
rm error2.dat
rm bisp_error1.dat
rm bisp_error2.dat
rm inv_bispec.log
# 20= (old) weights: inv.  prop. to radius 
#runs inv_bispec 2 12,0,0.008,220,0.,0.7 b_modsq b_real b_imag bisp1_100
runs inv_bispec 2 12,0,0.008,220,0.,0.8 modsq_1000 b_real b_imag bisp1_1000
#runs inv_bispec 2 12,20,0.008,220,0.,0.7 b_modsq b_real b_imag bisp1_200
#runs inv_bispec 0 12,20,0.008,220,0.,0.7 b_real b_imag 
runs fft2 -1 ref n imf test_100 testi
goto end
runs plot1 < plot_err2.com
runs plot1 < plot_err1.com
#print inv_bispec.log
runs fft2 -1 rei n imi testi_100 testi
#
goto fft
mv error2.dat error_100.dat
lpr inv_bispec.log
#goto end
runs bisp_ima 12 bisp2 weight_100
runs pstcopy0 weight_100 /AUTOSCALE
lpr pstcopy.tmp
goto end
fft:
runs pstcopy0 test_100
lpr pstcopy.tmp
runs fft2 -1 rei n imi testi_100 testi
runs pstcopy0 testi_100
lpr pstcopy.tmp
#
profile:
rm profile1.log
rm modsq.pro
cp modsq_200.mt modsq.mt
runs profile1 N modsq_centred.par
$EXEC/plot1.exe < plot1.com
end:
