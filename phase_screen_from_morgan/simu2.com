# Command procedure for simulations
#
unalias lpr 
alias lpr lpr -h 
unalias rm
alias rm rm -f
#
# BEGIN:
#
goto inv1:
#goto skip2 
#
# seed = s
$EXEC/simu2.exe -n 10 -p 200 -o object.bdf -s 1234 -r 8 -a
goto end
runs pstcopy0 endfrm
lpr pstcopy.tmp
goto skip1
runs pstcopy0 aper1
lpr pstcopy.tmp
runs pstcopy0 psf1
lpr pstcopy.tmp
skip1:
runs pstcopy0 obj_conv
lpr pstcopy.tmp
# Generating b_real, b_modsq, b_bisp,...
runs bispect 12 obj_conv
skip2:
rm jlp_lu5.tmp 
echo 12 > jlp_lu5.tmp
$EXEC/decode_simu2.exe < jlp_lu5.tmp
rm jlp_lu5.tmp
#goto move
runs pstcopy0 long 
lpr pstcopy.tmp
move:
cp simu2.data1 simu2_si.data1
mv long.bdf long_si.bdf
mv modsq.bdf modsq_si.bdf
mv bisp1.bdf bisp1_si.bdf
inv1:
# 20,0.1 initially
runs inv_bispec 2 12,60,0.1,0.000001 modsq_si b_real b_imag bisp1_si
goto end
lpr inv_bispec.log
fft:
runs fft2 -1 ref n imf test_si testi
runs pstcopy0 test_si
lpr pstcopy.tmp
end:
