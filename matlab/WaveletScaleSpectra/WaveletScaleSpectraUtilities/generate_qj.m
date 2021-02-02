function qj=generate_qj(dy,J,nSmax)
%
sj=(2*dy)*2.^(J-2:-1:0);   %DWT range
sj=sj(nSmax+1:end);      %Scale Spectrum
qj=1./sj;
return

