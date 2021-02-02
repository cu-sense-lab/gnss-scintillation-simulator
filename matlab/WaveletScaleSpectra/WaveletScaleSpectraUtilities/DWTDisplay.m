dx=Tseg; nsamp=length(I);
[DWT_Spec,wc,J,yp,qmf]=ComputeDWT(I,'Symmlet');
%Compute spatial wavenumber scales & extension to nearest power of 2
sj=(2*dx)*2.^(J-2:-1:0);
qScale=1./sj;
%Extended scale for nearest power of 2
[~,nsamp_ext]=size(DWT_Spec);
x_data_ext=dx*(0:nsamp_ext-1);
xscale    =x_data_ext(1:nsamp);
if DISP~=0
    MaxOffset=10; DynRange=100;
    DWT_max=10*ceil(max(dB10(DWT_Spec(:)))/10)-MaxOffset;
    nfig=figure;   %Display DWT segmentation
    imagesc(xscale,log10(qScale),dB10(DWT_Spec(:,1:nsamp)))
    grid on
    xlabel('span')
    caxis([DWT_max-DynRange, DWT_max])
    colorbar
    title([PLOTID, ' DWT'])
end
