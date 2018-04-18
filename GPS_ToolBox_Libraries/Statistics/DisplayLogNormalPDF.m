% LogNormal PDF & CDF
%
clear
newColors=colors;
blk{1}='k-';
blk{2}='k-.';
blk{3}='k--';
blk{4}='kx-';
blk{5}='k*-';
blk{6}='k+-';
blk{7}='kp-';
I=linspace(0.01,10,500);
I_dB=dB10(I);

figure
SI=0.25:0.25:1.75;
for m=1:length(SI)
    plot( 0, 10, 'Color',newColors{2*m,2});
    hold on
    PDF=LogNormalPDF(I,SI(m));
    subplot(2,1,1)
    plot( I_dB, PDF, 'Color',newColors{2*m,2});
    %plot( I_dB, PDF, blk{m});
    hold on
    subplot(2,1,2)
    for ii=1:length(PDF)
        CDF(ii)=sum(PDF(1:ii));
    end
    CDF=CDF.*[0,diff(I)];
    plot( I_dB, CDF, 'Color',newColors{2*m,2});
    %plot( I_dB, CDF,blk{m});
    hold on
end
subplot(2,1,1)
legend('0.25','0.5','0.75','1.0','1.25','1.5','1.75')
grid on
title('LogNormal')
ylabel('PDF')
axis([-20 10 0 2])
subplot(2,1,2)
grid on
ylabel('CDF')
xlabel('I/<I> - dB')
axis([-20 10 0 1])