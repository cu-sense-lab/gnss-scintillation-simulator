%Display Rician
p=linspace(0,100,500);
[err,n80]=min(abs(p-80));
figure
for b=[0,5,10,20,30,40]
    pdf=Rician(p,b);
    hold on
    plot(p,dB10(pdf)/10,'r')
    hold on
    text(p(n80),dB10(pdf(n80))/10,num2str(b))
    hold on
end
grid on
axis([0 100 -40 0])
xlabel('Neff*p')
ylabel('log10(PDF(Neff*p)/Neff)')
title('Normalized Rician--Neff*SNR')
bold_fig