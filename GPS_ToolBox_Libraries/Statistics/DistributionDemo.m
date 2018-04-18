%Distribution Demo
x=randn([1,10000]);
x_mean=mean(x);
x_std =std(x);
fprintf('mean=%f8.4 std=%f8.4 \n',x_mean,x_std)
[PDF,CDF,yN]=Distribution(x,100);
dy=diff(yN(1:2));

figure
subplot(2,1,1)
plot(yN,CDF,'b')
hold on
plot(yN,(erf(yN/sqrt(2))+1)/2,'r')
ylabel('CDF')
grid on
subplot(2,1,2)
plot(yN,PDF/dy,'b')
hold on
plot(yN,exp(-yN.^2/2)/sqrt(2*pi),'r')
ylabel('PDF')
grid on 
