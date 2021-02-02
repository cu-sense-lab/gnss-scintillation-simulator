%DisplayRadarHorizon
clear
close all
dtr=pi/180;
%%%%   Dam Neck Radar    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rlat= 36.7788;                        %36 46.73'  N
Rlon=-75.9573;                        %75 57.44'  W                       
Rheight=10000;

source0_llh=[Rlat; Rlon; 0]*dtr+[0; 0; +Rheight];

bscan=linspace(0,2*pi,500);
[rangeH, sthetaH, cthetaH]=GHorizon(bscan,source0_llh,0);
v_tcs=[sthetaH'.*cos(bscan); sthetaH'.*sin(bscan); cthetaH'];
figure
plot(rangeH'.*v_tcs(1,:)/1000,rangeH'.*v_tcs(2,:)/1000,'b')
grid on
hold on
plot(0,0,'mp')
xlabel('X-km')
ylabel('y-km')
title(['Range to radar horizon from',num2str(Rheight),'-m height'])
