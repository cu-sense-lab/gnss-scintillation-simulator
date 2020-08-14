%**************************************************************************
%Demo GPS Coordinate Trahsformations
%Chuck Rino
%Rino Consulting
%crino@comcast.net
clear
close all
dtr=pi/180;
%**********************************************************
%Generate a source at Rheight flying at a constant altitude 
nsec=200;
Rlat= 36.7788;                      
Rlon=-75.9573; 
Rheight=4000;

dsource_llh = 4*[.00001; .00002; 0];  
source_llh = zeros(3,nsec);
source_llh(:,1)=[Rlat*dtr; Rlon*dtr; Rheight];
for n=2:nsec
    source_llh(:,n)=source_llh(:,n-1)+dsource_llh;
end

%Generate offset target at same height
target_llh = zeros(3,nsec);
target_llh(:,1)=source_llh(:,nsec);
target_llh(3,1)=0;
dtarget_llh = -4*[.00002; .00001; 0];
for n=2:nsec
    target_llh(:,n)=target_llh(:,n-1)+dtarget_llh;
end
%************************************************************

%Calculate source positions in target_tcs system
source_tcs=llh2tcsT(source_llh,target_llh);
range1=sqrt(source_tcs(1,:).^2+source_tcs(2,:).^2+source_tcs(3,:).^2);
bearing1=atan2(source_tcs(1,:),source_tcs(2,:));

%Calculate target positions in source_tcs system
target_tcs=llh2tcsT(target_llh,source_llh);
range2=sqrt(target_tcs(1,:).^2+target_tcs(2,:).^2+target_tcs(3,:).^2);
bearing2=atan2(target_tcs(1,:),target_tcs(2,:));

figure
subplot(2,1,1)
plot(range1/1000,'b.')
hold on
plot(range2/1000,'ro')
grid on
ylabel('Target range--km')
title('Switch observer & target')
subplot(2,1,2)
plot(bearing1/dtr,'b.')
hold on
plot(bearing2/dtr,'ro')
grid on
ylabel('Target bearing--deg')
title('Bearing changes 180--deg')

source1_XYZ=tcs2ecfT(source_tcs,target_llh);
source1_llh=ecf2llhT(source1_XYZ);
source2_llh=tcs2llhT(source_tcs,target_llh); 

target1_XYZ=tcs2ecfT(target_tcs,source_llh);
target1_llh=ecf2llhT(target1_XYZ);

figure
plot(target_llh(1,:)/dtr,target_llh(2,:)/dtr,'r')
hold on
plot(source_llh(1,:)/dtr,source_llh(2,:)/dtr,'b')
hold on
plot(target1_llh(1,:)/dtr,target1_llh(2,:)/dtr,'r.')
hold on
plot(source1_llh(1,:)/dtr,source1_llh(2,:)/dtr,'b.')
hold on
plot(source2_llh(1,:)/dtr,source2_llh(2,:)/dtr,'k.')
grid on
legend('target','source')
xlabel('latitude--deg')
ylabel('longitude--deg')
title('Latitude and longitude from tcs')