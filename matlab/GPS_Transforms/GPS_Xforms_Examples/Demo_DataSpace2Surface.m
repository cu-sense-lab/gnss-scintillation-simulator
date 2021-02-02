%**************************************************************************
%Program to DEMO DataSpace2Surface Utility
clear
close all
dtr=pi/180;
rtd=180/pi;
%%%%  Radar coordinates  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rlat= 36.7788;                        %36 46.73'  N
Rlon=-75.9573;                        %75 57.44'  W                           %height

iStatic=input('Input 1 for static source, else CR ');
type{1}='Static ';
type{2}='Moving';
if isempty(iStatic)
    iStatic=2;
end
%Generate geodetic coordinates of source flying at constant altitude Rheight 
nsec=500;
Rheight=4000;
source0_llh = [Rlat*dtr; Rlon*dtr; Rheight];
dsource_llh = 4*[.00001; .00002; 0];  
source_llh = zeros(3,nsec);
source_llh(:,1)=source0_llh;
for n=2:nsec
    source_llh(:,n)=source_llh(:,n-1)+dsource_llh;
end

%Generate a moving ground target track
target_llh = zeros(3,nsec);
target_llh(:,1)=source_llh(:,nsec);
target_llh(3,1)=0;
dtarget_llh = -4*[.00002; .00001; 0];
for n=2:nsec
    target_llh(:,n)=target_llh(:,n-1)+dtarget_llh;
end

if iStatic==1
   source_llh=source_llh(:,1);
end

figure
plot3(target_llh(1,:)*rtd,target_llh(2,:)*rtd,target_llh(3,:),'r')
hold on
plot3(target_llh(1,1)*rtd,target_llh(2,1)*rtd,target_llh(3,1),'r>')
hold on
plot3(source_llh(1,:)*rtd,source_llh(2,:)*rtd,source_llh(3,:),'b')
hold on
plot3(source_llh(1,1)*rtd,source_llh(2,1)*rtd,source_llh(3,1),'b>')
grid on
xlabel('Latitude--deg')
ylabel('Longitude--deg')
zlabel('Height--meters')
title([type{iStatic},'Source & Target geodetic system'])

%Locate the target in a tcs system fixed on the aircraft 
target_tcs=llh2tcsT(target_llh,source_llh);     %Target as seen from moving source
%Now calculate the range and bearing to the target 
range=sqrt(target_tcs(1,:).^2+target_tcs(2,:).^2+target_tcs(3,:).^2);
bearing=atan2(target_tcs(1,:),target_tcs(2,:));
bearing=mod(bearing,2*pi);
%Calculate the range to the visible horizon at each bearing
rangeH=GHorizon(bearing,source_llh,0);
rangeH=rangeH';
iOK=find(range<rangeH);
figure
plot(bearing*rtd,range/1000,'go')
hold on
plot(bearing(iOK)*rtd,range(iOK)/1000,'b.')
grid on
xlabel('Bearing--deg')
ylabel('Range--km')
legend('All Ranges','Visible Ranges')
title([type{iStatic},'Radar Data'])

%Now reconstruct target positions from radar reports
[stheta, ctheta]=Source2Ellipsoid(range,bearing,source_llh,0);
v_tcs=[stheta'.*sin(bearing); stheta'.*cos(bearing); ctheta'];
tgt_tcs=repmat(range,3,1).*v_tcs;

figure
plot3(target_tcs(1,:)/1000,target_tcs(2,:)/1000,target_tcs(3,:)/1000,'b')
hold on
plot3(tgt_tcs(1,:)/1000,tgt_tcs(2,:)/1000,tgt_tcs(3,:)/1000,'mo')
grid on
xlabel('x--km')
ylabel('y--km')
zlabel('z--km')
title([type{iStatic},'Source TCS coordinates of target'])
legend('Truth','Reports')
