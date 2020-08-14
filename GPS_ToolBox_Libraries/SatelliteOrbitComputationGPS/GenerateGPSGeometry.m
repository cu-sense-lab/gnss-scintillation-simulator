function        satGEOM_struct=GenerateGPSGeometry(tsec_data,GPSTime,UT_time,eph,origin_llh,StationID,h_intercept,varargin)
%USAGE:   satGEOM_struct=GenerateGPSGeometry(tsec_data,UT_time,eph,origin_llh,StationID,h_intercept,Drift)
%
%INPUTS:
%       tsec_data  = UT in sec 
%       GPSTime  = GPS Time in [week,time] format  
%       UT_time    =  6 element date time for IGRF 
%       eph             = RINEX navagation format ephemeris file
%       origin_llh    = station  [latitude (rad); longitude( rad); height (m)]
%       stationID    = Station ID
%       Drift            = [downward,eastward,southward] drift mps  
%
%OUTPUTS:
%*************ecf Coordinates**********************************************************
%xsat_ecf, vsat_ecf = satellite state vector in ecf coordinates from spg4         (3XN)
%************GPS Coordinates***********************************************************
%sat_llh            = satellite geodetic coordinates                                                   (3XN)
%************Station TCS coordinates***************************************************
%sat_tcs,  vsat_tcs = satellite state vector in receiver tcs system at origin_llh (3XN)
%sat_rng,  sat_rdot = satellite rang & range rate (=> sat_tcs)                            (3XN)
%sat_elev, sat_phi  = satellite elevation & true bearing                                       (NX1)
%*************Propagation Reference Coordinates at penetration point*******************
%xyzp               = propagation coordinates with origin at h_intercept                 (3xN)
%thetap, phip       = polar angles wrt x  (phip  cw from y-axis                              (NX1)
%                                          theta cw from x-axis
%rngp               = range from receiver to intercept point                                      (Nx1)
%uk_xyzp         = unit vector pointing along propagation direction                    (3XN)
%s                     = unit magnetic field vector xp,yp,zp system                              (3xN)
%thetaB,psiB   = polar angles wrt xp                                                                    (Nx1)
%vp                 = penetration point velocity <= satellite motion 
%vk                 = apparent velocity in measurement plane
%sat_utsec          = time (sec)
%
%Libraries:   GPS_Transforms IGRF_Compston 

%Author:
%Charles Rino
%Rino Consulting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%VERSION: February 13, 2016%%%%%%%%%%%%%%%%%%%%

dtr=pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
    Drift=[0;0;0];
else
    Drift=varargin{1};
end

nsamps=length(tsec_data);
xsat_ecf=zeros(3,nsamps);
vsat_ecf=zeros(3,nsamps);
ns=size(GPSTime);
if ns(1)==2
    GPSTime=GPSTime(1,:);
else
    error('GPSTime format error ')
end
for nsamp=1:nsamps            %Calculate ecf position & velocity    
    [ xsat_ecf(:,nsamp), vsat_ecf(:,nsamp)]=satposvel(GPSTime(1,nsamp),eph);  % changed by Joy
end
sat_llh=ecf2llhT(xsat_ecf);                %ECF to geodetic (llh)  
sat_tcs=llh2tcsT(sat_llh,origin_llh);  %llh to tcs at origin_llh
sat_elev=atan2(sat_tcs(3,:),sqrt(sat_tcs(1,:).^2+sat_tcs(2,:).^2));

%%%%%%%%%%Truncate to visible segment containing UTC_time%%%%%%%%%%%%%%%%%% 
nVIS=find(sat_tcs(3,:)>=0);
%%%%%%%%=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Satellite range elevation azimuth (from North)
sat_rnge=sqrt(sat_tcs(1,:).^2+sat_tcs(2,:).^2+sat_tcs(3,:).^2);
sat_phi =atan2(sat_tcs(1,:),sat_tcs(2,:));  

%Satellite velocity & range rate
D=Rotate_ecf2tcs(origin_llh);
vsat_tcs=D*vsat_ecf;
usat_tcs=sat_tcs./repmat(sat_rnge,3,1);
sat_rdot=sum(vsat_tcs.*usat_tcs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('Calculating %6.2f km penetration points \n',h_intercept/1000)
%Intercept point geometry
satp_llh   = findIntercept(h_intercept,usat_tcs,sat_rnge,origin_llh);
satp_tcs = llh2tcsT(satp_llh,origin_llh);

%Convert satellite TCS to xyzp
xyzp_tcs=llh2tcsT(origin_llh,satp_llh);
xyzp(1,:)=-xyzp_tcs(3,:);  %Component 1 is -z tcs (Downward)  
xyzp(2,:)= xyzp_tcs(1,:);  %Component 2 is  x tcs (Eastward)
xyzp(3,:)=-xyzp_tcs(2,:);  %Component 3 is -y tcs (Southward)
rngp  =sqrt(xyzp(1,:).^2+xyzp(2,:).^2+xyzp(3,:).^2);
uk_xyzp  = xyzp./repmat(rngp,3,1);
%Propagation angles
thetap=atan2(sqrt(xyzp(2,:).^2+xyzp(3,:).^2),xyzp(1,:));
phip  =atan2(xyzp(3,:),xyzp(2,:));

vSF=rngp./sat_rnge;  %Velocity scale vactor at penetration point
for nsegs=1:nsamps
    D=Rotate_ecf2tcs(satp_llh(:,nsegs));
    vp(1,:)=-vsat_tcs(3,:).*vSF;
    vp(2,:)= vsat_tcs(1,:).*vSF;
    vp(3,:)=-vsat_tcs(2,:).*vSF;
end

%fprintf('Calculating Magnetic field at %6.2f km penetration points \n',...
%         h_intercept/1000)
%Magnetic field at penetration points 
%Inputs to igrf are deg/km!!
time=datenum(UT_time(1),UT_time(2),UT_time(3));
xyzt=zeros(3,nsamps);
[xyzt(1,:), xyzt(2,:), xyzt(3,:)] =igrf(time,...
               satp_llh(1,:)/dtr, satp_llh(2,:)/dtr, satp_llh(3,:)/1000);  

%Unit vector along B in xyzp system
Bmag=sqrt(xyzt(1,:).^2+xyzt(2,:).^2+xyzt(3,:).^2);
s(1,:)=  xyzt(3,:)./Bmag;    %xyzt(3,:) Dowward
s(2,:)=  xyzt(2,:)./Bmag;    %xyzt(2,:) East
s(3,:)= -xyzt(1,:)./Bmag;    %xyzt(1,:) North
clear xyzt

%s(1,:)=uk_xyzp(1,:);
%s(2,:)=uk_xyzp(2,:);
%s(3,:)=uk_xyzp(3,:);

thetaB=atan2(sqrt(s(2,:).^2+s(3,:).^2),s(1,:));
phiB  =atan2(s(3,:),s(2,:));
cosBP=abs(dot(uk_xyzp,s));

%Penetration point velocity= negative of Equation (4.11) "Theory of Scintillation"
vkyz=[Drift(2)-vp(2,:)+tan(thetap).*cos(phip).*(-Drift(1)+vp(1,:));...
      Drift(3)-vp(3,:)+tan(thetap).*sin(phip).*(-Drift(1)+vp(1,:))];
gam_b=0; a=50; b=1;
thetaB=thetaB+pi/2;    %Change magnetic angle ref to horizontal.  CLR 2/13/2016
    [A,B,C]=ABC(thetap,phip,thetaB,phiB,gam_b,a,b);

veff=sqrt((C.*vkyz(1,:).^2-B.*vkyz(1,:).*vkyz(2,:)+A.*vkyz(2,:).^2)./(A.*C-B.^2/4));
tsec_data=tsec_data(:)';
satGEOM_struct=struct('tsec_data',tsec_data,'eph',eph,...
    'UT_time',UT_time,'origin_llh',origin_llh','StationID',StationID,...
    'h_intercept',h_intercept','a',a,'b',b,'gam_b',gam_b,'Drift',Drift,...
    'xsat_ecf',xsat_ecf,'vsat_tcs',vsat_tcs,'sat_llh',sat_llh, 'sat_tcs',sat_tcs,...
    'sat_rnge',sat_rnge, 'sat_rdot',sat_rdot, 'sat_elev',sat_elev, 'sat_phi',sat_phi,...
    'satp_tcs',satp_tcs, 'satp_llh',satp_llh, 'xyzp',xyzp,'uk_xyzp',uk_xyzp,...
    'vkyz',vkyz,'rngp',rngp, 'thetap',thetap, 'phip',phip,'veff',veff,'A',A,'B',B,'C',C,...
    's',s,'thetaB',thetaB,'phiB',phiB,'cosBP',cosBP,'nVIS',nVIS);
return