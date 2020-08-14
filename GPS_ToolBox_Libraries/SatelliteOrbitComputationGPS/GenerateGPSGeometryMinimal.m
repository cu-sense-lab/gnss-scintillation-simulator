function        [sat_rnge,sat_rdot,sat_elev,sat_phi,sat_llh,sat_tcs]=GenerateGPSGeometryMinimal(tsec_data,GPSTime,eph,origin_llh)
%USAGE:   [sat_rnge,sat_rdot,sat_elev,sat_phi,sat_llh,sat_tcs]=GenerateGPSGeometryMinimal(tsec_data,GPSTime,eph,origin_llh)
%
%INPUTS:
%       tsec_data  = GPS time in seconds 
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
%Libraries:   GPS_Transforms 

%Author:
%Charles Rino
%Rino Consulting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%VERSION: February 13, 2016%%%%%%%%%%%%%%%%%%%%
nsamps=length(tsec_data);
xsat_ecf=zeros(3,nsamps);
vsat_ecf=zeros(3,nsamps);
for nsamp=1:nsamps            %Calculate ecf position & velocity    
    [ xsat_ecf(:,nsamp), vsat_ecf(:,nsamp)]=satposvel(GPSTime(2,nsamp),eph);  % changed by Joy
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
return