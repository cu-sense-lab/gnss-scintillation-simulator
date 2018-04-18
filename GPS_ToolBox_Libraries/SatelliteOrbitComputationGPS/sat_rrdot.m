function [sat_rng,sat_rdot]=sat_rrdot(UTsec,eph,UT2eph,origin_llh)
%
nsamps=length(UTsec);
xsat_ecf=zeros(3,nsamps);
vsat_ecf=zeros(3,nsamps);
for nsamp=1:nsamps                 %Calculate ecf position & velocity
    [ xsat_ecf(:,nsamp), vsat_ecf(:,nsamp)]=satposvel(UTsec(nsamp)+UT2eph,eph);
end
sat_llh=ecf2llhT(xsat_ecf);            %ECF to geodetic (llh)  
sat_tcs=llh2tcsT(sat_llh,origin_llh);  %llh to tcs at origin_llh
notVIS=find(sat_tcs(3,:)<0, 1);
if ~isempty(notVIS)
    error('Visibility ')
end
sat_rng=sqrt(sat_tcs(1,:).^2+...
    sat_tcs(2,:).^2+sat_tcs(3,:).^2);

%Satellite range elevation azimuth (from North)
sat_rnge=sqrt(sat_tcs(1,:).^2+sat_tcs(2,:).^2+sat_tcs(3,:).^2);

%Satellite velocity & range rate
D=Rotate_ecf2tcs(origin_llh);
vsat_tcs=D*vsat_ecf;
usat_tcs=sat_tcs./repmat(sat_rnge,3,1);
sat_rdot=sum(vsat_tcs.*usat_tcs);