function    R=RTCSP(psi,phi)
%             R=RTCSP(psi,phi)
%              Rotate from TCS to x-aligned system
%             psi=polar angle from vertical
%             phi=polar angle from yz plane

cpsi=cos(psi); spsi=sin(psi);
cphi=cos(phi); sphi=sin(phi);
R=zeros(3,3,length(psi));
R(1,1,:)=cpsi.*cphi;
R(1,2,:)=cpsi.*sphi;
R(1,3,:)=spsi;
R(2,1,:)=sphi;
R(2,2,:)=-cphi;
R(2,3,:)=zeros(size(cpsi));
R(3,1,:)=spsi.*cphi;
R(3,2,:)=spsi.*sphi;
R(3,3,:)=-cpsi;
R=squeeze(R);
return
  