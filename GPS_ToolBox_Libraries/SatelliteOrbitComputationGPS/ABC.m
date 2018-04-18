function  [A,B,C,Sphik,Cphik,tanthk]=ABC(theta_k,phi_k,theta_b,phi_b,gam_b,a,b)
%USAGE:   [A,B,C,Sphik,Cphik,tanthk]=ABC(theta_k,phi_k,theta_b,phi_b,gam_b,a,b)   
%    Compute A,B,C parameters
%
%    Principal propagation direction
%    k=[cos(theta_k),sin(theta_k)cos(phi_k),sin(theta_k)sin(phi_k)]
%
%    Principal axis of elongation
%    s=[cos(theta_b),sin(theta_b)cos(phi_b),sin(theta_b)sin(phi_b)]
%    a=principal axis elongation
%    b=transverse axis elongation
%
    if ~a>=1 || ~b>=1
        error('a,b range')
    end
    Sthetab=sin(theta_b); Cthetab=cos(theta_b);
    Sphib=sin(phi_b);     Cphib=cos(phi_b);
    
    Sgamb=sin(gam_b); Cgamb=cos(gam_b);
    
	c11  = Cgamb.*Cthetab;
	c12  = -(Cgamb.*Sthetab.*Cphib+Sgamb.*Sphib);
	c13  = -(Cgamb.*Sthetab.*Sphib-Sgamb.*Cphib);
	c21  = Sthetab;
	c22  = Cthetab.*Cphib;
	c23  = Cthetab.*Sphib;
	c31  = -Sgamb.*Cthetab;
	c32  = (Sgamb.*Sthetab.*Cphib-Cgamb.*Sphib);
	c33  = (Sgamb.*Sthetab.*Sphib+Cgamb.*Cphib);
    
	Chat11  = c11.^2+(c21*a).^2+(c31*b).^2;
	Chat22  = c12.^2+(c22*a).^2+(c32*b).^2;
	Chat33  = c13.^2+(c23*a).^2+(c33*b).^2;
	Chat12  = c11.*c12+a^2*c21.*c22+b^2*c31.*c32;
	Chat13  = c11.*c13+a^2*c21.*c23+b^2*c31.*c33;
	Chat23  = c12.*c13+a^2*c22.*c23+b^2*c32.*c33;
    
    tanthk=tan(theta_k); Sphik=sin(phi_k); Cphik=cos(phi_k);
    A  = Chat22+tanthk.^2.*Cphik.^2.*Chat11+2*tanthk.*Cphik.*Chat12;
	B  = 2*(Chat23+tanthk.^2.*Sphik.*Cphik.*Chat11...
            +tanthk.*(Chat12.*Sphik+Chat13.*Cphik));
	C  = Chat33+tanthk.^2.*Sphik.^2.*Chat11+2*tanthk.*Sphik.*Chat13;
return