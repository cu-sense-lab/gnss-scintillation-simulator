function  [Type,ID]=SpecFit(xxScale,yyScale)
%USAGE    [Type,ID]=SpecFit(xxScale,yyScale)
%Two-component log-linear fit, classification, and parameter extraction
%INPUT:
%         xxScale=log10(q)            log10 units!  
%         yyScale=10log10(SDF)        dB    units!
%OUTPUT:
%         log-linear fit=>
%         yyScale ~ cs1*q^ps1 for q<q(ns)
%                   cs2*q^ps2 for q>q(ns)
%         Type.Cs1=cs1_dB;
%         Type.ps1=ps1;
%         Type.Cs2=cs2__dB;
%         Type.ps2=ps2;
%         Type.ns =ns;
%         Type.err=err;
%
%Written by Charles Rino
%Rino Consulting
%November 2012
%crino@comcast.com
%http://chuckrino.com/wordpress
%

%Generate 2-component fits with varying small-large partition
nb_min=4; ns_max=length(xxScale);
%Arrays truncated to search )
err=zeros(2,ns_max-2*nb_min+2);
if isempty(err)
    Type=[];
    ID=[];
    return
end
cs1_dB=err(1,:); cs2_dB=err(1,:); ps1=err(1,:); ps2=err(1,:);
for ns=nb_min:ns_max-nb_min+1;
    P1=polyfit(xxScale(1:ns),       yyScale(1:ns),1);
    P2=polyfit(xxScale(ns+1:ns_max),yyScale(ns+1:ns_max),1);
    yy1 = polyval(P1,xxScale(1:ns));
    yy2 = polyval(P2,xxScale(ns+1:ns_max));
    err(1,ns-nb_min+1)=sum((yy1-yyScale(1:ns)).^2)/ns;
    err(2,ns-nb_min+1)=sum((yy2-yyScale(ns+1:ns_max)).^2)/(ns_max-ns);
    %Generate Parameters
    ps1(ns-nb_min+1)=-P1(1)/10;   ps2(ns-nb_min+1)=-P2(1)/10;   %Compensate for dB cs units 
    cs1_dB(ns-nb_min+1)=P1(2); cs2_dB(ns-nb_min+1)=P2(2);       %dB unitsfigure
end
%Parameter extraction
    sig=sum(err(1:2,:));
    [~,ns_min]=min(sig);
    Type.Cs1=cs1_dB(ns_min);
    Type.ps1=ps1(ns_min);
    Type.Cs2=cs2_dB(ns_min);
    Type.ps2=ps2(ns_min);
    if isempty(err)
        ID=[];
    elseif ps1(ns_min)> ps2(ns_min)
        ID='Noise';
    else
        ID='Two-Component';
    end
    ns=ns_min+nb_min;     %ns in non truncated array     
    Type.ns =ns;
    Type.sig=sig;
return