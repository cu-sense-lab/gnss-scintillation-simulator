function DisplayGeometry(satGEOM,PlotID,varargin)
%USAGE:  DisplayGeometry(satGEOM,PlotID,Plot_directory)

if ~isempty(varargin)
    SAVE_PLOTS=1;
    fig_dir=varargin{1};
else
    SAVE_PLOTS=0;
end
dtr=pi/180;

t_secG=satGEOM.tsec_data;
origin_llh=satGEOM.origin_llh;

UTsec=t_secG;
sat_phi  =satGEOM.sat_phi/dtr;
nNeg=find(sat_phi<0);

%Summarize Propagation Geometry
thetap  =satGEOM.thetap;
phip     =satGEOM.phip;
nG180=find(phip/dtr>180);
phip(nG180)=phip(nG180)-2*pi;
cosBP=satGEOM.cosBP;

figure
subplot(3,1,1)
plot(UTsec/60,thetap/dtr,'r')
grid on
ylabel('\theta_p-deg')
title(PlotID)
subplot(3,1,2)
plot(UTsec/60,phip/dtr,'r')
grid on
ylabel('\phi_p-deg')
subplot(3,1,3)
plot(UTsec/60,cosBP,'r')
grid on
ylabel('cosBP')
xlabel('t-min')
bold_fig
if SAVE_PLOTS==1
    saveas(gcf,[fig_dir,'\satp_angBP'],'jpg')
end

%%%%%%%%%%3D views%%%%%%%%%%%%%%

%Pierce point in llh coordinates
satp_llh=satGEOM.satp_llh;
nG180=find(satp_llh(2,:)>180*dtr);
satp_llh(2,nG180)=satp_llh(2,nG180)-360*dtr;


figure
plot3(satp_llh(2,:)/dtr,satp_llh(1,:)/dtr,satp_llh(3,:)/1000,'r')
hold on
plot3(origin_llh(2)/dtr,origin_llh(1)/dtr,origin_llh(3)/1000,'mp')
grid on
xlabel('longitude-deg (x)')
ylabel('latitude-deg (y)')
zlabel('height-km (z)')
title([PlotID,' Pierce Point'])
bold_fig
if SAVE_PLOTS==1
    saveas(gcf,[fig_dir,'\satp_llh'],'jpg')
end
    
%Pierce point in propagation coordinate system (PCS)
satp_tcs=satGEOM.satp_tcs;
uk_xyzp =satGEOM.uk_xyzp;
rngp      =satGEOM.rngp;
%Direction of magnetic field 
s=satGEOM.s;  %Magnetic field vector in PCS system

nsamps=length(UTsec);
%Decimate for plot
nskip=ceil(nsamps/10);
ndec=(1:nskip:nsamps);
ndec=[ndec,nsamps];
figure
plot3(0,0,0,'r')
hold on
plot3(0,0,0,'k')
hold on
plot3(0,0,0,'b')
plot3(satp_tcs(3,:)/1000,satp_tcs(2,:)/1000,satp_tcs(1,:)/1000,'r')
hold on
eps=1.0e-1;
for nn=ndec
    line([satp_tcs(3,nn),satp_tcs(3,nn)-eps*rngp(nn)*uk_xyzp(1,nn)]/1000,...            %xp=>-z
        [satp_tcs(2,nn),satp_tcs(2,nn)-eps*rngp(nn)*uk_xyzp(3,nn)]/1000,...               %zp=>-y
        [satp_tcs(1,nn),satp_tcs(1,nn)+eps*rngp(nn)*uk_xyzp(2,nn)]/1000,'color','k')  %yp=>x
    hold on
    line([satp_tcs(3,nn),satp_tcs(3,nn)+eps*rngp(nn)*s(3,nn)]/1000,...
        [satp_tcs(2,nn),satp_tcs(2,nn)+eps*rngp(nn)*s(2,nn)]/1000,...
        [satp_tcs(1,nn),satp_tcs(1,nn)+eps*rngp(nn)*s(1,nn)]/1000,'color','b')
    grid on
end
xlabel('Eastward (x)')
ylabel('Northward (y)')
zlabel('Upward (z)')
title([PlotID,' Pierce Point'])
legend('Origin','To Receiver','B-field')
bold_fig
if SAVE_PLOTS==1
    saveas(gcf,[fig_dir,'\satp_tcs'],'jpg')
end

A=satGEOM.A;
B=satGEOM.B;
C=satGEOM.C;
D=A.*C-B.^2/4;
a=satGEOM.a;
b=satGEOM.b;
GFAC=a*b*sec(thetap)./sqrt(D);
veff=satGEOM.veff;

figure
subplot(3,1,1)
plot(t_secG/60,cosBP,'r')
grid on
ylabel('cosBP')
title(PlotID)
subplot(3,1,2)
plot(t_secG/60,GFAC,'r')
grid on
ylabel('GFAC')
subplot(3,1,3)
plot(t_secG/60,veff,'r')
grid on
ylabel('veff-mps')
xlabel('t-min')
bold_fig
if SAVE_PLOTS==1
    saveas(gcf,[fig_dir,'\prop_BP-GFAC-veff'],'jpg')
end

if 0
figure
scale=1;
for nn=ndec
    Z=satp_tcs(2,nn);
    Y=satp_tcs(3,nn);
     DisplayAnisotropy(A(nn),B(nn),C(nn),Y/1000,Z/1000,scale)
     hold on
     plot(Y/1000,Z/1000,'rp')
end
xlabel('Eastward')
ylabel('Northward')
grid on
title([PlotID,' Anisotropy'])
axis equal
bold_fig
if SAVE_PLOTS==1
    saveas(gcf,[fig_dir,'\Anisotropy'],'jpg')
end
end
