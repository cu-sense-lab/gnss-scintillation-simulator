function   [Path_struct,UT_StartTime]=Generate_Path_struct(t_seg,satGEOM_struct)

tsec_GEOM=satGEOM_struct.tsec_data;
nsamps=length(tsec_GEOM);
UT_StartTime=satGEOM_struct.UT_StartTime;
%fprintf('year=%4i  %2i/%2i  min=%s \n',UT_StartTime(1:3),num2str(t_seg(1)/60))
satp_tcs=satGEOM_struct.satp_tcs;
s=satGEOM_struct.s;
nsegs=length(t_seg);
nskip=floor(nsamps/nsegs);
nseg_locs=[nskip/2:nskip:nsamps];
ntG=zeros(size(t_seg));
for n=1:nsegs
    [~,ntG(n)]=min(abs(t_seg(n)-tsec_GEOM));
end
dr=diff(satp_tcs')';
path=zeros(1,nsamps);
for n=1:nsamps-1
    path(n+1)=path(n)+[sqrt(dr(1,n).^2+dr(2,n).^2+dr(3,n).^2)];
end
path_seg=zeros(1,nsegs);
for n=1:nsegs-1
    path_seg(n)=path(ntG(n+1))-path(ntG(n));
end
path_seg(nsegs)=path(ntG(nsegs)+10)-path(ntG(nsegs-1));  %Ad Hoc fix

Path_struct.path=path;
Path_struct.tsec_GEOM=tsec_GEOM;
Path_struct.ntG=ntG;
Path_struct.path_seg=path_seg;
return
    