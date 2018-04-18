function  OMEGA=tcs2ecfR(origin)
%  Generate tcs2ecf rotation matrix
%  v*OMEGA'=tcs2ecf_vel(v,origin);
%
YAW_TYPE = 1;
PITCH_TYPE = 2;
ROLL_TYPE = 3;

DC1 = Get_DirCos_ForwardT (pi/2, ROLL_TYPE);
DC2 = Get_DirCos_ForwardT (pi/2, PITCH_TYPE);
DC21=DC2*DC1;
DC3 = Get_DirCos_ForwardT (-origin(1,:), ROLL_TYPE);
DC4 = Get_DirCos_ForwardT (-origin(2,:), YAW_TYPE);

sz=size(origin);
if sz(2)==1
    DC = DC3*DC21;
    DCINV = DC^(-1);
    OMEGA=DC4*DCINV;
else
    OMEGA=zeros(3,3,sz(2));
    for n=1:sz(2)
        DC=DC3(:,:,n)*DC21;
        DCINV=DC^(-1);
        OMEGA(:,:,n)=DC4(:,:,n)*DCINV;
    end
end
return