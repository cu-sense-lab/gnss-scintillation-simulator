function  r=randp(y,F,nsamp)
%  Generate nsamp realizations of rv with cumulative density F 
%    F(y)     =integral of non-negative function with domain ymin<=y<=ymax
%        
%    nsamp =number of realizations
%
%    METHOD: Inverse map uniform variates onto range of CDF
%                 y=g(u)
%             p(y)=p(u)/g'(u)
%             p(u)=1 & u=CDF(y) => p(y)=int(CDF)                
%
%
method=2;
%Eliminate duplicates  
[Fu,iF]=unique(F);
yu=y(iF);
%%%%%%%%%%%%%%%%%%%%%%
u=rand(nsamp,1);
nLow= u<Fu(1);
u(nLow)=Fu(1);
nHigh= u>Fu(end);
u(nHigh)=Fu(end);
r=interp1(Fu,yu,u);
return