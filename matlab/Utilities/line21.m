function x=line21(y,xy1,xy2)
%
% Non horizontal line versus y passing through xy1 and xy2
%
if xy1(2)==xy2(2)
    error('Use line12')
else
    x=(xy2(1)-xy1(1))*(y-xy1(2))/(xy2(2)-xy1(2))+xy1(1);   
end
return