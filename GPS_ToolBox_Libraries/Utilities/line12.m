function y=line12(x,xy1,xy2)
%
%  Non vertical line versus x passing through xy1 and xy2
%
if xy1(1)==xy2(1)
    error('use line21');
else
    y=(xy2(2)-xy1(2))*(x-xy1(1))/(xy2(1)-xy1(1))+xy1(2);   
end
return