function fracMom=generate_fracMom(I)
%  Generate Fractional Moments
fracMom=zeros(1,5);
fracMom(1)=mean(I(:));
for m=2:5;
        fracMom(m)=mean(I(:).^m);
        fracMom(m)=fracMom(m)/fracMom(1)^m;
end
return