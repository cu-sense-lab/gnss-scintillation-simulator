function [intensity_d,intensity_bar]=Intensity_Detrend(intensity,Jmin)
%
[nfreqs,~]=size(intensity);
intensity_d=zeros(size(intensity));
intensity_bar=zeros(size(intensity));
for nfreq=1:nfreqs
    intensity_bar(nfreq,:)=WaveletFilter(intensity(nfreq,:),Jmin);
    intensity_d(nfreq,:)=intensity(nfreq,:)./intensity_bar(nfreq,:);
end
return