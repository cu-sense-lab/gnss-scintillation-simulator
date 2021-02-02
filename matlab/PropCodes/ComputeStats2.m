function       [SDF,q,ACF,CDF,PDF,yU]=ComputeStats2(data,data_span,nseg,nstep,nfft,units)
% USAGE: [SDF,q,ACF,CDF,PDF,yU]=ComputeStats2(data,data_span,nseg,nstep,nfft,units)
%,,
%Intensity Statistical Analysis
%INPUT:
%                    data                     = Intensity data uniform spatial sampling (1,nsamps)
%                    data_span span = time or distance                                       (1,nsamps)
%                    nfft    =fft size >   = nste
%
%OUTPUT:  
%                    SDF=spectral intensity function of data/<data>                   
%                               Normalization: sim(SDF)=sum(abs(data).^2);
%                         q= frequency 
%                    ACF=Autocorrelation function
%                    CDF=cumulative distribution
%                    yU    =CDF ordinate
%

nfft_span=-nfft/2:nfft/2-1;                 %Natural order fft index span                                            
nrecs=length(nseg);                         %Number of N-sample records 

SDF=zeros(nrecs,nfft);
q=zeros(nrecs,nfft);
ACF=SDF;
CDF=zeros(nrecs,nstep-1);
PDF=CDF;
yU=PDF;
ds=diff(data_span(1:2));;
nsamps=length(data_span);
span=nsamps*ds;
if units==1
    dq=1/span;
else
    dq=2*pi/span;
end

for nrec=1:nrecs
     n1=(nrec-1)*nstep+1; n2=min(n1+nstep-1,nsamps);
    q(nrec,:)=nfft_span'*2*pi/span;
    SDFnrec=abs(fft(data(n1:n2),nfft)).^2/nfft;
    SDF(nrec,:)=SDFnrec';
    ACF(nrec,:)=ifft(SDFnrec)';
    [PDFtemp,CDFtemp,yUtemp]=Distribution(data(n1:n2),20);
    CDF(nrec,1:length(CDFtemp))=CDFtemp;
    PDF(nrec,1:length(CDFtemp))=PDFtemp;
    yU(nrec,1:length(CDFtemp))=yUtemp;
end

return