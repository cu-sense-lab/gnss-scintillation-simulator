function [psd,nsg] = pdgm(z,nfft,noff)
%USAGE:  [psd,nsg] = pdgm(z,nfft,noff)
%
%PURPOSE: Computs raw periodogam = 1/N*|fft(z,nfft)|^2
%
%INPUT:     
%     z    = complex data array
%     nfft = size of fft 
%     noff = offset for successive nfft-length segments (<= nfft)
%          = 0 for no overlap
%OUTPUT:
%     psd  = Power spectral density
%     nsg  = number of segments averaged 
%
      nz  = length(z);
      nsg = fix((max(nz,nfft)-nfft+1)/(nfft+noff))+1;
      sz  = size(z);
      if sz(1) > 1 
         sz(1) = nfft;
      else
         sz(2) = nfft;
      end
      psd = zeros(sz);
      for nn=1:nsg
         nfrst = (nn-1)*noff + 1;
         nlast = nfrst + min(nz,nfft )- 1;
         psd = psd + abs(fft(z(nfrst:nlast),nfft)).^2/nfft;
      end
      psd = psd/nsg;
      return         