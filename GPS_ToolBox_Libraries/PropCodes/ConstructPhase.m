function   [phase_int,nint,diffW,Wdiff]=ConstructPhase(psi,varargin)
%USAGE:  [phase,nint]=ConstructPhase(psi)
%PURPOSE: Construct unwrapped phase from complex field with check for unwrapping errors
%INPUT:      psi = samplex complex field
%
%OUTPUT   phase=complex filed phase: psi=|psi|*exp(i*phase)
%                           sampled at nint 
%                int = interpolation required for error-free unwrapping
%

if ~isempty(varargin)
    verbose=1;
else
    verbose=0;
end
nint_max=10;
nint=0; err=10;
while err>1 && nint<nint_max
    nint=nint+1;
    psi_nint=fftInterp(psi,nint);                                                   %Interpolate psi
    [phase_int,diffW,Wdiff]=unwrapPhase(wrapPhase(psi_nint));                  %Construct phase_ing
    phase=phase_int(1:nint:end); 
    %Extract original samples
    if nint==1
        phaseL=phase;
    else
        DP=phaseL-phase;
        err=rms(DP);
        if verbose
           fprintf('nint=%5i  err=%12.4f \n',nint,err)
        end
        phaseL=phase;
    end
end
return