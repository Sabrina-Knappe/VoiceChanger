function [iftloc, iftphase, iftval]= interpolatedvalues(r, phi, N, zp, ftloc, ftval)
%function [iftloc, iftphase, iftval]= interpolatedvalues . . .
% (r, phi, N, zp, ftloc, ftval)
%
%==> computation of the interpolated values
% and phase (linear interpolation)
%
% data:
% iftloc: interpolated location (bin)
% iftval: interpolated magnitude
% iftphase: interpolated phase
% f tloc : peak locations (bin)
% f tval : peak magnitudes
% r: magnitude of the FFT
% phi : phase of the FFT
% N: size of the FFT
% zp : zero-padding multiplicative coefficient
%-- calculate interpolated peak position in bins (--i-f-t-l-o c)
leftftval = r((ftloc-1).*((ftloc-1)>0)+((ftloc-1)<=0).*1);
disp(leftftval);
rightftval= r((ftloc+1).*((ftloc+1)<N/2)+((ftloc+1)>=N/2).*(N/2));
leftftval= 20*log10(leftftval);
rightftval= 20*log10(rightftval);
ftval= 20*log10(ftval);
iftloc= ftloc+ .5*(leftftval-rightftval)./(leftftval-2*ftval+rightftval);

%--- interpolated ftloc --------------------------------
iftloc= (iftloc>=1).*iftloc+(iftloc<1).*1;
iftloc= (iftloc>N/2+1).*(zp/2+1)+(iftloc<=N/2+1).*iftloc;

%--- calculate interpolated phase (iphase)
disp(iftloc)
leftftphase= phi(floor(iftloc));
rightftphase= phi(floor(iftloc)+1);
intpfactor= iftloc-ftloc;
intpfactor= (intpfactor>0).*intpfactor+(intpfactor<0).*(1+intpfactor);
diffphase= unwrap2pi(rightftphase-leftftphase);
iftphase= leftftphase+intpfactor.*diffphase;

%--- calculate interpolate amplitude (iftval) ----------------
iftval= ftval- .25*(leftftval-rightftval).*(iftloc-ftloc);