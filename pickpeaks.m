function [loc, val]=pickpeaks(spectrum, nPeaks, minspace)
%function [loc, val]= pickpeaks(spectrum, nPeaks, minspace)
%
%==>peaking the nPeaks highest peaks in the given spectrum from the
%greatest to the lowest
%data:
%   loc: bin number of peaks (if loc(i)==0, no peak detected)\
%   val: amplitude of the given spectrum
%   spectrum: spectrum (abs(fft(signal))
%   nPicks: number of peaks to pick
%   minspace: minimum of space between two peaks

[r, c]= size(spectrum)
rmin=min(spectrum)-1

%--- find a peak, zero out the data around the peak, and repeat
val= ones(nPeaks, c)*NaN;
loc= zeros(nPeaks, c);

for k=1:c  %--- find all local peaks
    difference= diff([rmin; spectrum(:,k); rmin]); %derivate
    iloc= find(difference(1:r)>=0 & difference(2:r+1) <=0);
        %peak locations
    ival= spectrum(iloc, k) %peak values
    
    for p=1:nPeaks
        [val(p,k),l]= max(ival);    %find current maximum
        loc(p,k)= iloc(l);          %save value and location
        ind= find(abs(iloc(l)-iloc)>minspace);
            %find peaks which are far away
        if(isempty(ind))
            break         %no more local peaks to pick
        end
        ival= ival(ind);    %shrink peaks value and location array
        iloc= iloc(ind);
    end
end
