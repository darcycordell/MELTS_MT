function [sig] = MAL3(sigc,sigm,sigg,phic,phim,phig,mm,mg)
% Solve Modified Archie's Law for 3-phases (Glover, 2010)
%
% Usage: [sig] = MAL3(sigc,sigm,sigg,phic,phim,phig,mm,mg)
%
% Inputs:
%       sigc = crystal conductivity
%       sigm = melt conductivity
%       sigg = volatile conductivity
%       phic = crystal fraction
%       phim = melt fraction
%       phig = volatile fraction
%       mm = melt connectivity
%       mg = volatile connectivity
%       
if phic == 0
    phic = 10^-4;
end

Gm = phim.^mm;
Gg = phig.^mg;

mc = log(1-Gg-Gm)./log(phic);

Gc = phic.^mc;

Gchecksum = Gc + Gm + Gg;

sig = sigc.*Gc + sigm.*Gm + sigg.*Gg;

if abs(Gchecksum-1) > 10^-3
   error('G not equal 1. Something wrong.')
end

if imag(1/sig)>10^-3
    error('Imaginary conductivies. Something wrong.')
else
    sig = real(sig);
end
    

%     phichecksum = phic + phim + phig;
%     if abs(phichecksum-1) > 10^-3
%         disp('phi not equal 1')
%     end






end