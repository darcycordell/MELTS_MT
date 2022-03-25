function [results] = MELTS_MT_Rho(H2O_sys,T,P,varargin)
% MELTS_MT_Rho script computes the bulk resistivity as a function of outputs
% from a rhyolite parameterization of MELTS thermodynamic modelling.
%
% Usage: MELTS_MT_Rho(H2O_sys,T,P,varargin)
%
% Required Inputs: 
%         H2O_sys = the total system water content in wt%
%         T = the temperature in Celsius
%         P = the pressure in MPa
%
% Options:
%
% MELTS_MT_Rho(H2O_sys,T,P,c)
%
%   c = melt connectivity parameter for Modified Archie's Law (c > 1.0)
%       (default = 1.15)
%
%
% MELTS_MT_Rho(H2O_sys,T,P,[],rhoh)
%
%   rhoh = set the fixed crystal resistivity
%           (default sets rhoh based on temperature from Hashim et al. (2013))
%           Note that rhoh = 0 (or []) sets default
%   
% MELTS_MT_Rho(H2O_sys,T,P,[],[],volflag,concentration,cf)
%
%   volflag = calculate bulk resistivity in the saturated regime with
%           volatiles present (true) or not (false) (default = false)
%
%   concentration = NaCl-equivalent in wt% for volatile phase
%           (default = 1 wt%)
%
%   cf = connectivity parameter for the volatile phase using 3-phase MAL
%           (default = 3.0; poorly connected)
%
%
% (C) 2022 Darcy Cordell, Georgia Institute of Technology. Please contact
% dcordell6@gatech.edu for more information. See license.txt for more
% information.
%
%
%
%

%Set defaults
c = 1.15;
rhoh = 0;
volflag = false;
concentration = 1;
cf = 3;

if nargin < 3
    error('This function requires H2O_sys, T, and P inputs. See documentation')
end


if nargin > 3
    if ~isempty(varargin{1})
        c = varargin{1};
    end
end

if nargin > 4
    if ~isempty(varargin{2})
        rhoh = varargin{2};
    end 
end

if nargin > 5
    if nargin < 8
        errortext = ['When calculating bulk resistivity in the saturated regime, you must specify\n' ...
             'the true/false flag, the concentration, and the connectivity\n' ...
             'See documentation for more details'];
        error(errortext)
    end
    
    volflag = varargin{3};
    concentration = varargin{4};
    cf = varargin{5};
end
    

%Get all the outputs from MELTS using the MELTS parameterization
[eps_h,eps_m,eps_f,H2O_m] = volfrac(P,T,H2O_sys);

if eps_f > 0 && ~volflag
    disp('Warning: Melt is saturated but you are ignoring the volatiles!')
    rhom = NaN; rhof = NaN; rhob = NaN;  
else
    
    %Calculate melt resistivity from Guo et al. (2016) using MELTS outputs
    rhom = 1./melt_rho(P,T,H2O_m,1);

    %Volatile resistivity computed using Sinmyo and Keppler (2017)
    rhof = 1./Keppler(P,T,concentration,1);

    %If not specified, then calculate crystal resistivity usign Hashim et al. (2013)
    if rhoh == 0 || isempty(rhoh)
        rhoh = exp(-log(2.9)+8581.8/(T+273));
    end

    %Calculate bulk resistivity using three-phase MAL
    rhob = 1./MAL3(1/rhoh,1/rhom,1/rhof,eps_h,eps_m,eps_f,c,cf);

end

results.meltfrac = eps_m;
results.volfrac = eps_f;
results.xfrac = eps_h;
results.meltH2O = H2O_m;
results.rhoh = rhoh;
results.rhom = rhom;
results.rhof = rhof;
results.rhob = rhob;




% 
% %Get the maximum and minimum values of melt fraction, temperature,
% %dissolved water content, etc. along the fixed bulk resistivity contour
% for rhocount = 1:length(rhob_fixed)
%     v = [-1 rhob_fixed(rhocount) 10000];
%     if rhocount == 1
%         lt = '-';
%     else
%         lt = '--';
%     end
%     figure(2);
%     A = contour(H2O_sys,T,(rhob),v,[lt,'w'],'ShowText','off','LineWidth',3);
% 
%     A(:,A(1,:)==rhob_fixed(rhocount))=[];
% 
%     finalx = []; finalmelt = []; finalvolatile = []; finalwater = []; finaltemp = [];
%     for i = 1:length(A(1,:))
% 
%         [finalx(i,rhocount),finalmelt(i,rhocount),finalvolatile(i,rhocount),finalwater(i,rhocount)] = volfrac(P,A(2,i),A(1,i));
%         finaltemp(i,rhocount) = A(2,i);
% 
%     end
%     
%     if rhocount==1
%         figure(1); subplot(1,3,3);plot(finalwater(finalvolatile<=0,rhocount),finaltemp(finalvolatile<=0,rhocount),'-w','LineWidth',3)
%     end
%     
%     if case_study == 1
%     
%         indtemp = intersect(find(finaltemp(:,rhocount)>T1), find(finaltemp(:,rhocount)<T2));
%         indwater = intersect(find(finalwater(:,rhocount)>w1), find(finalwater(:,rhocount)<w2));
%         indpoly = intersect(indtemp,indwater);
% 
%         [maxmelt(rhocount) indmaxmelt] = max(finalmelt(indpoly,rhocount));
%         maxmeltT(rhocount) = finaltemp(indpoly(indmaxmelt),rhocount);
%         maxmeltw(rhocount) = finalwater(indpoly(indmaxmelt),rhocount);
% 
%         [minmelt(rhocount) indminmelt] = min(finalmelt(indpoly,rhocount));
%         minmeltT(rhocount) = finaltemp(indpoly(indminmelt),rhocount);
%         minmeltw(rhocount) = finalwater(indpoly(indminmelt),rhocount);
% 
%     else
%         
%         indpoly = 1:size(A,2);
%         [maxmelt(rhocount) indmaxmelt] = max(finalmelt(indpoly,rhocount));
%         maxmeltT(rhocount) = finaltemp(indpoly(indmaxmelt),rhocount);
%         maxmeltw(rhocount) = finalwater(indpoly(indmaxmelt),rhocount);
%         maxmeltv(rhocount) = finalvolatile(indpoly(indmaxmelt),rhocount);
% 
%         [minmelt(rhocount) indminmelt] = min(finalmelt(indpoly,rhocount));
%         minmeltT(rhocount) = finaltemp(indpoly(indminmelt),rhocount);
%         minmeltw(rhocount) = finalwater(indpoly(indminmelt),rhocount);
%         minmeltv(rhocount) = finalvolatile(indpoly(indminmelt),rhocount);
%     end
%     
% end
% 
% 
% 
% 
% 
% 
