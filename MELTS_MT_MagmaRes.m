function [magmaRes] = MELTS_MT_MagmaRes(rhob_fixed,H2O_sys,T,P,varargin)
% MELTS_MT_MagmaRes script estimates the range of magma reservoir conditions
% given a bulk resistivity value using a rhyolite parameterization of MELTS 
% thermodynamic modelling to constrain parameters.
%
% Usage: MELTS_MT_MagmaRes(H2O_sys,T,P,varargin)
%
% Outputs: magmaRes is a cell variable containing matrices of total system
% water, temperature, melt fraction, volatile fraction, and dissolved water
% content, and bulk resistivity at a pressure. The cell index is the
% pressure index. The bulk resistivity will be near rhob_fixed within a
% fixed tolerance.
%
% Required Inputs: 
%         rhob_fixed = The given bulk resistivity value in Ohm m
%         H2O_sys = Range of total system water contents to consider wt%
%                   e.g. [0:0.1:7]
%         T = Range of temperatures to consider in Celsius
%                   e.g. [700:10:1000]
%         P = Range of pressures to consider in MPa
%                   e.g. [100 200]
%
% Options:
%
% MELTS_MT_MagmaRes(rhob_fixed,H2O_sys,T,P,c)
%
%   c = melt connectivity parameter for Modified Archie's Law (c > 1.0)
%       (default = 1.15)
%
%
% MELTS_MT_MagmaRes(rhob_fixed,H2O_sys,T,P,[],rhoh)
%
%   rhoh = set the fixed crystal resistivity
%           (default sets rhoh based on temperature from Hashim et al. (2013))
%           Note that rhoh = 0 (or []) sets default
%   
% MELTS_MT_MagmaRes(rhob_fixed,H2O_sys,T,P,[],[],volflag,concentration,cf)
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
%
% NOTE:
% The code is currently written to brute force the solution by exploring
% the entire parameter space. As a result, the total number of calculations
% required is length(H2O_sys)*length(T)*length(P). If this is >10^4, the
% calculation takes quite a while. This could be made more sophisticated
% using e.g. Bayesian sampling.
%
% (C) 2022 Darcy Cordell, Georgia Institute of Technology. Please contact
% dcordell6@gatech.edu for more information. See license.txt for more
% information.
%
%
%
%

totalCalc = length(H2O_sys)*length(P)*length(T);

if totalCalc > 10^4
    mn = menu('Your input vectors require greater than 10^4 calculations and this may take a while. Proceed?','Yes','No');
    
    if mn ~= 1
        magmaRes = NaN;
        return
    end
end

%Set defaults
c = 1.15;
rhoh = 0;
volflag = false;
concentration = 1;
cf = 3;

if nargin < 4
    error('This function requires rhob_fixed, H2O_sys, T, and P inputs. See documentation')
end


if nargin > 4
    if ~isempty(varargin{1})
        c = varargin{1};
    end
end

if nargin > 5
    if ~isempty(varargin{2})
        rhoh = varargin{2};
    end 
end

if nargin > 6
    if nargin < 9
        errortext = ['When calculating bulk resistivity in the saturated regime, you must specify\n' ...
             'the true/false flag, the concentration, and the connectivity\n' ...
             'See documentation for more details'];
        error(errortext)
    end
    
    volflag = varargin{3};
    concentration = varargin{4};
    cf = varargin{5};
end
    
eps_m = nan(length(T),length(H2O_sys),length(P));
eps_f = nan(size(eps_m)); eps_h = nan(size(eps_m));
H2O_m = nan(size(eps_m));
rhom = nan(size(eps_m)); rhof = nan(size(eps_m));
rhob = nan(size(eps_m));
for i = 1:length(T)
    for j = 1:length(H2O_sys)
        for k = 1:length(P)
            
        [results] = MELTS_MT_Rho(H2O_sys(j),T(i),P(k),c,rhoh,volflag,concentration,cf);
        
        eps_m(i,j,k) = results.meltfrac;
        eps_f(i,j,k) = results.volfrac;
        eps_h(i,j,k) = results.xfrac;
        H2O_m(i,j,k) = results.meltH2O;
        rhom(i,j,k) = results.rhom;
        rhof(i,j,k) = results.rhof;
        rhob(i,j,k) = results.rhob;
        
        end
    end
    disp(['T = ',num2str(T(i)),'C: Complete'])
end

%pcolor(H2O_sys, T, squeeze(eps_m(:,:,1))); hold on
magmaRes = {};
for K = 1:length(P)
    C = contourc(H2O_sys, T, squeeze(rhob(:,:,K)),[rhob_fixed, rhob_fixed]);
    
    C(:,1) = [];
    for J = 1:size(C,2)
        [results] = MELTS_MT_Rho(C(1,J),C(2,J),P(K),c,rhoh,volflag,concentration,cf);
        C(3,J) = results.meltfrac;
        C(4,J) = results.volfrac;
        C(5,J) = results.meltH2O;
        C(6,J) = results.rhob;
    end
    
    magmaRes{K} = C';
end





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
