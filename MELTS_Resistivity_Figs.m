%MELTS RESISTIVITY FIGURES
%
% This script reproduces figures for three case studies in manuscript
% entitled "Constraining magma reservoir conditions by integrating 
% thermodynamic petrological models and bulk resistivity from \
% magnetotellurics" by Cordell, D., Naif, S., Troch, J., Huber, C.
%
%
% The scripts link bulk resistivity to MELTS thermodynamic modelling
%

clearvars; close all


case_study = 2; % 0 = synthetic; 1 = Mono Basin; 2 = Newberry; 3 = LdMVF


if case_study == 0
    %Synthetic
    w = 0.25:0.01:6; %*dissolved water content* in wt% (in the melt)
    H2O_sys = w; %Total system water content in wt%
    T = 700:1:1000; %Temperature in C
    P = 100; %Pressure in MPa
    
    rhoh_fixed = 1000; %Crystal or host resistivity in Ohm m
    c = 1.5; % Modified Archie's Law connectivity parameter 1.0 < c < 1.5 for melt
    labflag = 1; %Flag for empirical melt resistivity relationship (1 = Guo et al. (2016))

    % T1, w1, T2, and w2 are melt and temperature ranges from petrological
    % analyses of past eruptive products. For the synthetic case study we
    % don't have any constraints so they are NaN.
    T1 = NaN; T2 = NaN;
    w1 = NaN; w2 = NaN;
    
    % Compute the range of melt resistivity values given the w and T
    % constraints. For the synthetic case they are NaN.
    q1 = 1/melt_rho(P, T1, w1, labflag);
    q2 = 1/melt_rho(P, T2, w2, labflag);
    
    % This is the bulk resistivity of the anomaly in Ohm m
    rhob_fixed = 2;

elseif case_study == 1
    %Mono Basin
    w = 0.25:0.01:7;
    H2O_sys = 0.25:0.01:3;
    T = 700:1:950;
    P = 250;
    
    rhoh_fixed = 1000;
    c = 1.05;
    labflag = 1;

    T1 = 750; T2 = 850;
    w1 = 4; w2 = 5;
    q1 = 1/melt_rho(P, T1, w1, labflag);
    q2 = 1/melt_rho(P, T2, w2, labflag);
    
    rhob_fixed = [6.5 3 10];

elseif case_study == 2
    %Newberry
    w = 0.25:0.01:7;
    H2O_sys = w;
    T = 650:1:950;
    P = 100;
    
    rhoh_fixed = 1000;
    c = 1.5;
    labflag = 1;

    T1 = 850; T2 = 850;
    w1 = 1.5; w2 = 1.5;
    q1 = 1/melt_rho(P, T1, w1, labflag);
    q2 = 1/melt_rho(P, T2, w2, labflag);
    
    rhob_fixed = [38 25 50];

elseif case_study == 3
    %LdMVF
    w = 0.25:0.01:7;
    H2O_sys = w;
    T = 650:1:1100;
    P = 100;
    
    rhoh_fixed = 1000;
    c = 1.5;
    labflag = 1;

    T1 = 760; T2 = 1000;
    w1 = 4; w2 = 5;
    q1 = 1/melt_rho(P, T1, w1, labflag);
    q2 = 1/melt_rho(P, T2, w2, labflag);
    
    rhob_fixed = 0.3;
   
end


%Compute melt conductivity for the w-T parameter space------------------

sigm = nan(length(w),length(T)); %Melt conductivity in S/m
for i = 1:length(w)
    for j = 1:length(T)
        sigm(i,j) = melt_rho(P, T(j), w(i), labflag);
    end
end


%Compute bulk conductivity for the chi-sigm parameter space------------
sigm_fixed = 1./logspace(-1,1,200); %Melt resistivity
chi = 0:0.005:1; %melt fraction

sigb = nan(length(chi),length(sigm_fixed));
for i = 1:length(sigm_fixed)
    for j = 1:length(chi)
        %Using Modified Archie's Law (MAL)
        sigb(j,i) = MAL(0.001,sigm_fixed(i),c,chi(j));
      
    end
end


%Compute melt fraction for w-T parameter space at a fixed rhob
phi = nan(length(w),length(T));
for i = 1:length(w)
    for j = 1:length(T)

        %For a fixed rhob, the melt conductivity is a "middleman" between
        %melt fraction and w-T parameter space. You can skip the middleman
        %and just compute phi as a function of w-T parameter space
        sigm_middleman = melt_rho(P, T(j), w(i), labflag);
        phi(i,j) = MAL_solve_phi(rhob_fixed(1),rhoh_fixed,1./sigm_middleman,c);


    end

end

%Compute volatile resistivity as a function of temperature
concentration1 = 1;
concentration2 = 16;
rhof_param = nan(length(T),2);
for i = 1:length(T)
    rhof_param(i,1) = 1./Keppler(P,T(i),concentration1,1);
    rhof_param(i,2) = 1./Keppler(P,T(i),concentration2,1);
end

%-----------------------------------------------------------------------
% FIGURES FOR W-T PARAMETER SPACE: MELT RESISTIVITY, BULK RESISTIVITY, AND MELT FRACTION
markersize = 12;
screensize=get(groot,'Screensize');
fig = figure(1);
clf
set(fig,'Position',[0.1*screensize(3) 0.1*screensize(4) 0.8*screensize(3) 0.4*screensize(4)])


%W-T Parameter Space for Melt Resistivity++++++++++++++++++++++++++++++
subplot(1,3,1)
v = 10.^(-4:1:4);
contour(w,T,1./sigm',v,'-k','LineWidth',2,'ShowText','on'); hold on

v = (2:9)'*(10.^(-4:1:4));
contour(w,T,1./sigm',v(:),'--k','LineWidth',0.5,'ShowText','on'); hold on

ylabel('Temperature (°C)')
xlabel('Dissolved Water Content (wt%)')
title('(a) Melt Resistivity (\Omegam)')
set(gca,'FontSize',14)

%Plot a patch of petrological constraints
p = patch([w1 w2 w2 w1 w1],[T1 T1 T2 T2 T1],[0.5 0.5 0.5]);
p.FaceAlpha = 0.2;

%Plot max and min melt resistivity
plot(w1,T1,'ok','MarkerFaceColor','b','MarkerSize',markersize)
plot(w2,T2,'ok','MarkerFaceColor','r','MarkerSize',markersize)

%Rhom-Phi Parameters Space for Rhob+++++++++++++++++++++++++++++++++++
subplot(1,3,2)
v = 10.^(-4:1:4);
contour(1./sigm_fixed,chi,1./sigb,v,'-k','LineWidth',2,'ShowText','on'); hold on

v = (2:9)'*(10.^(-4:1:4));
contour(1./sigm_fixed,chi,1./sigb,v(:),'--k','LineWidth',0.5,'ShowText','off')
set(gca,'XScale','log')

%For the case studies, we find the max and min melt fraction based on
% melt resistivity constraint
if case_study ~= 0 
    maxmeltorig = nan(length(rhob_fixed),1);
    minmeltorig = nan(size(maxmeltorig));
    for i = 1:length(rhob_fixed)
        if i==1
            lw = 4; %Linewidth variable
        else
            lw = 2;
        end

        v = [-1 rhob_fixed(i) 10000];
        h = contour(1./sigm_fixed,chi,1./sigb,v,'-m','LineWidth',lw,'ShowText','off'); hold on

        melt1 = InterX(h,[q1 q1; 0 1]);
        melt2 = InterX(h,[q2 q2; 0 1]);

        if sum([~isempty(melt1) ~isempty(melt2)])==2
            
            if i == 1
                subplot(1,3,2)
                plot(melt1(1),melt1(2),'ok','MarkerFaceColor','b','MarkerSize',markersize)
                plot(melt2(1),melt2(2),'ok','MarkerFaceColor','r','MarkerSize',markersize)
            end
            %Max and Min melt fraction estimate using the "original" or
            %"standard" interpretation with no thermodynamic constraint
            maxmeltorig(i) = melt1(2);
            minmeltorig(i) = melt2(2);
        end



    end
end

xlabel('Melt Resistivity (\Omegam)')
ylabel('Melt Fraction')
title('(b) Bulk Resistivity (\Omegam)')
set(gca,'FontSize',14)

plot([q1 q1],[0 1],'-b','LineWidth',2)
plot([q2 q2],[0 1],'-r','LineWidth',2)

phi(phi==1) = 1.1;

cmap = parula(22);
cmap = cmap(2:2:end,:);

%Reparameterizing the above parameter spaces++++++++++++++++++++++++++++
subplot(1,3,3)
pcolor(w,T,phi'); colormap(cmap); hold on
hcb = colorbar; shading flat
hcb.Label.String = 'Melt Fraction';
hcb.Ticks = 0:0.1:1;
hcb.TickLabels = num2cell(0:0.1:1);
xlabel('Dissolved H2O wt%')
ylabel('Temperature (°C)')
caxis([0 1.1])

title(['(c) Fixed \rho_b = ',num2str(rhob_fixed(1)),' \Omegam'])

v = (0:0.1:1);
contour(w,T,phi',v,':k','LineWidth',1,'ShowText','on')

set(gca,'FontSize',14)

p = patch([w1 w2 w2 w1 w1],[T1 T1 T2 T2 T1],[0.5 0.5 0.5]);
p.FaceAlpha = 0.1;

plot(w1,T1,'ok','MarkerFaceColor','b','MarkerSize',markersize)
plot(w2,T2,'ok','MarkerFaceColor','r','MarkerSize',markersize)

if case_study == 3
    %Plot volatile resistivity as a function of T+++++++++++++++++++++++++++
    figure(3);
    plot(T,rhof_param(:,1),'-k','LineWidth',3); hold on
    plot(T,rhof_param(:,2),'--k','LineWidth',3);
    xlabel('Temperature (°C)')
    ylabel('Resistivity of Volatile Phase (\Omegam)')
    set(gca,'FontSize',20)
    legend('1 wt% NaCl Equivalent','16 wt% NaCl Equivalent')
    grid on
    axis([min(T) max(T) 0 10])

end

%----------------------------------------------------------------------
% MELTS COUPLING----------------------------------------------------------

[X,Y] = meshgrid(H2O_sys,T);

%Set up all the connectivity parameters for the melt and volatile phase
if case_study == 1 || case_study == 0
    
    %For Mono Basin we don't compute the saturated regime
    concentration = 0;
    mm = c; %Melt connectivity
    mf = 10; %Ignore volatile phase
    
elseif case_study == 2
    
    %For Newberry we have 1 wt% NaCl equivalent volatile
    concentration = 1;
    mm = c; %Melt connectivity
    mf = 2; %Poorly-connected volatile (e.g. bubbles)
    
elseif case_study == 3
    
    %For LdMVF we have a 16 wt% NaCl equivalent volatile
    concentration = 16;
    mm = c; %Melt connectivity
    mf = 1; %Well-connected volatile (e.g. brine lens)
    
end

%Compute crystal frac (eps_c), melt frac (eps_m), fluid volatile frac
%(eps_g) and melt dissolved H2O (H2O_m)
rhob = nan(length(T),length(H2O_sys));
rhom = nan(size(rhob));
eps_h = nan(size(rhob)); %Crystal volume fraction from MELTS
eps_m = nan(size(rhob)); %Melt volume fraction from MELTS
eps_f = nan(size(rhob)); %Volatile volume fraction from MELTS
H2O_m = nan(size(rhob)); %Dissolved melt content wt% from MELTS
rhof = nan(length(T),1); %Volatile resistivity
for i = 1:length(T)
    
    %Volatile resistivity computed using Sinmyo and Keppler (2017)
    rhof(i) = 1./Keppler(P,T(i),concentration,1);
    
    for j = 1:length(H2O_sys)

        %Get all the outputs from MELTS using the MELTS parameterization
        [eps_h(i,j),eps_m(i,j),eps_f(i,j),H2O_m(i,j)] = volfrac(P,T(i),H2O_sys(j));
        
        %Calculate melt resistivity from Guo et al. (2016) using MELTS outputs
        rhom(i,j) = 1./melt_rho(P,T(i),H2O_m(i,j),labflag);

        %Calculate bulk resistivity using three-phase MAL
        rhob(i,j) = 1./MAL3(1/rhoh_fixed,1/rhom(i,j),1/rhof(i),eps_h(i,j),eps_m(i,j),eps_f(i,j),mm,mf);
        %rhob(i,j) = 1./MAL(1/rhoh_fixed,1/rhom(i,j),c,eps_m(i,j));
    end
end

%Post-processing to remove invalid melt fractions <0
cutoff = 0;
rhob(eps_m<=cutoff) = NaN;
eps_f(eps_m<=cutoff) = NaN;
eps_m(eps_m<=cutoff) = NaN;

if case_study <= 1 %For synthetic case and Mono Basin ignore saturated regime
    rhob(eps_f>0) = NaN;
end

%Ignore dissolved melt in saturated regime (i.e. assume it is at the
%saturation level)
H2O_m(eps_f>0) = NaN;

%Set "liquidus" to 0.99 melt
ind0 = find(eps_h<=0.01);
eps_h(ind0) = 0;
eps_m(ind0) = 1 - eps_f(ind0);

eps_m(eps_m==1) = 1.1;


%PLOT MELTS PARAMETER SPACE--------------------------------------------
figure(2)
pcolor(H2O_sys,T,(eps_m)); view([0 90]); hold on

% v = [-1 0.01 1000];
% contour(H2O_sys,T,eps_m,v,'--w','LineWidth',3)
colormap(cmap);
hcb = colorbar; shading flat
hcb.Label.String = 'Melt Fraction';
hcb.Ticks = 0:0.1:1;
hcb.TickLabels = num2cell(0:0.1:1);
xlabel('Total System H2O wt%')
ylabel('Temperature (°C)')
axis([min(H2O_sys) max(H2O_sys) min(T) max(T)])
caxis([0 1.1])

%--------------------SATURATION CURVE----------------------
v = [-1 10^-4 10];
S = contour(H2O_sys,T,eps_f,v,'-k','LineWidth',2,'ShowText','off');

%--------------------LIQUIDUS CURVE--------------------
v = [-1 10^-2 10];
contour(H2O_sys,T,eps_h,v,'--k','LineWidth',2,'ShowText','off')

set(gca,'FontSize',14)

v = 0:10;
[C,h] = contour(H2O_sys,T,(H2O_m),v,'--','ShowText','on','LineColor',[0.2 0.2 0.2]);
clabel(C,h,'Color',[0.2 0.2 0.2],'FontName','Courier')

if case_study ~= 0 
v = [-1 w1 1000];
[w1contour] = contour(H2O_sys,T,(H2O_m),v,'LineStyle','none');
w1contour = w1contour(:,2:end)';

v = [-1 w2 1000];
[w2contour] = contour(H2O_sys,T,(H2O_m),v,'LineStyle','none');
w2contour = w2contour(:,2:end)';

end

[C1,h1] = contour(H2O_sys,T,(eps_f),'--','ShowText','on','LineColor',[1 0 0]);
clabel(C1,h1,'Color',[1 0 0],'FontName','Courier')

v = 10.^(-2:2)'*(1:9);
v = v(:);
contour(H2O_sys,T,(rhob),v,'-k','ShowText','on')


%Get the maximum and minimum values of melt fraction, temperature,
%dissolved water content, etc. along the fixed bulk resistivity contour
for rhocount = 1:length(rhob_fixed)
    v = [-1 rhob_fixed(rhocount) 10000];
    if rhocount == 1
        lt = '-';
    else
        lt = '--';
    end
    figure(2);
    A = contour(H2O_sys,T,(rhob),v,[lt,'w'],'ShowText','off','LineWidth',3);

    A(:,A(1,:)==rhob_fixed(rhocount))=[];

    %finalx = []; finalmelt = []; finalvolatile = []; finalwater = []; finaltemp = [];
    for i = 1:length(A(1,:))

        [finalx(i,rhocount),finalmelt(i,rhocount),finalvolatile(i,rhocount),finalwater(i,rhocount)] = volfrac(P,A(2,i),A(1,i));
        finaltemp(i,rhocount) = A(2,i);

    end
    
    if rhocount==1 && case_study ~= 2
        figure(1); subplot(1,3,3);plot(finalwater(finalvolatile<=0,rhocount),finaltemp(finalvolatile<=0,rhocount),'-w','LineWidth',3)
    end
    
    if case_study == 1
    
        indtemp = intersect(find(finaltemp(:,rhocount)>T1), find(finaltemp(:,rhocount)<T2));
        indwater = intersect(find(finalwater(:,rhocount)>w1), find(finalwater(:,rhocount)<w2));
        indpoly = intersect(indtemp,indwater);

        [maxmelt(rhocount) indmaxmelt] = max(finalmelt(indpoly,rhocount));
        maxmeltT(rhocount) = finaltemp(indpoly(indmaxmelt),rhocount);
        maxmeltw(rhocount) = finalwater(indpoly(indmaxmelt),rhocount);

        [minmelt(rhocount) indminmelt] = min(finalmelt(indpoly,rhocount));
        minmeltT(rhocount) = finaltemp(indpoly(indminmelt),rhocount);
        minmeltw(rhocount) = finalwater(indpoly(indminmelt),rhocount);

    else
        
        indpoly = 1:size(A,2);
        [maxmelt(rhocount) indmaxmelt] = max(finalmelt(indpoly,rhocount));
        maxmeltT(rhocount) = finaltemp(indpoly(indmaxmelt),rhocount);
        maxmeltw(rhocount) = finalwater(indpoly(indmaxmelt),rhocount);
        maxmeltv(rhocount) = finalvolatile(indpoly(indmaxmelt),rhocount);

        [minmelt(rhocount) indminmelt] = min(finalmelt(indpoly,rhocount));
        minmeltT(rhocount) = finaltemp(indpoly(indminmelt),rhocount);
        minmeltw(rhocount) = finalwater(indpoly(indminmelt),rhocount);
        minmeltv(rhocount) = finalvolatile(indpoly(indminmelt),rhocount);
    end
    
end

% figure(2)
% if case_study == 0  
%     text(2.5,750,'Saturated Partial Melt','FontWeight','bold')
%     text(1.5,850,'Under-Saturated','FontWeight','bold')
%     text(1.7,835,'Partial Melt','FontWeight','bold')
%     text(2.2,950,'Super-liquidus Melt','FontWeight','bold')
%     text(4.5,905,'Saturated','FontWeight','bold')
%     text(4.3,890,'Super-Liquidus','FontWeight','bold')
%     text(4.6,875,'Melt','FontWeight','bold')
% end




