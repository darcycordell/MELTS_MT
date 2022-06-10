%MELTS RESISTIVITY SUPPLEMENTARY FIGURES
%
% This script reproduces supplementary figures for in the manuscript
% entitled "Constraining magma reservoir conditions by integrating 
% thermodynamic petrological models and bulk resistivity from \
% magnetotellurics" by Cordell, D., Naif, S., Troch, J., Huber, C.
%
%
% The scripts examine the effects of various parameters on the final
% results in the paper
%

% SUPPLEMENTARY FIGURE 1: COMPARISON GUO AND GAILLARD-------------------
clearvars;

T = 700:1000; %Temperature in C
P = 200; % Pressure in MPa
w = 0:0.1:7; %Dissolved water content in wt%

sigf1 = nan(length(w),length(T)); sigf2 = nan(size(sigf1));
for i = 1:length(w)
    for j = 1:length(T)
        sigf1(i,j) = melt_rho(P, T(j), w(i), 1); %Guo et al. 2016
        sigf2(i,j) = melt_rho(P, T(j), w(i), 5); %Gaillard, 2004
    end
end

% Percent Difference in Resistivity between two relations
dsig = (1./sigf2 - 1./sigf1)./(1./sigf1);

figure(1)
subplot(2,1,1)
v = 10.^(-4:1:4);
contour(w,T,1./sigf1',v,'-k','LineWidth',2,'ShowText','on'); hold on

v = (2:9)'*(10.^(-4:1:4));
contour(w,T,1./sigf1',v(:),'--k','LineWidth',0.5,'ShowText','on'); hold on

v = 10.^(-4:1:4);
[C1,h1] = contour(w,T,1./sigf2',v,'-r','LineWidth',2,'ShowText','on'); hold on
clabel(C1,h1,'Color','r')

v = (2:9)'*(10.^(-4:1:4));
[C2,h2] = contour(w,T,1./sigf2',v(:),'--r','LineWidth',0.5,'ShowText','on'); hold on
clabel(C2,h2,'Color','r')

ylabel('Temperature (^oC)')
xlabel('Dissolved Water Content (wt%)')

title('Melt Resistivity (\Omegam)')
set(gca,'FontSize',14)

manual_legend('Guo et al. (2016)','-k','Gaillard (2004)','-r');

subplot(2,1,2)
v = 100*[-4:0.5:4];
contour(w,T,100*dsig',v,'-k','LineWidth',2,'ShowText','on'); hold on

v = 100*[-4:0.1:1];
contour(w,T,100*dsig',v,'--k','LineWidth',1,'ShowText','on'); hold on

ylabel('Temperature (^oC)')
xlabel('Dissolved Water Content (wt%)')
grid on

title('% Difference in Melt Resistivity (\Omegam)')
set(gca,'FontSize',14)

%% SUPPLEMENTARY FIGURE 2: PRESSURE EFFECTS GUO--------------------------
clearvars;

T = 850;
w = 0:0.1:7;

%Change pressure from 100 MPa to 300 MPa
sigf1 = melt_rho(100, T, w, 1);
sigf2 = melt_rho(300, T, w, 1);

figure(2); 
plot(w,1./sigf1,'LineWidth',2); hold on
plot(w,1./sigf2,'LineWidth',2)
xlabel('Dissolved Water Content (wt%)')
ylabel('Melt Resistivity (\Omegam)')
legend('100 MPa','300 MPa')
grid on
set(gca,'FontSize',14)
set(gca,'XLim',[0 7])

%% SUPPLEMENTARY FIGURE 3 AND 7: CRYSTAL RESISTIVITY AND EFFECTS----------
clearvars; 

T = 700:1000;
chi = 0:0.01:1;
w = 4;
P = 100;

rhom = exp(-log(2.9)+8581.8./(T+273));
rhof = 1./melt_rho(P,T,w,1);

%Plot of crystal resistivity as a function of T
figure(3);
semilogy(T,rhom,'-k','LineWidth',2)
xlabel('Temperature (°C)')
ylabel('Crystal Resistivity (\Omegam)')
set(gca,'FontSize',14)
grid on

sigb1 = nan(length(rhom),length(chi)); sigb2 = nan(size(sigb1));
for i = 1:length(rhom)
    for j = 1:length(chi)
        
        %Compute bulk resistivity with variable crystal rho
        sigb1(i,j) = MAL(1./rhom(i),1./rhof(i),1.5,chi(j));
        
        %Compute bulk resistivity with fixed crystal rho
        sigb2(i,j) = MAL(1./1000,1./rhof(i),1.5,chi(j));
    end
end

%Plot the effect of crystal resistivity
figure(7)
v = 10.^(-4:1:4);
contour(T,chi,1./sigb1',v,'-k','LineWidth',2,'ShowText','on'); hold on
contour(T,chi,1./sigb2',v,'--r','LineWidth',2,'ShowText','on'); hold on

v = (2:9)'*(10.^(-4:1:4));
contour(T,chi,1./sigb1',v(:),'-k','LineWidth',0.5,'ShowText','off'); hold on
contour(T,chi,1./sigb2',v(:),'--r','LineWidth',0.5,'ShowText','off'); hold on
%set(gca,'XScale','log')
xlabel('Temperature (°C)')
ylabel('Melt Fraction')
set(gca,'FontSize',14)
grid on
manual_legend('Variable \rho_c','-k','Fixed \rho_c','--r');

%% SUPPLEMENTARY FIGURES 5 AND 6: COMPARE MIXING RELATIONS-----------------
%This section compares three mixing relations: Modified Archie's Law with 
% m = 1.15, MAL with m = 1.5 and Hashin-Shtrikman upper bound (HS+)

clearvars;

plot_fig = 1;

%Varying the connectivity in MAL from 1.15 to 1.5
c = 1.5;
c1 = 1.15;
rhom_fixed = 1000;
sigf_fixed = 1./logspace(-1,1,201);
sigb_fixed = 1./logspace(0,2,200);

phi = nan(length(sigf_fixed),length(sigb_fixed));
phihs = nan(size(phi));
philow = nan(size(phi));

for i = 1:length(sigf_fixed)
    for j = 1:length(sigb_fixed)
        
        if sigf_fixed(i) > sigb_fixed(j)
        
            phi(i,j) = MAL_solve_phi(1./sigb_fixed(j),rhom_fixed,1./sigf_fixed(i),c);
            [phihs(i,j)] = HS_solve_phi(1./sigb_fixed(j),rhom_fixed,1./sigf_fixed(i));
            philow(i,j) = MAL_solve_phi(1./sigb_fixed(j),rhom_fixed,1./sigf_fixed(i),c1);
            
            
        end
        
    end
end

%Difference between each mixing relation 
% MAL with m = 1.15, MAL with m = 1.5 and Hashin-Shtrikman (HS).
d1 = (phi-phihs);

d2 = (phi-philow);

d3 = (philow-phihs);

%Get absolute difference between each combo
d(:,:,1) = abs(d1);
d(:,:,2) = abs(d2);
d(:,:,3) = abs(d3);


length(find(d3<0))

%Get the color-coding id where 1 is difference between MAL(m=1.5) and HS+,
%2 is difference between MAL(m=1.5) and MAL(m=1.15), and 3 is difference
%between MAL(m=1.15) and HS+.
[maxval,id] = max(d,[],3);

id(isnan(phi)) = NaN;

%Plot absolute difference color-coded by the combo id
figure(6)
pcolor(1./sigf_fixed,1./sigb_fixed,id'); hold on; shading flat
v = 0:0.01:1;
contour(1./sigf_fixed,1./sigb_fixed,maxval',v,'-k','LineWidth',2,'ShowText','on'); hold on
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('Melt Resistivity (\Omegam)')
title(['Maximum Difference in Melt Fraction = ',num2str(max(maxval(:)))])
ylabel('Bulk Resistivity (\Omegam)')
set(gca,'layer','top')

hcb = colorbar;
colormap(jet(3));
hcb.Ticks = 0:1:2;
caxis([0 2]);
hcb.TickLabels = {['MAL(',num2str(c1),')-HS+'],['MAL(',num2str(c),')-MAL(',num2str(c1),')'],'MAL(1.5)-HS+'};

set(gca,'FontSize',14)

grid on
grid minor

%
%Another way to show the difference is to plot as contours at specific
%fixed bulk resistivity values as a function of melt fraction and melt
%resistivity
chi = 0:0.005:1;

sigb = nan(length(chi),length(sigf_fixed));
sigbHS = nan(size(sigb));
sigb_lowc = nan(size(sigb));
% sigb_varc = nan(size(sigb));
for i = 1:length(sigf_fixed)
    for j = 1:length(chi)
        sigb(j,i) = MAL(1/rhom_fixed,sigf_fixed(i),c,chi(j));
        sigbHS(j,i) = HS(1/rhom_fixed,sigf_fixed(i),chi(j));
        sigb_lowc(j,i) = MAL(1/rhom_fixed,sigf_fixed(i),c1,chi(j));
        
        if chi(j) > 0.4
            
            mvar(j,i) = 1;
            sigb_varc(j,i) = MAL(1/rhom_fixed,sigf_fixed(i),1,chi(j));
        else
            
            mvar(j,i) = -2.75*chi(j) +2.1;
            sigb_varc(j,i) = MAL(1/rhom_fixed,sigf_fixed(i),mvar(j,i),chi(j));
        end
        
    end
end


figure(5)
v = 10.^(-4:1:4);
contour(1./sigf_fixed,chi,1./sigb,v,'-k','LineWidth',2,'ShowText','on'); hold on
contour(1./sigf_fixed,chi,1./sigb_lowc,v,'--r','LineWidth',2,'ShowText','on'); hold on
contour(1./sigf_fixed,chi,1./sigbHS,v,':b','LineWidth',2,'ShowText','on'); hold on
contour(1./sigf_fixed,chi,1./sigb_varc,v,'-g','LineWidth',2,'ShowText','on'); hold on


set(gca,'XScale','log')

xlabel('Melt Resistivity (\Omegam)')
ylabel('Melt Fraction')
title('Bulk Resistivity (\Omegam)')
manual_legend(['MAL (m = ',num2str(c),')'],'-k',['MAL (m = ',num2str(c1),')'],'--r',['MAL (m = variable)'],'-g','HS+',':b');

set(gca,'FontSize',14)

grid on
grid minor


%% SUPPLEMENTARY FIGURE 4: PRESSURE EFFECT MVP RESISTIVITY----------------
clearvars; 
T = 650:1000;

col = {'k','b','r'};
figure(4);
rhog1 = []; rhog2 = [];
for j = [1 3]
concentration1 = 1;
concentration2 = 16;
P = 100;

for i = 1:length(T)
    rhog1(i,j) = 1./Keppler(P,T(i),concentration1,j);
    rhog2(i,j) = 1./Keppler(P,T(i),concentration2,j);
end


semilogy(T,rhog1(:,j),['-',col{j}],'LineWidth',3); hold on
semilogy(T,rhog2(:,j),['--',col{j}],'LineWidth',3);

end

xlabel('Temperature (^oC)')
ylabel('Resistivity of Volatile Phase (\Omegam)')
set(gca,'FontSize',20)
legend('1 wt% SK17','16 wt% SK17','1 wt% KK20','16 wt% KK20')
grid on

