
function [eps_x,eps_m,eps_g,H2O_m] = volfrac(P,T,H2O_sys)

%Other variables required for calculations
rho_x = 2900;   % density of bulk crystal load in kg/m3
rho_m = 2400;   % density of melt in kg/m3

%Calculate r_melt: Volume ratio melt/(melt+crystals)
A = -0.01006*erfc(-0.7464*(H2O_sys-2.841)) - 0.00639;
B = 250.4*exp(-0.9687*H2O_sys) + 726.7*exp(0.003698*H2O_sys);
r_melt = 0.5*erfc(A*(T-B));

%Calculate maximum water solubility
H2O_sat = (354.94*P^0.5 +9.623*P - 1.5223*P^1.5)/(T+273.15) + 0.0012439*P^1.5; % from Liu et al 2005

if r_melt == 0
    H2O_hyp = H2O_sys; %hypothetical melt water content if all water is hosted in melt
else
    H2O_hyp = H2O_sys* (r_melt * rho_m + (1-r_melt)*rho_x)/(r_melt*rho_m);
end


%Calculate vol fractions
rho_g = (-112.528*T^-0.381 + 127.811*(P*10)^-1.135 + 112.04*T^-0.411*(P*10)^0.033)*1e3; % EOS water

% initial estimate of volume fractions (no gas)
eps_m = r_melt;
eps_x = 1-r_melt;

% check for fluid saturation and calculate eps g
if H2O_hyp >= H2O_sat  
    H2O_m = H2O_sat;
    eps_g_guess = 0.0001;
    eps_x = r_melt*(eps_g_guess-1)+1-eps_g_guess;
    eps_m = 1-eps_g_guess - eps_x;
    rho_sys = eps_m*rho_m + eps_x*rho_x + eps_g_guess * rho_g;
    H2O_total_est = eps_g_guess * rho_g/rho_sys + eps_m*H2O_m/100* rho_m/rho_sys;
    while H2O_total_est < H2O_sys/100
        eps_g_guess = eps_g_guess + 0.0001;
        eps_x = r_melt*(eps_g_guess-1)+1-eps_g_guess;
        eps_m = 1-eps_g_guess - eps_x;
        rho_sys = eps_m*rho_m + eps_x*rho_x + eps_g_guess * rho_g;
        H2O_total_est = eps_g_guess * rho_g/rho_sys + eps_m*H2O_m/100* rho_m/rho_sys;
    end
    eps_g = eps_g_guess;     
else
    H2O_m = H2O_hyp;
    eps_g = 0;
end

eps_x = r_melt*(eps_g-1)+1-eps_g;
eps_m = 1-eps_g - eps_x;

if eps_m < 10^-3
    %eps_m = 0;
    H2O_m = 0;
end

end


