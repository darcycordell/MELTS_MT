function [phi] = MAL_solve_phi(rho_b,rho_h,rho_f,m)
% Function which solves for melt fraction using Modified Archie's Law
% given a known bulk resistivity, matrix resistivity, fluid resistivity,
% and connectivity exponent.
%
% Usage: [phi,i] = solve_phi_MAL(rho_b,rho_m,rho_f,m)
%
%   Inputs: rho_b = bulk resistivity (Ohm m)
%           rho_h = host rock resistivity (Ohm m)
%           rho_f = fluid resistivity (Ohm m)
%           m = connectivity parameter (m>0, unitless)
%
%   Outputs: phi = melt fraction (or porosity)
%            i = number of iterations to reach solution
%
% Modified Archie's Law (MAL) formula:
%
%   rho_b*rho_h*phi^m + rho_b*rho_f(1-phi)^p = rho_f*rho_h
%   
%       p = log10(1-phi^m)/log10(1-phi)
%
% The formula cannot be solved algebraically so the Newton-Raphson method
% is used to solve the equation numerically. Since 0<phi<1, a hard-coded
% starting guess of 0.5 is used.

% Check for input errors first:
flag = 0; indnan = [];
if any(rho_b > rho_h)
    disp('Bulk resistivity cannot be greater than matrix resistivity')
    flag = 1;
    indnan = [indnan find(rho_b>rho_h)];
end

if any(rho_b < rho_f)
    disp('Bulk resistivity cannot be less than fluid resistivity')
    flag = 1;
    indnan = [indnan find(rho_b<rho_f)];
end

if any(rho_f > rho_h)
    disp('Fluid resistivity cannot be greater than matrix resistivity')
    flag = 1;
    indnan = [indnan find(rho_f>rho_h)];
end

sigb = 1./rho_b;
sigf = 1./rho_f;
sigh = 1./rho_h;

phi = ((sigb-sigh)./(sigf-sigh)).^(1/m);


if flag == 1
    phi(indnan) = NaN;
end