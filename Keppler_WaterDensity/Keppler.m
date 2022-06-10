function [sigf] = Keppler(P,T,con,flag)
%
% Function to compute the resistivity of high-pressure and high-temperature
% fluids using different empirical results from Hans Keppler's research
% group
%
% Usage: [sigf] = keppler(P,T,con,flag)
%
%   Inputs: P = Pressure in MPa
%           T = Temperature in Celsius
%           con = Concentration in wt% (can be NaCl or HCl wt% equivalent)
%           flag = 1, 2, or 3 to choose which empirical relation to use
%
%   flag = 1 is for Sinmyo and Keppler (2017) which investigated
%   NaCl-bearing aqueous fluids to 1000 MPa and 600 C. Their figures
%   extrapolate between 250 C to 1400 C.
%
%   flag = 2 is for Guo and Keppler (2019) which investigated NaCl-bearing
%   aqueous fluids to 5000 MPa and 900 C. Their figures range from 1000 MPa
%   to 5000 MPa and temperatures from 300 C to 900 C
%
%   flag = 3 is for Klumbach and Keppler (2020) which investigated
%   HCl-bearing aqueous fluids to 1000 GPa and 700 C. Their figures
%   extrapolate to 100 MPa and temperatures between 150 C to 800 C.
%
%   Pers. comm. with Hans Keppler suggested using either Sinmyo and Keppler
%   or Klumbach and Keppler in volcanic environments and he suggested lower
%   pressures are valid.
%
% Outputs:
%       sigf = Conductivity (linear) in S/m
%
%

%Constants
molar_mass = 18.02; %molar mass of H2O in g/mol
molar_concentration_NaCl = (con*10)/58.44; %in mol/L; molar mass of NaCl = 58.44
molar_concentration_HCl = (con*10)/36.458; %in mol/L; molar mass of HCl = 36.458

c_hcl = con*(molar_concentration_HCl/molar_concentration_NaCl); %NaCl-equivalent HCl concentration

%Pitzer Sterner (1994) H2O Equation of State from pers. comm. with Hans Keppler
% See: https://publish.uwo.ca/~awither5/fugacity/index.htm
% Codes can be found here: https://publish.uwo.ca/~awither5/fugacity/fugacity.py
% Script is in Python which is wrapped here in MATLAB:
mod = py.importlib.import_module('fugacity_PitzerSterner');

%The function to compute density/volume of pure H2O is PSvolume:
%       PSvolume(P,T), with P in bars and T in Kelvin
%   Output is cubic centimeters per mol
cc_per_mol = double(mod.PSvolume(P*10,T+273)); %volume of pure H2O
rho = molar_mass*(1/cc_per_mol); %Convert to density of H2O in g/cc

if flag == 1
    %Sinmyo and Keppler, 2017
    lambda0 = 1573 - 1212*rho + 537062/(T+273) - 208122721/(T+273).^2; %for NaCl
    sigf = 10.^(-1.7060 - 93.78/(T+273)+3.0781*log10(rho)+0.8075*log10(con)+log10(lambda0));

elseif flag == 2
    %Guo and Keppler, 2019
    lambda0 = 1573 - 1212*rho + 537062/(T+273) - 208122721/(T+273).^2; %for NaCl
    sigf = 10.^(-0.919-872.5/(T+273)+7.61*log10(rho)+0.852*log10(con)+log10(lambda0));

elseif flag == 3
    %Klumbach and Keppler, 2020
    lambda0 = 2550.14 -505.1*rho - 429437/(T+273); %for HCl
    sigf = 10.^(-2.032+205.8/(T+273)+3.888*log10(rho)+0.895*log10(c_hcl)+log10(lambda0));

end

