function sigf = melt_rho(P,T,w,flag)
% P input in MPa
% T input in C
% w input in wt%

if flag == 1 %RHYOLITE Guo et al. (2016)
    sigf = (10.^(2.983-0.0732*w-(3528-233.8*w+(0.763-0.0075*w.^2)*(P))./(T+273)));
    
elseif flag == 2 %DACITE Laumonier et al. (2015)
    dv = -0.176*w-0.388;
    Ea = 6146*w-88400;
    s1 = (-0.064*w+5.96)+P*10*(0.0000106*w+0.0000249);
    sigf = exp(s1+(Ea+P*10*dv)./(8.314*(T+273))); %Laumonier 2015

    
elseif flag == 3 %ANDESITE Laumonier et al. (2017)
    s1 = exp((-0.34*w+8.96)+P*10*(-0.00000807*w+0.000167));
    Ea = -9627*w+125000;
    dv = -0.146*w+2.462;
    sigf = s1*exp(-((Ea+P*10*dv)/(8.314*(T+273))));    
    
elseif flag == 4 %ANDESITE Guo et al. (2017)
    
    sigf = 10.^(5.23-0.56*w.^0.6-(8130.4-1462.7*w.^0.6+(581.3-12.7*w.^2)*(P/1000))./(T+273)); %Guo 2017
 
elseif flag == 5 %RHYOLITE Gaillard (2004)
    
    a = -78.9*log(w)+754;
    b = -2925*log(w)+64132;
    sigf = a.*exp((-b+2*P)./(8.3144621*(T+273)));
    
elseif flag == 6 %BASALT Ni (2011)
    
    sigf = 10^(2.172-(860.82-204.46*sqrt(w))/(T+273-1146.8));

end