
%clear; clc

K = 2;


offsetDays = 0;

mu_s = 132712401800;
mu_e = 000000398600;
aukm = 149600000;
T_e = ((2*pi)/sqrt(mu_s))*(aukm^(3/2));
Ve = sqrt(mu_s*(2/aukm - 1./aukm));



% Setup Launch Orbit
T_f = T_e.*K + (offsetDays*86400);
a = ((sqrt(mu_s)/(2*pi))*T_f)^(2/3);

Vp = getVel(mu_s, aukm, a);
Vinflaunch_f = Vp - Ve;




function v = getVel(mu_bdy, r_dist, sma)
    v = sqrt(mu_bdy*(2/r_dist - 1/sma));  
end
