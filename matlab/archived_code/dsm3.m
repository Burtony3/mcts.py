    
clear; clc;


K = 2;          % Earth Orbit Revolutions
L = 1;          % K:L(M) L: S/C Orbit Revolutions
p = 1;         % Crossing before perihelion (-1) or after (1)
tol = 10e-5;    % Itteration Tolerance

% ____________________________________________________________________________________
% Constants
mu_s = 132712401800;
mu_e = 000000398600;
aukm = 149600000;
T_e = ((2*pi)/sqrt(mu_s))*(aukm^(3/2));
a_e = aukm;
e_e = 0.0;


T = T_e.*K;

a = ((sqrt(mu_s)/(2*pi))*T)^(2/3);

Vp = sqrt(mu_s*(2/aukm - 1./a));
Ve = sqrt(mu_s*(2/aukm - 1./aukm));
Vinflaunch = Vp - Ve;


rp = aukm;
rai = ((rp^2)*(Vp^2))/(2*mu_s - rp*(Vp^2));
Va = Vp*rp/rai;
ar = ((2/rai) - (Va^2 / mu_s))^(-1);


gammaminusIn = 0:1:45;

for i=1:length(gammaminusIn)

Vminus = Vp;
gammaminus = gammaminusIn(i)*(pi/180);


Vinfre = (Vminus^2 + Ve^2 - 2*Vminus*Ve*cos(gammaminus))^(0.5);
beta = asin((Vminus/Vinfre)*sin(gammaminus));
delta = 2*asin(1/(1+(((6378+200)*Vinfre^2)/mu_e)));
deltadeg = delta*180/pi;
Vplus = (Ve^2 + Vinfre^2 - 2*Ve*Vinfre*cos(delta+beta))^(0.5);
gammaplus = asin((Vinfre/Vplus)*sin(beta+delta));
anew = 1/(2/aukm - (Vplus^2)/mu_s);
e = (((aukm*(Vplus^2)/mu_s)-1)^2 * cos(gammaplus) + (sin(gammaplus))^2)^(0.5);
ra = anew*(1+e);
ra = ra/aukm;

gp(i) = gammaplus*(180/pi);
vplus(i) = Vplus;
raNew(i) = ra;
end


plot(gammaminusIn, raNew)


