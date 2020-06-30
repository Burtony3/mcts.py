% Deep Space Maneuver Theory (1)
% C: 29JUN20


% inclusion of parking orbit will help model flyby altitude and reduces DV 



clear; close; format long g; clc;

% ____________________________________________________________________________________
% Inputs
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

% Leveraging Orbit Elements
T = T_e.*K + (6.5*86400);

a = ((sqrt(mu_s)/(2*pi))*T)^(2/3);

Vp = sqrt(mu_s*(2/aukm - 1./a));  
Ve = sqrt(mu_s*(2/aukm - 1./aukm));
Vinflaunch = Vp - Ve;


rp = aukm;
ra = ((rp^2)*(Vp^2))/(2*mu_s - rp*(Vp^2));
e = (ra - rp)/(ra + rp);

Va = Vp*rp/ra;

Pl = 2*pi*sqrt((((ra + rp)/2)^3) / mu_s);  % Eqn. (4) from 1997 Paper is Incorrect
                                           % Should be identical to T from above

%dVap = 0.08827500; %dThetaStar = 1; % 1:1 rev
dVap = 0.35; %dTheta = 1;

itt = 1;
%while abs(dTheta) > tol
for i=1
    Var = Va - dVap;    % dVap = DV applied at the aphelion for leveraging orbit

    ar = ((2/ra) - (Var^2 / mu_s))^(-1);
    Pr = ((2*pi)/sqrt(mu_s))*(ar^(3/2));
    er = (ra/ar) - 1;
    rpr = ar*(1-er);

    
    Ver = sqrt(mu_s * ((2/aukm) - (1/ar)));
    fpa = acos((ra*Var)/(aukm*Ver));
    fpadeg = fpa*(180/pi);
    Vminus = Ver;
    Vinfre = (Vminus^2 + Ve^2 - 2*Vminus*Ve*cos(fpa))^(0.5);
    beta = asin((Vminus/Vinfre)*sin(fpa));
    %delta = 60*pi/180;
    delta = 2*asin(1/(1+(((6378+200)*Vinfre^2)/mu_e)));
    deltadeg = delta*180/pi
    Vplus = (Ve^2 + Vinfre^2 - 2*Ve*Vinfre*cos(delta+beta))^(0.5)

    thetae = acos(((ar*(1-er^2)/aukm)-1)/er);

    tanEe2 = sqrt((1-er)/(1+er))*tan(thetae/2);
    Ee = 2*atan(tanEe2);

    tep = sqrt((ar^3)/mu_s)*(Ee - er*sin(Ee));

    Te = (Pl/2) + (Pr/2) + p*tep;
    Te = Te/T_e;

    tauE = Te*2*pi;
    

    thetaE = tauE - 2*pi*K;

    dTheta = thetae - thetaE;
    
    dVap = dVap + 0.00005;
    itt = itt + 1;
end

Vinflaunch
fpadeg
disp('_____________________')
dTheta
disp('_____________________')
Te
%itt
%dVap


