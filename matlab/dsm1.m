% Deep Space Maneuver Theory (1)
% C: 29JUN20


% inclusion of parking orbit will help model flyby altitude and reduces DV 



clear; clf; format long g; clc;

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

% Leveraging Orbit Elements         ***VINF SET SIMULAR TO JUNO MISSION***
T = T_e.*K + (37.5*86400);

a = ((sqrt(mu_s)/(2*pi))*T)^(2/3);

Vp = sqrt(mu_s*(2/aukm - 1./a));  
Ve = sqrt(mu_s*(2/aukm - 1./aukm));
Vinflaunch = Vp - Ve

% Vary launch Vinf value to see how trend changes in DSM DV


rp = aukm;
ra = ((rp^2)*(Vp^2))/(2*mu_s - rp*(Vp^2));
e = (ra - rp)/(ra + rp);

Va = Vp*rp/ra;

Pl = 2*pi*sqrt((((ra + rp)/2)^3) / mu_s);  % Eqn. (4) from 1997 Paper is Incorrect
                                           % Should be identical to T from above

%dVap = 0.08827500; %dThetaStar = 1; % 1:1 rev
%dVap = 0.396; %dTheta = 1;

dsmDV = 0.0001:0.005:1.00;
itt = 1;
%while abs(dTheta) > tol
for i=1:length(dsmDV)
    
    dVap = dsmDV(i);
    
    Var = Va - dVap;    % dVap = DV applied at the aphelion for leveraging orbit

    ar = ((2/ra) - (Var^2 / mu_s))^(-1);
    Pr = ((2*pi)/sqrt(mu_s))*(ar^(3/2));
    er = (ra/ar) - 1;
    rpr = ar*(1-er);

    
    Ver = sqrt(mu_s * ((2/aukm) - (1/ar)));
    fpa = acos((ra*Var)/(aukm*Ver));
    fpadeg(i) = fpa*(180/pi);
    % ----------------------------- Sims 1994 -----------------------------
    Vminus = Ver;
    Vinfre = (Vminus^2 + Ve^2 - 2*Vminus*Ve*cos(fpa))^(0.5);
    beta = asin((Vminus/Vinfre)*sin(fpa));
    %delta = 60*pi/180;
    delta = 2*asin(1/(1+(((6378+200)*Vinfre^2)/mu_e)));
    deltadeg = delta*180/pi;
    Vplus = (Ve^2 + Vinfre^2 - 2*Ve*Vinfre*cos(delta+beta))^(0.5);
    % ---------------------------------------------------------------------
    
    costhetae = ((ar*(1-er^2)/aukm)-1)/er;
    thetae = acos(costhetae);
    if p < 1
        thetae = -thetae;
    end
    thetae;
    
    tanthetae = sqrt((1-costhetae)/(1+costhetae));
    tanEe2 = sqrt((1-er)/(1+er))*tanthetae;
    
    Ee = 2*atan(tanEe2);

    tep = sqrt((ar^3)/mu_s)*(Ee - er*sin(Ee));

    Te = (Pl/2) + (Pr/2) + p*tep;
    Te = Te/T_e;
    Tesaved(i) = Te;

    tauE = Te*2*pi;
    

    thetaE = tauE - 2*pi*K;

    dTheta(i) = thetae - thetaE;
    
    %dVap = dVap + 0.00005;
    %itt = itt + 1;
end

figure(1)
name = ['DSM \DeltaV versus \Delta\theta^* ','K=',num2str(K),' P=',num2str(p),' L=',num2str(L)];
yyaxis left
plot(dsmDV,dTheta)
line([0 1],[0 0],'color','k')
title(name,'fontsize',18)
yyaxis left
ylabel('\Delta\theta^* = \theta^*_{e} - \theta^*_{E}','fontsize',16)
yyaxis right
plot(dsmDV,fpadeg)
ylabel('Flight Path Angle \gamma','fontsize',16);
xlabel('DSM \DeltaV (km/s)','fontsize',16);
grid on; box on; set(gcf,'color','w');

figure(2)
title('EGA Encounter Time','fontsize',18)
yyaxis left
plot(dsmDV,Tesaved)
ylabel('T_e TOF to EGA from Launch (years)','fontsize',16);
yyaxis right
line([0.5 1],[Vinflaunch Vinflaunch],'color','k')
ylabel('Launch Vinf','fontsize',16)
xlabel('DSM \DeltaV (km/s)','fontsize',16);
grid on; box on; set(gcf,'color','w');


%Vinflaunch
%fpadeg
%disp('_____________________')
%dTheta
%disp('_____________________')
%Te
%itt
%dVap


