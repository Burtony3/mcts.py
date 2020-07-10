function dTheta = dvreq(dsmDV,vinf)
%DVREQ Summary of this function goes here
%   Detailed explanation goes here

    K = 3;          % Earth Orbit Revolutions
    %L = 1;          % K:L(M) L: S/C Orbit Revolutions
    p = 1;         % Crossing before perihelion (-1) or after (1)

    mu_s = 132712401800;
    mu_e = 000000398600;
    aukm = 149600000;
    T_e = ((2*pi)/sqrt(mu_s))*(aukm^(3/2));
    a_e = aukm;
    e_e = 0.0;


    Ve = sqrt(mu_s*(2/aukm - 1./aukm));
    Vp = vinf + Ve;

    rp = aukm;
    ra = ((rp^2)*(Vp^2))/(2*mu_s - rp*(Vp^2));
    e = (ra - rp)/(ra + rp);

    Va = Vp*rp/ra;

    Pl = 2*pi*sqrt((((ra + rp)/2)^3) / mu_s);
    a = ((sqrt(mu_s)/(2*pi))*Pl)^(2/3);

    dVap = dsmDV;
    Var = Va - dVap;
    
    ar = ((2/ra) - (Var^2 / mu_s))^(-1);
    Pr = ((2*pi)/sqrt(mu_s))*(ar^(3/2));
    er = (ra/ar) - 1;
    rpr = ar*(1-er);

    
    Ver = sqrt(mu_s * ((2/aukm) - (1/ar)));
    fpa = acos((ra*Var)/(aukm*Ver));
    
    costhetae = ((ar*(1-er^2)/aukm)-1)/er;
    thetae = acos(costhetae);
    if p < 0
        thetae = -thetae;
    end
    
    tanthetae = sqrt((1-costhetae)/(1+costhetae));
    tanEe2 = sqrt((1-er)/(1+er))*tanthetae;
    
    Ee = 2*atan(tanEe2);
    %if p > 0
    %    Ee = -Ee;
    %end    

    tep = sqrt((ar^3)/mu_s)*(Ee - er*sin(Ee));

    Te = (Pl/2) + (Pr/2) + p*tep;
    Te = Te/T_e;

    tauE = Te*2*pi;
    
    thetaE = tauE - 2*pi*K;

    dTheta = abs(thetae - thetaE);
     
end

