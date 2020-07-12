

clear; clc; close; format long g; firstPC = 0; PC = 0; firstMac = 1; Mac = 1; 

%% Initialize Local Orbital Mechanics Source Files
if firstPC == 1
    cd2 = 'C:\Users\Aurora\dev\dev_orbitalmechanics\orbitalmechanics_rp';
    addpath(genpath(cd2))
    cd3 = 'C:\Users\Aurora\dev\dev_orbitalmechanics\research';
    addpath(genpath(cd3))   
end
if firstMac == 1
    cd2 = '/Users/rohanpatel/dev/dev_orbitalmechanics/orbitalmechanics_rp';
    addpath(genpath(cd2))
    cd3 = '/Users/rohanpatel/dev/dev_orbitalmechanics/research';
    addpath(genpath(cd3))
end
if PC == 1
    cspice_kclear;
    nf009 = 'C:\Users\Aurora\dev\dev_orbitalmechanics\orbitalmechanics_rp\src\SPKs\naif0009.tls';
    de438 = 'C:\Users\Aurora\dev\dev_orbitalmechanics\orbitalmechanics_rp\src\SPKs\de438.bsp';
    pck = 'C:\Users\Aurora\dev\dev_orbitalmechanics\orbitalmechanics_rp\src\SPKs\pck00010.tpc';
    cspice_furnsh({de438,nf009,pck})
end
if Mac == 1
    cspice_kclear;
    nf009 = '/Users/rohanpatel/dev/dev_orbitalmechanics/orbitalmechanics_rp/src/SPKs/naif0009.tls';
    de438 = '/Users/rohanpatel/dev/dev_orbitalmechanics/orbitalmechanics_rp/src/SPKs/de438.bsp';
    pck = '/Users/rohanpatel/dev/dev_orbitalmechanics/orbitalmechanics_rp/src/SPKs/pck00010.tpc';
    cspice_furnsh({de438,nf009,pck})   
end

clear firstPC firstMac PC Mac

%% new sec dsm
% Inputs
K = 2;          % Earth Orbit Revolutions
p = 1;          % Crossing before perihelion (-1) or after (1)

% ____________________________________________________________________________________
% Constants
mu_s = 132712401800;
mu_e = 000000398600;
r_e = 6378;
alt_min = 200;
aukm = 149600000;
T_e = ((2*pi)/sqrt(mu_s))*(aukm^(3/2));
a_e = aukm;
e_e = 0.0;

% Leveraging Orbit Elements

T = T_e.*K + (10*86400);

a = ((sqrt(mu_s)/(2*pi))*T)^(2/3);


Vp = getVel(mu_s, aukm, a);
Ve = getVel(mu_s, aukm, a_e);
Vinflaunch = Vp - Ve;


xi = [0; -aukm; 0; Ve; 0; 0];
ki = conv_carKep(mu_s, xi, 0);

xL = [0; -aukm; 0; Vp; 0; 0];
kL = conv_carKep(mu_s, xL, 0);

Vd1 = getVel(mu_s, kL.ra, a);
xD1 = [0; kL.ra; 0; -Vd1; 0; 0];





% DSM DeltaV and DTime Calculations
% -------------------------------------------------------------------------
thetaIndeg = 28;
ddti_val =  79.4963388442993;
% -------------------------------------------------------------------------

if p>0
    thetaIn = (thetaIndeg-90)*p*(pi/180);
    xIn = [aukm*cos(thetaIn); aukm*sin(thetaIn); 0; Ve*cos(thetaIn); -Ve*sin(thetaIn); 0]
else
    thetaIn = (thetaIndeg+90)*p*(pi/180);
    xIn = [aukm*cos(thetaIn); aukm*sin(thetaIn); 0; -Ve*cos(thetaIn); Ve*sin(thetaIn); 0];
end

dti = K*T_e/2 + p*T_e*(thetaIndeg/100);
if p<0
    dti = (dti + ddti_val*86400);
else
    dti = (dti - ddti_val*86400);
end

lambcall = l0(2,xD1,xIn,dti,mu_s);
dvDsm = norm(lambcall(1,1:3));


% Body Intercept State and PLANAR Flyby Calculations
xD2 = [0; kL.ra; 0; -Vd1+dvDsm; 0; 0];
xIn2 = [aukm*cos(thetaIn); aukm*sin(thetaIn); 0;lambcall(end,1);lambcall(end,2);lambcall(end,3)];
vinf1 = xIn2(4:5) - xIn(4:5);
vinf1mag = norm(vinf1);
vinf1hat = vinf1/norm(vinf1);
vphat = xIn2(4:5)/norm(xIn2(4:5));
vehat = xIn(4:5)/norm(xIn(4:5));
bending = acos(dot(vehat,vinf1hat));
bendingdeg = bending*180/pi;
r_car = (mu_e/(vinf1mag^2))*(-1 + (1/sin(bending/2)));
% Actual Bending Achieved due to Min. FB Constraint
if r_car < r_e+alt_min
    r_car = r_e+alt_min;
    bending = 2*asin(1/(1+(r_car*vinf1mag^2/mu_e)));
    bendingdeg = bending*180/pi
end
vinf1ang = -acos(dot(vinf1hat,[1;0]));
vinf1angdeg = vinf1ang*180/pi;
if p<0
    vinf2ang = vinf1ang - bending;
else
    vinf2ang = vinf1ang + bending;
end
vinf2 = vinf1mag.*[cos(vinf2ang); sin(vinf2ang)];
vinf2hat = vinf2/norm(vinf2);
vp2 = xIn(4:5)+vinf2;
xOut1 = [aukm*cos(thetaIn); aukm*sin(thetaIn); 0; vp2(1); vp2(2); 0];
out = conv_carKep(mu,x,p)

% Debugging Flyby Vectors
if false
   figure
   hold on
   quiver(0,0,xIn(4),xIn(5),'autoscale','off')
   quiver(0,0,xIn2(4),xIn2(5),'autoscale','off')
   quiver(xIn(4),xIn(5),vinf1(1),vinf1(2),'autoscale','off')
   quiver(xIn(4),xIn(5),vinf2(1),vinf2(2),'autoscale','off')
   quiver(0,0,vp2(1),vp2(2),'autoscale','off')
   hold off
   legend('Ve','Vsc','Vinf1','Vinf1','Vsc2')
   axis equal; grid on;
end



% Integration of new trajectory
addtlproptime = 0;
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
postManeuverState = tbp(xD2,dti+addtlproptime,mu_s,0,options);
postFBState = tbp(xOut1,3.5*365*86400, mu_s,0,options);
preManeuverState = tbp(xL,(K*T_e/2),mu_s,0,options);

r_postMS = postManeuverState(end,1:3)';
r_xIn = xIn(1:3);
r_diff = norm(r_postMS - r_xIn);




if true
    figure
    hold on
    scatter(0,0,'MarkerEdgeColor','r')
    scatter(xi(1),xi(2),'MarkerEdgeColor','b')
    scatter(xD1(1), xD1(2),'MarkerEdgeColor', 'k')
    scatter(xIn(1), xIn(2),'MarkerEdgeColor', [0.5 0.5 0.5])
    plot(postManeuverState(:,1),postManeuverState(:,2))
    plot(preManeuverState(:,1),preManeuverState(:,2))
    plot(postFBState(:,1),postFBState(:,2))
    pltCirc(0,0,aukm)  
    pltCirc(0,0,778.6e6)
    pltCirc(0,0,1.433e9)
    hold off
    axis equal; grid on; 
end


function v = getVel(mu, r, a)
    v = sqrt(mu*(2/r - 1/a));  
end

function out = conv_carKep(mu,x,p)
% KEPELEM: Given a mu and state vector x (6x1), calculates orbital elements
%          output: Structure with fields: a; e; i; raan; aop; ta; t; rp; ra]
%               h    = Angular Momentum [norm]
%               a    = Semi Major Axis (distance)
%               e    = Eccentricity (none)
%               e_vec= Eccentricity Vector 
%               i    = inclination (radians)
%               raan = RAAN (radians)
%               aop  = Argument of Periapsis (radians)
%               ta   = True Anomaly (radians)
%               t    = Period (seconds)
%               rp   = Apoapsis (distance)
%               ra   = Periapsis (distance)
%           Source: Curtis 4.4 - alg4.2 p197bk

    if nargin < 3
        p = 0;
    end


    r_ = x(1:3); v_ = x(4:6);
    vr = dot(v_,r_)/norm(r_);
    h = cross(r_,v_); h_ = norm(h);
    
    i = acos(h(3)/h_);
    if (pi/2 < i) && (i < 3*pi/2)
        i_ = 'retrograde';
    else
        i_ = 'prograde';
    end
    
    N = cross([0;0;1],h);
    omega = acos(N(1)/norm(N));
    if N(2) < 0
        omega = 2*pi -omega;
    end
    
    e = (1/mu)*((norm(v_)^2 - mu/norm(r_))*r_ - norm(r_)*vr*v_);
    e_ = sqrt(dot(e,e));
    
    w = acos(dot(N,e)/(norm(N)*e_));
    if e(3) < 0
        w = 2*pi -w;
    end
    
    theta = acos(dot(e,r_)/(e_*norm(r_)));
    if vr<0
       theta = 2*pi - theta; 
    end
    if isnan(theta)
        disp('Theta NaN, changed to 0');
        theta = 0;
    end
    
    rp = (h_^2/mu)*(1/(1+e_));
    ra = (h_^2/mu)*(1/(1-e_));
    a = 0.5*(rp+ra);
    T = ((2*pi)/sqrt(mu))*a^(3/2);

    out = struct;
    out.h = h_;
    out.h_vec = h;
    out.a = a;
    out.e = e_;
    out.e_vec = e;
    out.i = i;
    out.raan = omega;
    out.aop = w;
    out.ta = theta;
    out.t = T;
    out.rp = rp;
    out.ra = ra;
    
    if p == 1

        
        disp('-------------------------------------------')
        disp('Input State Vector: [x;y;z;vx;vy;vz]')
        disp(' ')
        disp(x)
        disp('-------------------------------------------')
        disp(['Angular Momentum:    ',num2str(h_), ' km^2/s'])
        disp(' ')
        disp(['Inclination:         ',num2str(i*(180/pi)), ' deg'])
        disp(['Orbit Type:          ',i_])
        disp(' ')
        disp(['Eccentricity:        ',num2str(e_)])
        disp(['RAAN:                ',num2str(omega*(180/pi)), ' deg'])
        disp(['Argument of Perigee: ',num2str(w*(180/pi)), ' deg'])
        disp(['True Anomaly:        ',num2str(theta*(180/pi)), ' deg'])
        disp(' ')
        if e<1
            T_hrs = T/3600;
            disp(['Orbit Period:        ',num2str(T_hrs), ' hours'])
        else
            disp(['Orbit Period:        ','HYPERBOLIC'])
        end
        disp(['Periapsis Radius:    ',num2str(rp), ' km'])
        disp(['Apoapsis Radius:     ',num2str(ra), ' km'])
        disp(['Semimajor Axis:      ',num2str(a), ' km'])
        disp('-------------------------------------------')
    end
    
end

function pltCirc(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'b');
    hold off
end







