function out = calcDsm(K,thetaInt,pltFBvecs,pltLevOrb,offsetDays)
%CALCDSM Summary of this function goes here
%   Detailed explanation goes here

if K<0; p=-1; else; p=1; end; K = abs(K);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

% _________________________________________________________________________
% Defining Constants
mu_s = 132712401800;
mu_e = 000000398600;
r_e = 6378;
alt_min = 200;
aukm = 149597870.7;
T_e = ((2*pi)/sqrt(mu_s))*(aukm^(3/2));
a_e = aukm;
e_e = 0.0;
% _________________________________________________________________________
% Pre-DSM Leveraging Orbit Elements
T = T_e.*K + (offsetDays*86400);
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
% _________________________________________________________________________
% Find DSM Dt Value for Returning Leg
fixedintercepttime = (K*T_e/2 + p*(thetaInt/360)*T_e);
thetaIndeg = thetaInt;

% _________________________________________________________________________
% Given a Dt Value, Compute Return Leg of Leveraging Maneuver and State
ddti_val =  fixedintercepttime;
outputDiff = getLambert(ddti_val,thetaIndeg,K,p,mu_s,aukm);
thetaIn = outputDiff.thetaIn;
tof_toEGA = outputDiff.dti/86400;   % from DSM to EGA
dvDsmvec = outputDiff.xD2(4:6) - xD1(4:6)
dvDsm = norm(dvDsmvec);

% _________________________________________________________________________
% Body Intercept State
xIn = outputDiff.xIn;       % Earth's State at Intercept
xD2 = outputDiff.xD2;       % SC State right after DSM DV maneuver
xIn2 = outputDiff.xInSC;    % SC State at Earth Intercept

% _________________________________________________________________________
% PLANAR Flyby Calculations
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
    bendingdeg = bending*180/pi;
end
%bendingdeg
%r_car
vinf1ang = acos(dot(vinf1hat,[1;0]));
%vinf1angdeg = vinf1ang*180/pi
if p<0
    vinf2ang = vinf1ang - bending;
else
    vinf2ang = -vinf1ang + bending;
end
%vinf2ang*180/pi
vinf2 = vinf1mag.*[cos(vinf2ang); sin(vinf2ang)];
vinf2hat = vinf2/norm(vinf2);
vp2 = xIn(4:5)+vinf2;
xOut1 = [aukm*cos(thetaIn); aukm*sin(thetaIn); 0; vp2(1); vp2(2); 0];
kOut1 = conv_carKep(mu_s,xOut1,0);

disp(' ')
disp('DV DSM')
disp(dvDsm)
disp('Vinf launch')
disp(Vinflaunch)
disp('Total DV')
disp(Vinflaunch + dvDsm)
disp('vinf mag')
disp(norm(vinf2))


% Debugging Flyby Vectors Plot
if pltFBvecs
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

% _________________________________________________________________________
% Integration of new trajectory and Plotting
preManeuverState = tbp(xL,T/2,mu_s,0,options);
postManeuverState = tbp(xD2,fixedintercepttime,mu_s,0,options);
postFBState = tbp(xOut1,12*365*86400, mu_s,0,options);

% DSM Trajectory Visual
if pltLevOrb
    clf
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
    pltCirc(0,0,2.872e9)
    hold off
    axis equal; grid on; 
end

% _________________________________________________________________________
% calcDsm Results
out = struct;
out.solfound = true;
out.K = K;
if p<0
    out.thetaInt = -thetaInt;
else
    out.thetaInt = thetaInt;
end
out.dsmDV = dvDsm;
out.dsmDVvec = dvDsmvec;
out.departOrbitTime = (T/2)/(86400);
out.levOrbitTime = tof_toEGA + out.departOrbitTime;
out.vinflaunch = Vinflaunch;
out.vinf1 = vinf1;
out.vinf2 = vinf2;
out.bendingdeg = bendingdeg;
out.r_car = r_car;
out.dvegaRA = kOut1.ra;


% Itterating Function for Dti and DV_DSM
function outputDiff = getLambert(ddti_val,thetaIndeg,K,p,mu,aukm)
    if p>0
        thetaIn_f = (thetaIndeg-90)*p*(pi/180);
        xIn_f = [aukm*cos(thetaIn_f); aukm*sin(thetaIn_f); 0; -Ve*sin(thetaIn_f); Ve*cos(thetaIn_f); 0];
    else
        thetaIn_f = (thetaIndeg+90)*p*(pi/180);
        xIn_f = [aukm*cos(thetaIn_f); aukm*sin(thetaIn_f); 0; -Ve*sin(thetaIn_f); Ve*cos(thetaIn_f); 0];
    end

    lambcallf = l0(2,xD1,xIn_f,ddti_val,mu);
    xD2f = [0; kL.ra; 0; lambcallf(5,1); lambcallf(5,2); lambcallf(5,3)];

    outputDiff = struct;
    outputDiff.xIn = xIn_f;
    outputDiff.thetaIn = thetaIn_f;
    outputDiff.dti = ddti_val;
    outputDiff.lambert = lambcallf;
    outputDiff.xD2 = xD2f;
    outputDiff.xInSC = [xIn_f(1); xIn_f(2); xIn_f(3); lambcallf(6,1); lambcallf(6,2); lambcallf(6,3)];

end

end

% Aux Req. Functions
function v = getVel(mu_bdy, r_dist, sma)
    v = sqrt(mu_bdy*(2/r_dist - 1/sma));  
end

function out = conv_carKep(mu,x,debug)
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
        debug = 0;
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
        %disp('Theta NaN, changed to 0');
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
    
    if debug == 1 
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