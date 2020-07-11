function [r_diff] = dvreq2(ddti,thetadeg)
%DVREQ2 Summary of this function goes here
%   Detailed explanation goes here
    global dvDsm

    K = 2;          % Earth Orbit Revolutions
    p = -1;         % Crossing before perihelion (-1) or after (1)


    mu_s = 132712401800;
    %mu_e = 000000398600;
    aukm = 149600000;
    T_e = ((2*pi)/sqrt(mu_s))*(aukm^(3/2));
    a_e = aukm;
    %e_e = 0.0;


    % Leveraging Orbit Elements

    T = T_e.*K + (10*86400);

    a = ((sqrt(mu_s)/(2*pi))*T)^(2/3);


    Vp = getVel(mu_s, aukm, a);
    Ve = getVel(mu_s, aukm, a_e);
    %Vinflaunch = Vp - Ve;


    xi = [0; -aukm; 0; Ve; 0; 0];
    %ki = conv_carKep(mu_s, xi, 0);

    xL = [0; -aukm; 0; Vp; 0; 0];
    kL = conv_carKep(mu_s, xL, 0);

    Vd1 = getVel(mu_s, kL.ra, a);
    xD1 = [0; kL.ra; 0; -Vd1; 0; 0];



    thetaIndeg = thetadeg;
    if p>0
        thetaIn = (thetaIndeg-90)*p*(pi/180);
    else
        thetaIn = (thetaIndeg+90)*p*(pi/180);
    end
    xIn = [aukm*cos(thetaIn); aukm*sin(thetaIn); 0; Ve*cos(thetaIn); -Ve*sin(thetaIn); 0];


    dti = K*T_e/2 + p*T_e*(thetaIndeg/100);
    if p<0
        dti = (dti + ddti*86400);
    else
        dti = (dti - ddti*86400);
    end
    

    lambcall = l0(2,xD1,xIn,dti,mu_s);

    dvDsm = norm(lambcall(1,1:3));
    xD2 = [0; kL.ra; 0; -Vd1+dvDsm; 0; 0];


    % Integration of new trajectory
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    postManeuverState = tbp(xD2,dti,mu_s,0,options);
    %preManeuverState = tbp(xL,(K*T_e/2));

    r_postMS = postManeuverState(end,1:3)';
    r_xIn = xIn(1:3);

    r_diff = abs(norm(r_postMS - r_xIn));

    
    function v = getVel(mu, r, sma)
        v = sqrt(mu*(2/r - 1/sma));  
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

end

