


clear; clc; clf; format long g;


%{
K = 2;          % Earth Orbit Revolutions
L = 1;          % K:L(M) L: S/C Orbit Revolutions
p = -1;         % Crossing before perihelion (-1) or after (1)
Vp = 34.9345817481521
Ve = 29.7844755752744
Vinflaunch = 5.15010617287766
min(abs(dTheta))
ans = 0.000325887253415313

j = 1;
for i=1:length(vinf_in)
    
    vinf = vinf_in(i);
    aF = @(DVap)dvreq(DVap,vinf);
    
    [dsmDVreq, fVal] = fminsearch(aF,dsmDVic);

    if fVal < tol
        v(j,1) = vinf;
        v(j,2) = dsmDVreq;
        v(j,3) = vinf + dsmDVreq;
        v(j,4) = fVal;
        j = j+1;
    end
end

clf;
figure(1)
plot(v(:,1),v(:,4))
title('Final d\theta Value')
xlabel('Vinf (km/s)')
ylabel('fVal')
grid on; box on; set(gcf,'color','w')

figure(2)
plot(v(:,2),v(:,1))
title('DSM DV Required')
ylabel('Departure V_{\infty} (km/s)')
xlabel('DSM DV (km/s)')
xlim([0 5]);
grid on; box on; set(gcf,'color','w')
%}




%thetacase31minus = 0.711658617463823*(180/pi)



thetadeglist = 5:0.5:45;

j=1;
for i=1:length(thetadeglist)
    ddti = 40;
    thetadeg = thetadeglist(i);
    %options = optimset('PlotFcns',@optimplotfval);

    aF = @(ddti)dvreq2(ddti,thetadeg);
    [ddti_val, fVal] = fminsearch(aF,ddti)
    
    if fVal < 100000
        ddti_store(j) = ddti_val;
        global dvDsm;

        dvReq(j) = dvDsm;
        j=j+1;
        disp(thetadeg)
    end

end


%plot(thetadeglist, ddti_store)