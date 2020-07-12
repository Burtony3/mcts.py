%% Deep Space Maneuver Required Dt Calculation
%  C: 11JUL20

clear; clc; close all; format long g; firstPC = 0; PC = 0; firstMac = 1; Mac = 1; 

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

%% DSM Calculation Loop
%thetacase31minus = 0.711658617463823*(180/pi); % Jon Sims 1997 Paper Case

thetadeglist = 5:0.5:45;

Klist = [3];

% Looped Calculation for DDTI and DVReq.
for k=1:length(Klist)
    K = Klist(k);

    j=1;
    for i=1:length(thetadeglist)
        ddti = 40;
        thetadeg = thetadeglist(i);
        %options = optimset('PlotFcns',@optimplotfval);
        aF = @(ddti)dvreq2(ddti,thetadeg,K);
        [ddti_val, fVal] = fminsearch(aF,ddti);

        if fVal < 100000
            ddti_store(j,k) = ddti_val;
            fVal_store(j,k) = fVal;
            thetaFound(j,k) = thetadeg;
            global dvDsm; dvReq(j,k) = dvDsm;
            j=j+1;
            %disp(thetadeg)
        else
            ddti_store(j,k) = NaN;
            fVal_store(j,k) = NaN;
            thetaFound(j,k) = NaN;
            dvReq(j,k) = NaN;
            j=j+1;
            %disp(thetadeg)           
        end
        
    end
    disp(K)
    
    
end

% Plot Data
if false
    figure
    hold on
    plot(thetaFound(:,1), dvReq(:,1))
    plot(thetaFound(:,2), dvReq(:,2))
    hold off   
    legend({'$3:1^-$','$3:1^+$'},'fontsize',16,'location','northwest', ...
        'Interpreter','latex')
    xlabel('$\theta$ Intercept True Anomaly ($^\circ$)','fontsize',16,'Interpreter','latex');
    ylabel('Required DSM $\Delta$V ($km/s$)','fontsize',16,'Interpreter','latex');
    if K<0
        title(['$',num2str(abs(K)),':1^-$ DSM $\Delta$V for fixed Intercept True Anomaly $\theta$'],'fontsize',20,'Interpreter','latex')
    else
        title(['$',num2str(abs(K)),':1^+$ DSM $\Delta$V for fixed Intercept True Anomaly $\theta$'],'fontsize',20,'Interpreter','latex')
    end        
    ax = gca; ax.FontSize = 16;
    set(gca,'TickLabelInterpreter','latex')
    grid on; box on; set(gcf,'color','w');
    set(gca,'LooseInset',get(gca,'TightInset'));
    

    savefilename = ['/Users/rohanpatel/Desktop/','DV_DSM_', num2str(K)];
    %print(gcf,'-dpng','-r200',savefilename)
    
end

% Figure of Trajectories
if true
    figure
    for i=1:length(thetaFound)
        if isnan(dvReq(i,1))
        else
            pltDSMorbits(Klist(1),thetaFound(i),ddti_store(i))
        end
    end
    pltCirc(0,0,149600000)
    if Klist(1)<0
        title(['$',num2str(abs(K)),':1^-$ Leveraging Orbit Trajectory'],'fontsize',20,'Interpreter','latex')
    else
        title(['$',num2str(abs(K)),':1^+$ Leveraging Orbit Trajectory'],'fontsize',20,'Interpreter','latex')
    end        
    xlabel('X Distance (km)','fontsize',16,'Interpreter','latex');
    ylabel('Y Distance (km)','fontsize',16,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')
    grid on; box on; set(gcf,'color','w'); axis equal;
    set(gca,'LooseInset',get(gca,'TightInset'));
    
    

end


function pltCirc(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'b');
    hold off
end