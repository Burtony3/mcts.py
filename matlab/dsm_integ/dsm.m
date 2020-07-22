%% Deep Space Maneuver Table Script
%  C: 11JUL20

clear; clc; format long g; firstPC = 0; PC = 0; firstMac = 1; Mac = 1; 


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

%% DSM Inputs
Klist = [-2 2 -3 3 -4 4];          % Leveraging Maneuver Type (K:1 w/ Earth, pos or neg.) 
thetaList = [2:5:70];      % List of intercept thetas
%Klist = -3;
%thetaList = 40.7;
pltFBvecs = false;        % Debug Flyby Velocity Vectors
pltLevOrb = false;        % Debug Trajectory Visual

%% DSM Calculation
allDsm = struct;
for i=1:length(Klist)
    for j=1:length(thetaList)
        K = Klist(i);
        theta = thetaList(j);
        
        % Save Variable Name
        if K < 0
            strucName = ['k',num2str(abs(K)),'m'];
            
        else
            strucName = ['k',num2str(abs(K)),'p'];
        end
        thetaName = matlab.lang.makeValidName(num2str(theta)); thetaName = thetaName(2:end);


        % Calculate DSM and Leveraging Orbit Properties
        out = calcDsm3(K,theta,pltFBvecs,pltLevOrb);
        %out = calcDsm2(K,theta,pltFBvecs,pltLevOrb,20);
        
        
        if out.solfound
            allDsm.(strucName).(['t',thetaName]) = out;
        else
            disp('Solution did not compute')
            %allDsm.(strucName).(['t',num2str(theta)]) = NaN;
        end 
    end
    disp(i)
end

% Nominal Launch Vinf Values
% 2: 5.0779
% 3: 6.9272
% 4: 7.9273


%% Processing Results
if false
    allDsmfld = fieldnames(allDsm);
    u=0;
    for i=1:length(allDsmfld)    
        outK = allDsm.(allDsmfld{i});
        fld = fieldnames(outK);

        for j=1:length(fld)     
            csvArray(u+j,1) = outK.(fld{j}).K;
            csvArray(u+j,2) = outK.(fld{j}).thetaInt;
            csvArray(u+j,3) = outK.(fld{j}).totalDV;
            csvArray(u+j,4) = outK.(fld{j}).dsmDV;
            csvArray(u+j,5) = outK.(fld{j}).levOrbitTime;
            csvArray(u+j,6) = outK.(fld{j}).vinflaunch;
            csvArray(u+j,7) = outK.(fld{j}).vinf1(1);
            csvArray(u+j,8) = outK.(fld{j}).vinf1(2);
            csvArray(u+j,9) = outK.(fld{j}).vinf2(1);
            csvArray(u+j,10) = outK.(fld{j}).vinf2(2);
            csvArray(u+j,11) = outK.(fld{j}).bendingdeg;
            csvArray(u+j,12) = outK.(fld{j}).r_car;
            csvArray(u+j,13) = outK.(fld{j}).dvegaRA;
        end
        u = u+j;
    end   
    % Writing Table and Exporting to .csv
    csvTable = array2table(csvArray,'VariableNames', ...
        {'K','Theta','TotalDV','DSMDV','dt','vinfLaunch','vinf1x','vinf1y', ...
        'vinf2x','vinf2y','bendang','fb_radius','dvegaRA'});
    writetable(csvTable,'dsm_ega.csv');
    
end


%% Figure of Leveraging Orbit
if true
plottitle = '$\Delta$VEGA Performance';

SMA_J = 778.6e6*1000;   % m         v_esc = 18.4635 km/s
SMA_S = 1.433e9*1000;   % m         v_esc = 13.6097 km/s
SMA_U = 2.872e9*1000;   % m         v_esc = 9.61344 km/s
SMA_N = 4.495e9*1000;   % m

levTypes = fieldnames(allDsm);
for i=1:length(levTypes)
    
    thetaValues = fieldnames(allDsm.(levTypes{i}));

    for j=1:length(thetaValues)
        % CalcDSM2 Function
        dsmdv = allDsm.(levTypes{i}).(thetaValues{j}).dsmDV;
        vinfdv = allDsm.(levTypes{i}).(thetaValues{j}).vinflaunch;
        dv(j,i) = dsmdv+vinfdv;
        ra(j,i) = (allDsm.(levTypes{i}).(thetaValues{j}).dvegaRA)/149600000;
        
        % CalcDSM3 Function
        %dv(j,i) = allDsm.(levTypes{i}).(thetaValues{j}).totalDV;
        %ra(j,i) = (allDsm.(levTypes{i}).(thetaValues{j}).dvegaRA)/149600000;
    end
end

hold on
for i=1:6
    if i==1 || i==3 || i==5
        plot(dv(:,i),ra(:,i),'r','linewidth',1.5)
    else
        plot(dv(:,i),ra(:,i),'b','linewidth',1.5)
    end
end

SMA = [778.6e6;1.433e9;2.872e9;4.495e9]./149600000;
for i=1:length(SMA)
   plot([0 10],[SMA(i) SMA(i)],'k','linewidth',2)  
end
hold off
set(gca, 'YScale', 'log')
xlim([4.5 10])


ax = gca; ax.FontSize = 14;
set(gca, 'TickLabelInterpreter','Latex');
xlabel('Total $\Delta$V ($km/s$)','fontsize',12,'Interpreter','Latex');
ylabel('Final Aphelion Radius (AU)','fontsize',14,'Interpreter','Latex');
title(plottitle,'fontsize',14,'Interpreter','Latex');
InSet = get(ax, 'TightInset');
set(ax, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
grid on; box on; set(gcf,'color','w');
% print(gcf,'-dpng','-r150',name);

end
    