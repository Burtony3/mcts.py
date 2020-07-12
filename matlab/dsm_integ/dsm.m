%% Deep Space Maneuver Table Script
%  C: 11JUL20


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

%% DSM Inputs
% Inputs
K = 2;          % Earth Orbit Revolutions
thetaList = 47.7; % List of Theta Values to Solve
pltFBvecs = false;
pltLevOrb = true;

out = calcDsm(K,thetaList,pltFBvecs,pltLevOrb)
