%% Main file for running PSF fit simulations
% specifies parameters and runs psfFitSimulations for several different dipole orientations
% Output is printed to command window!

% Add folder and subfolders to path
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

% signal  photon count
photon = 5*10^5;
% amount of background counts per pixel (mean)
noise = 100^3;
% angle_uncertainty specifies amount of uncertainty in orientation
% estimation
% 0  .. exact knowledge
% 1 .. normal distribution with std 2 (for theta and phi)
angle_uncertainty = 0;
reducedExcitation=0; % set to 0 (no reduced excitation) or 1 (reduced excitation)
phi = pi/4; % azimuthal angle
astigmatism = 0.7; % astigmatism coefficient

%================================================================

input = [phi, noise, astigmatism, angle_uncertainty,reducedExcitation];
disp('Parameters used')
disp(['Astigmatism: ', num2str(astigmatism)])
disp(['Photon count: ', num2str(photon)])
disp(['Background noise: ', num2str(noise), ' counts per pixel'])
disp(['Angle uncertainty: ' , num2str(angle_uncertainty)])
if reducedExcitation==0
    disp('Reduced excitation: no')
else
    disp('Reduced excitation: yes')
end

% Run simulations
disp('Simulation with theta = pi/2')
psfFitSimulations(photon,pi/2,input);

disp('Simulation with theta = pi/3')
psfFitSimulations(photon,pi/3,input);

disp('Simulation with theta = pi/6')
psfFitSimulations(photon,pi/6,input);