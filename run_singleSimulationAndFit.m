% Script for simulation and fitting of PSF for fixed dipole emitters
% Create PSF and calculate fit using psfFitMLE

% Add folder and subfolders to path
folder = fileparts(which(mfilename));
addpath(genpath(folder));

%% Set parameters
photon = 5000; % signal photon count
theta = pi/3; % dipole inclination angle
phi = pi/3; % dipole azimuthal angle

% x and y position
% 0 corresponds to the center of the center pixel, 100 units are 1 pixel
xpos = 0;
ypos = 0;

noise = 0; % background mean photon count
defocus = 0; % defocus in nm
astigmatism = 0; % astigmatism coefficient (0.7 corresponds to 75nm RMS wavefront error)
Nx = 21; % simulation field size in pixels
reducedExcitation = 0; % no reduced excitation

%% Simulate PSF
psf = PSF(photon, theta, phi,xpos,ypos,[noise,defocus,astigmatism,reducedExcitation,Nx],'y');

%% Calculate fit
x = psfFitMLE(psf, noise, [theta,phi,astigmatism]);


%% Plot psf
% indicate actual position with black 'o' marker
% indicate estimated position with red 'x' marker
size = 12;
width = 2.5;
center = (length(psf)+1)/2;

imagesc(psf);
hold on
plot(center+xpos/100,center+ypos/100,'marker','o','MarkerSize',size,'Color','black','LineWidth', width)
hold on
plot(center+x(1)/100,center+x(2)/100,'marker','x','MarkerSize',size,'Color','red','LineWidth', width)
axis equal
axis tight
cb = colorbar;
ylabel(cb,'Intensity','FontSize',15)