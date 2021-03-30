%% Simulations of PSF for dipoles

% authors: Philipp Zelger, Fabian Hinterer, Magdalena Schneider
% date: 18.06.2019
%--------------------------------------------------------------------------
% Simulation of point spread function for given dipole orientation of a
% fluorophore
% Cramer-Rao lower bound

function [psf,coordinates] = PSF(n_photon, theta, phi,xpos,ypos,parameters,shotnoise)
%% OUTPUT
% psf... matrix corresponding to PSF
% coordinates... lateral position of the fluorophore (x,y)

%% INPUT

% n_photon... photon count of one blink
% theta ... inclination angle of dipole
% phi ... azimuthal angle of dipole
% xpos ... horizontal position (0 = in the center of the center pixel, 100 units = 1 pixel)
% ypos ... vertical position
% shotnoise set to 'y' or 'n'

% function can be run as script, using the parameters set below
if nargin==0
    n_photon = 10^4;
    theta = pi/2;
    phi = 0;
    xpos = 0;
    ypos = 50;
    bgnoise = 0;
    def = 0;
    astigmatism = 0;
    exProb = 0;
    Nx = 17;
    shotnoise = 'y';
end


% additional input parameters
if nargin==7
    bgnoise = parameters(1);  % background noise (mean)
    def = parameters(2); % defocus in nm
    astigmatism = parameters(3); % astigmatism
    exProb = parameters(4); % reduced excitation, set to 1 (yes) or 0 (no)
    Nx = parameters(5); % field size in pixels
end


%% Set parameters
dipole = [theta,phi];

% Additional parameters
N=129;  % resolution in BFP, set to odd number!
lambda_0=680e-9; % wavelength
NA=0.7; % objective NA
RI=[1.33 1.46 1]; %refractive indices; RI=[RI_specimen, RI_intermed., RI_immoil]
d2=0;   %thickness of intermediate layer (layer 2)
z=0e-6;    % distance of dipole from medium 2
f=3e-3;   %focal length of objective
mu=1e-12;  %magnitude of dipole (arbitrary value)
dz=def*1e-9;  %defocus of objective lens (negative values mean moving the focus into the fluid)

% add oversampling factor for more accurate PSF
oversampling = 9 ;
ux=108/oversampling*1e-9; %resolution in focal space

%% Calculating defocus and aberrations
% Defocus function refers to refractive index n2
[~,Defocus,~] =fun_SA_RImismatch(N,RI(3),RI(3),NA,lambda_0,1);

uk=4*pi/lambda_0*NA/N; %unit in pupil space (k-space)
[~,~,R,pupil]=create_coord(N,1,'FFT'); % coordinate system

% considering aberrations
tilt = ux*uk*N/4;
Z_aberr2 = zeros(5,1);
Z_aberr2(1) = -1/100*(xpos)*oversampling*tilt; %x-tilt
Z_aberr2(2) = -1/100*(ypos)*oversampling*tilt; %y-tilt
Z_aberr2(5) = astigmatism; % vertical astigmatism
% (circle) Zernike polynomials
aberr=sum(ZernikeCalc([2 3 6], [Z_aberr2(1);Z_aberr2(2);Z_aberr2(5)], pupil,'NOLL'),3); %aberration file "coef.mat" must be loaded

%% Calculating BFP-fields for dipole orientation
[Ex_Px(:,:),Ey_Px(:,:)]=fun_dipole_imaging(N,lambda_0,NA,RI,dipole,d2,z,f,mu);


%% Calculating BFP images
mask=pupil.*exp(1i*aberr+1i*dz*Defocus);
%% Calculating PSF as seen on the camera
% czt2 chirp-z transformation
I_xx=abs(czt2(Ex_Px(:,:).*mask,uk,ux,Nx*oversampling)).^2;
I_yx=abs(czt2(Ey_Px(:,:).*mask,uk,ux,Nx*oversampling)).^2;
psf(:,:)=I_xx+I_yx;

% Calculate dipole excitation probability based on
% angle between electric field and dipole orientation
if exProb==0
    excitationProbability = 1;
else
    excitationProbability = (sin(theta))^2;
end

% sum over block matrices (to return to desired pixelsize)
psf = squeeze(sum(sum(reshape(psf,oversampling,Nx,oversampling,Nx),1),3));

C_norm=sum(sum(psf));  %normalization to total intensity
psf=psf(:,:)/C_norm*n_photon*excitationProbability;
I_xx = I_xx./C_norm*n_photon;
I_yx = I_yx./C_norm*n_photon;

coordinates = [Nx/2+xpos/100, Nx/2+ypos/100];
if strcmp(shotnoise,'y')
    psf = poissrnd(psf+bgnoise);
else
    psf = psf + bgnoise;
end

%% Plotting total PSF
if nargin==0
    figure
    imagesc(psf(:,:),[0 max(psf(:))]); colorbar; axis equal; axis tight; title(['Pixelsize = ', num2str(ux*oversampling*10^9), 'nm. ', 'Dipole = [', num2str(dipole(1)),',',num2str(dipole(2)) ,']']);
end
end