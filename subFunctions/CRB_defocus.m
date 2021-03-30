function [CR_x, CR_y, CR_d]=CRB_defocus(psf,pixelsize,n_photon,theta, phi,noise,def, ast, exProb)
%calculating the Cramer Rao lower bound for x/y/defocus estimation
%output is in nm^2

% INPUT
% psf ..  a 2d PSF matrix
% pixelsize .. pixelsize in nm
% n_photon .. total number of photons in psf
% theta .. inclination angle of dipole
% phi .. azimuthal angle of dipol
% bg .. background level (mean) per pixel
% def .. defocus in nm
% ast .. astigmatism coefficient
% exProb .. reduced excitation (0 or 1)

Nx=length(psf); % Nx .. simulation field size

noise = noise+0.1; % to avoid numerical instabilities

delta_def = 1; % size of difference quotient step (defocus)
delta_x = 1;  % size of difference quotient step (lateral position)

% compute slightly displaced PSFs
psf_def = PSF(n_photon, theta, phi,50,50,[noise,def+delta_def, ast,exProb,Nx],'n');
psf_x = PSF(n_photon, theta, phi,50+delta_x,50,[noise,def, ast,exProb,Nx],'n');
psf_y = PSF(n_photon, theta, phi,50,50+delta_x,[noise,def, ast,exProb,Nx],'n');

% compute differential quotient
Dx = (-psf+psf_x)./delta_x;
Dy = (-psf+psf_y)./delta_x;
Ddef = (-psf+psf_def)./delta_def;

LI1=1:(size(psf,1)-1);
LI2=1:(size(psf,2)-1);

Dx=Dx(LI1,LI2); %bring matrices into same size
Dy=Dy(LI1,LI2);
Ddef=Ddef(LI1,LI2);

mu=psf(LI1,LI2);

%calculating Fisher information matrix entries
FI_xx=squeeze(sum(sum(Dx.^2./mu,1),2));
FI_yy=squeeze(sum(sum(Dy.^2./mu,1),2));
FI_dd=squeeze(sum(sum(Ddef.^2./mu,1),2));

FI_xy=squeeze(sum(sum(Dx.*Dy./mu,1),2));
FI_xd=squeeze(sum(sum(Dx.*Ddef./mu,1),2));
FI_yd=squeeze(sum(sum(Dy.*Ddef./mu,1),2));

FI=zeros(3,3); %Fisher matrix
FI(1,1)=FI_xx;
FI(2,2)=FI_yy;
FI(3,3)=FI_dd;

FI(1,2)=FI_xy;
FI(2,1)=FI_xy;
FI(1,3)=FI_xd;
FI(3,1)=FI_xd;
FI(2,3)=FI_yd;
FI(3,2)=FI_yd;

% calculating CRB
%inverting Fisher matrix and taking diagonal values
FI_inv=(pixelsize/100)^2.*inv(FI);
CR_x=FI_inv(1,1);
CR_y=FI_inv(2,2);
CR_d=FI_inv(3,3);

end