function x = psfFitMLE(psf, noise, input)
% receives psf matrix + additional parameters and estimates position and defocus using MLE

% OUTPUT
% x ... vector of estimated parameter [xpos, ypos, defocus]
% position coordinates are in units of pixels/100, defocus is in nanometers

% INPUT
% psf .. quadratic matrix containing psf
% noise .. estimate of mean background noise
% input .. [dipole(1), dipole(2), astigmatism coefficient]  or [dipole(1), dipole(2), astigmatism coefficient, defocus in nm]

%%
if noise < 1e-5 % this is to avoid numerical instability of log(psf)
    noise=1e-5;
end

% estimate signal photon count
n_photon = sum(sum(psf-noise));

% normalize psf
pnorm = norm(psf,2);
psf = psf./pnorm;

% fit position and defocus
if length(input)==3
    % bounds for x, y, defocus
    lb = [-500, -500,-2000];
    ub = [500,500,2000];
    % initial guess
    x = [-100+200*rand(),-100+200*rand(),-500+1000*rand()];
end

% fit position only
if length(input)==4
    % bounds for x, y
    lb = [-250, -250];
    ub = [250,250];
    % initial guess
    x = [0,0];
end

% parameters
Nx = length(psf);
ux=108e-9; %115e-9; %resolution in focal space
% Additional parameters
N=129;  % resolution in BFP
lambda_0=680e-9;
NA=0.7; % objective NA
dipole = [input(1), input(2)];
RI=[1.33 1.46 1]; %refractive indices; RI=[RI_specimen, RI_intermed., RI_immoil]
d2=0*1e-3;   %thickness of intermediate layer (layer 2)
z=0e-6;    % distance of dipole from medium 2
f=3e-3;    %focal length of objective
mu=1e-12;  %magnitude of dipole (arbitrary value)


% do calculations that are required in every step
%% Calculating defocus and aberrations
% Defocus function refers to refractive index n2
[~,Defocus,~] =fun_SA_RImismatch(N,RI(3),RI(3),NA,lambda_0,1);
uk=4*pi/lambda_0*NA/N; %unit in pupil space (k-space)
[~,~,~,pupil]=create_coord(N,1,'FFT'); % create pupil

%% Calculating BFP-fields for all dipole orientations
[Ex_Px(:,:),Ey_Px(:,:)]=fun_dipole_imaging(N,lambda_0,NA,RI,dipole,d2,z,f,mu); %x-dipole

stuff = {Ex_Px, Ey_Px, pupil, Defocus, uk, ux,Nx,NA};

%% create function handle required by lsqcurvefit
fun = @(x_,xdata_)psfFunction(x_,n_photon,noise,input,stuff); %+ 10^(-2).*norm(x_)*eye(size(psf));
xdata = zeros(Nx,Nx);
options = optimoptions('lsqcurvefit','Algorithm', 'trust-region-reflective', 'OptimalityTolerance', 5e-7, 'Display','off');
x = lsqcurvefit(fun,x,xdata,psf,lb,ub,options);


%% create log pdf function handle required by mle
lnpdf = @(z,x_) sum(z.*log(psfFunction(x_,n_photon,noise, input, stuff))-psfFunction(x_,n_photon,noise,input,stuff)-log(gamma(z+1)) , 'all');
options = optimoptions(@fminunc, 'Display','off','StepTolerance',1e-10, 'OptimalityTolerance', 1e-10);
x = fminunc(@(x_) -lnpdf(psf,x_), x,options);

end