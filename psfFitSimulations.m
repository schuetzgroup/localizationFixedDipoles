function [xerror,yerror,derror] = psfFitSimulations(n,theta,input,derror_est)
% repeats PSF simulation + localization N times (with different realization
% of noise in each simulation) and compute mean deviation and std from true
% solution
% n .. photon count
% theta .. inclination angle
phi = input(1); % azimuthal angle
noise = input(2); % background noise (mean)
astigmatism = input(3); % astigmatism coefficient
UNC = input(4); % angle uncertainty
reducedExcitation = input(5); %1 or 0

N = 1; %1000; % number of simulations per data point
Nx = 17; % simulation field size (pick odd number!)

pixelsize = 108; % pixelsize (+BFP resolution) also needs to be changed in psfFunction.m and PSF.m
BFP_resolution = 129; % pick odd number

def_grid = 11; % amount of defocus simulation points
xerror = zeros(N,def_grid);
yerror = zeros(N,def_grid);
derror = zeros(N,def_grid);


%% calculate Cramer Rao Bound
crb = [];
for def=-500:100:500
    psf=PSF(n, theta,phi,50,50,[noise+0.1,def, astigmatism, reducedExcitation,Nx],'n');
    [a,b,~]=CRB_defocus(psf, pixelsize, n,theta,phi,noise,def,astigmatism, reducedExcitation);
    crb = [crb; a,b];
end
crb = sqrt(crb);
printData(crb(:,1)) % print CRB


%% Fit simulations
for defocus = 0:1000/(def_grid-1):1000
    index = defocus/(1000/(def_grid-1))+1;
    for k = 1:N
        % simulate estimation of noise from background measurement
        noise_estimate = mean(poissrnd(noise*ones(51,51)),'all');
        % displacement from center in %pixel
        xpos = -100+200*rand();
        ypos = -100+200*rand();
        
        % create PSF
        psf  = PSF(n,theta,phi,xpos,ypos,[noise,-500+defocus, astigmatism, reducedExcitation,Nx],'y');
        
        % add some noise to dipole angles if selected
        if UNC==0 % perfect estimation of orientation
            theta_noise = theta;
            phi_noise = phi;
        else % normal distribution with std 2/2°
            theta_noise = theta+ 2*randn()*pi/180;
            phi_noise = phi+ 2*randn()*pi/180;
        end
        
        % fit position and defocus
        if nargin == 3
            % calculate fit
            x = psfFitMLE(psf,noise_estimate, [theta_noise,phi_noise,astigmatism]);
            % calculate error
            xerror(k,index) = x(1) - xpos;
            yerror(k,index) = x(2) - ypos;
            derror(k,index) = x(3) + 500 - defocus;
        end
        
        % fit position only (if defocus estimate is provided)
        if nargin == 4
            x = psfFitMLE(psf,noise_estimate,[theta_noise,phi_noise,astigmatism,derror_est(index)-500+defocus]);
            xerror(k,index) = x(1) - xpos;
            yerror(k,index) = x(2) - ypos;
        end
    end
end

% account for pixel size
xerror = pixelsize/100*xerror;
yerror = pixelsize/100*yerror;


disp('-------------------')
disp('xprecision')
disp(num2str(std(xerror)))
disp('xbias')
disp(num2str(mean(xerror)))
disp('-------------------')
disp('yprecision')
disp(num2str(std(yerror)))
disp('ybias')
disp(num2str(mean(yerror)))
disp('====================================================')

end