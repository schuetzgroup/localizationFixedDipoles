function psf = psfFunction(x,n_photon,noise, input, stuff)
% Same as PSF.m but with some precomputed parts being passed by psfFitMLE

oversampling = 3; % oversampling factor in fitting procedure

xpos = x(1);
ypos = x(2);

if length(input) == 3 % defocus is also fitted
    astigmatism = input(3);
    def = x(3);
end

if length(input) == 4 % defocus provided as input
    astigmatism = input(3);
    def = input(4);
end

% precomputed variables
Ex_Px = stuff{1};
Ey_Px = stuff{2};
pupil = stuff{3};
Defocus = stuff{4};
uk = stuff{5};
ux = stuff{6}/oversampling;
Nx = stuff{7};
NA = stuff{8};
N=129;


%% Considering aberrations
dz=def*1e-9;  %defocus of objective lens (negative values mean moving the focus into the fluid)
tilt = ux*uk*N/4;
Z_aberr2 = zeros(5,1);
Z_aberr2(1) = -1/100*(xpos)*oversampling*tilt;
Z_aberr2(2) = -1/100*(ypos)*oversampling*tilt;
Z_aberr2(5) = astigmatism;
aberr = sum(ZernikeCalc([2 3 6], [Z_aberr2(1);Z_aberr2(2);Z_aberr2(5)], pupil,'NOLL'),3);


%% Calculating DFT
mask=pupil.*exp(1i*aberr+1i*dz*Defocus);
% czt2 chirp-z transformation
I_xx=abs(czt2(Ex_Px(:,:).*mask,uk,ux,Nx*oversampling)).^2;
I_yx=abs(czt2(Ey_Px(:,:).*mask,uk,ux,Nx*oversampling)).^2;
psf(:,:)=I_xx+I_yx;

% sum over block matrices (to return to desired pixelsize)
if oversampling>1
    psf = squeeze(sum(sum(reshape(psf,oversampling,Nx,oversampling,Nx),1),3));
end

C_norm=sum(sum(psf));  %normalization to total intensity
psf=psf./C_norm*n_photon + noise;

pnorm = norm(psf,2);
psf = psf./pnorm;

end