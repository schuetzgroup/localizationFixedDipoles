function [SA_out,Defocus,a] =fun_SA_RImismatch(N,n1,n2,NA,lambda_0,rd)
%outputs spherical aberration arising from a RI-mismatch, for z=1m inside
%the material of RI2. Output can be scaled to the existing z-value
%returned defocus phase is defined for medium n2.
%according to paper Jesacher et al. 2010. Opex
%the "defocus-term" is removed analytically
%inputs:
%n1...refractive index of immersion medium
%n2...refractive index of target material
%NA...numerical aperture (e.g. 1.4)
%lambda_0...vacuum wavelength
%a...coefficient of "Defocus" contained in SA
%rd...0 or 1, if set to 1, defocus is removed

k=2*pi/lambda_0;
[~,~,R,mask]=create_coord(N,2/(N-1),'FFT');

%mean value of spherical defocus
MW=@(RI) 2/3*k*(-(-1 + RI^2/NA^2)^(3/2)+(RI^2/NA^2)^(3/2))*NA;
%defocus function
Def=@(RI) real(k*NA*sqrt(RI^2/NA^2-R.^2)-MW(RI)).*mask; %spherical defocus in medium with refractive index RI (for 1m of defocus)

SA=-(Def(n2)-Def(n1)); %spherical aberration phase (complete, i.e including defocus) (for 1m of defocus), (see paper: Jesacher & Booth, Opex 2010)
Defocus=Def(n2);
a_spher=-1/72*k^2*pi*(72*n2^2-36*NA^2+(32*(-(-1 + n2^2/NA^2)^(3/2) + n2^3/NA^3)*(n1^3-n1^2*sqrt((n1 - NA)*(n1 + NA))+NA^2*sqrt((n1 - NA)*(n1 + NA))))/NA - (32*(n2^3 - n2^2*sqrt((n2 - NA)*(n2 + NA))+NA^2*sqrt((n2 - NA)*(n2 + NA)))^2)/NA^4+(9/NA^2)*(2*(-n1^3*n2 - n1*n2^3 + n1^2*sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA))+sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA))*(n2^2-2*NA^2)) - (n1^2 - n2^2)^2*(log((n1 - n2)^2)-log(n1^2 + n2^2 - 2*(NA^2 + sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA)))))));

def_norm=-((k^2*(16*n2^6 - 24*n2^4*NA^2 + 6*n2^2*NA^4 + NA^6 - 16*n2^5*sqrt(n2^2 - NA^2) + 16*n2^3*NA^2*sqrt(n2^2 - NA^2))*pi)/(18*NA^4));
SA_corr=SA-a_spher*Def(n2)/def_norm;
a=a_spher/def_norm; %coefficient of Defocus contained in SA (without normalization)

if rd==1
    SA_out=SA_corr;
else
    SA_out=SA;
end

end