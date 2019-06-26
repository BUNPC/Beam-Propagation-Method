function [ph_mask] = RandPhaseScreen_RealSpaceNM(sigma_x, ...
    seed_density, sigma_p, dx_pixel, N, lambda)
% RandPhaseScreen 
% simulating forward scattering using multi-slice model, assuming each
% slice is a thin phase mask.
% the profile of each phase mask is related to physical properties of
% biological samples according to the refs.

% Outputs:
%   ph_mask: simulated phase screen mask that satisfy required statistical
%   distribution

%
% Inputs:
%   sigma_x: speckle size as defined by the std of gaussian bump % can
%   change g using this
%   seed_density: seeding density size of random phase bumps, normally just
%   set equal to the pixel size
%   sigma_p: amplitude of each rand phase bump
%   dx_pixel: pixel size'
%   N = [N1, N2]: size of the phase screen
%   lambda: wavelength
%   Normalization after spatial averaging by Xiaojun 03/06/2019

%% general parameters
% Define Fourier operators

F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
% psd_type = 'gaussian';
k = 2*pi/lambda;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: each phase mask is created by a random matrix by multiplying its
% Fourier transform with a 2D Gaussian function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define the coordinates
% spatial sampling should ensure enough samples per PSF
% sigma_x = lambda*10; % half size of PSF/speckle size/ std of gaussian bump
% dx_pixel = sigma_x/2; % max sampling space by Nyquist

% space
x = [-N(2)/2:N(2)/2-1]*dx_pixel;
y = [-N(1)/2:N(1)/2-1]*dx_pixel;
[xx,yy] = meshgrid(x,y);

Nsame=round(seed_density/dx_pixel);

ph_mask=normrnd(0,sigma_p,[length(x),length(y)]);
if (Nsame>1)
   for ii=1:Nsame:N(1)-Nsame
       for jj=1:Nsame:N(2)-Nsame
      ph_mask(ii+1:ii+Nsame,jj+1:jj+Nsame) =ph_mask(ii,jj);
      
       end
   end
       
end
    
ph_kernel = exp(-(xx.^2+yy.^2)/2/sigma_x^2)/(2*pi*sigma_x^2);
ph_kernel=ph_kernel/sum(ph_kernel(:));   
%% assume uniform distribution of seedings
% define density parameter
% seed_density = 1;
% the index of ON pixels
% idx_on = randsample(prod(N),round(prod(N)*seed_density));
% mask_on = zeros(N);
% mask_on(idx_on) = 1;
% tmp = (rand(N)-0.5)/0.5/(sigma_x/dx_pixel).*mask_on;
% TMP = F(tmp);
% 
% % number of speckles
% Nspeckle = prod(N)*(dx_pixel/sigma_x)^2*seed_density;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% converting it in the fourier space
% % some consequence of doing this:
% % ref: Goodman, Speckle phenomena in optics
% % 1. Taking the fourier transform of the seed picture is equivalent to look
% % at the phasor sum of random complex variables. So,
% % 1.1) the Real/Imag part satisfy Gaussian dist:
% %      p(R,I) = 1/2/pi/sigma^2*exp(-(R^2+I^2)/2/sigma^2)
% %      sigma^2 = sum(E(a_n^2)/2)
% % 1.2) the abs A satisfies  Rayleigh dist:
% %      p(A) = A/sigma^2*exp(-A^2/2/sigma^2)
% % 1.3) the ph phi satisfies uniform dist in [-pi, pi]
% % 1.4 abs^2 I satisfies neg-exp dist
% %      p(I) = (1/I0)*exp(-I/I0)
% %      mean I0 = 2*sigma^2
% % for uniform distribution a_n in [a,b]
% %   E(a_n^2) = (a^2+ab+b^2)/3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % note:
% %  1. when applying Fourier-scaling here, only scaling factor in terms of
% %  pixels (sigma_psd/du) matters as FFT does not have units
% %  2. N^2 factor comes from def in FFT

%%%%%%%%%%%%%%%%%%%%
dx_obj=dx_pixel;
um_obj= 1/dx_obj /2;

Xsize = N(2)*dx_pixel;
du = 1/(Xsize);
umax = 1/(2*dx_pixel);
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);
% propogator
k2 = pi*lambda*(U.^2+V.^2);
eva = double(k2<pi/lambda);
%eva=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 ph_mean1=mean(ph_mask(:));
 ph_std1=std(ph_mask(:));
% cut off evanecent wave in ph_mask
 ph_mask=real((Ft((F(ph_kernel).*F(ph_mask).*eva))));
 ph_mean2=mean(ph_mask(:));
 ph_std2=std(ph_mask(:));
% 
ph_mask=(ph_mask-ph_mean2)*ph_std1/ph_std2+ph_mean1;



end



