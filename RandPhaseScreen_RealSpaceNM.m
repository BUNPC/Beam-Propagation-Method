function [ph_mask, g_anisotropy] = RandPhaseScreen_RealSpaceNM(psd_type, dx_speckle, ...
    rho_speckle, max_ph, dx_pixel, N, lambda)
% RandPhaseScreen 
% simulating forward scattering using multi-slice model, assuming each
% slice is a thin phase mask.
% the profile of each phase mask is related to physical properties of
% biological samples according to the refs.
% ref: Schott, Bertolotti, Leger, Bourdieu, Gigan, OE, 2015

% Outputs:
%   ph_mask: simulated phase screen mask that satisfy required statistical
%   distribution
%   g_anisotropy= [g_sim, g_theory1, g_theory2]: estimated anisotropy 
%   parameters based three different methods:
        % g_sim: based on definition of g and direct calculation
        % g_theory1,2: two kind of theory predictions
%
% Inputs:
%   psd_type: power spectral density (PSD) type of the random phase distribution
%   dx_speckle: speckle size as defined by the std of gaussian bump % can
%   change g using this
%   rho_speckle: seeding density of random phase bumps
%   max_ph: amplitude of each rand phase bump
%   dx_pixel: pixel size'
%   N = [N1, N2]: size of the phase screen
%   lambda: wavelength
%   Normalization after spatial averaging by Xiaojun 03/06/2019

%% general parameters
% Define Fourier operators
% F = @(x) fftshift(fft2(ifftshift(x)));
% Ft = @(x) fftshift(ifft2(ifftshift(x)));


F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
% psd_type = 'gaussian';

% % focal length f calculated from tube length and magnification
% mag = 20;
% NA = 0.4;
% z_tube = 200e3;
% f = z_tube/mag;
% lambda = 0.5;
k = 2*pi/lambda;
% N = 1024*2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: each phase mask is created by a random matrix by multiplying its
% Fourier transform with a 2D Gaussian function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define the coordinates
% spatial sampling should ensure enough samples per PSF
% dx_speckle = lambda*10; % half size of PSF/speckle size/ std of gaussian bump
% dx_pixel = dx_speckle/2; % max sampling space by Nyquist

% space
x = [-N(2)/2:N(2)/2-1]*dx_pixel;
y = [-N(1)/2:N(1)/2-1]*dx_pixel;
[xx,yy] = meshgrid(x,y);

Nsame=round(rho_speckle/dx_pixel);

ph_mask=normrnd(0,max_ph,[length(x),length(y)]);
if (Nsame>1)
   for ii=1:Nsame:N(1)-Nsame
       for jj=1:Nsame:N(2)-Nsame
      ph_mask(ii+1:ii+Nsame,jj+1:jj+Nsame) =ph_mask(ii,jj);
      
       end
   end
       
end
    
    

% r2 = xx.^2+yy.^2;

% umax = 1/dx_pixel/2;
% du = 1/dx_pixel/N(2);
% dv = 1/dx_pixel/N(1);
% u = [-umax:du:umax-du];
% v = [-umax:dv:umax-dv];
% [uu,vv] = meshgrid(u,v);
% ww2 = uu.^2+vv.^2;

% % evanescent field mask
%eva = double(sqrt(uu.^2+vv.^2)<1/lambda);

%% assume uniform distribution of seedings
% define density parameter
% rho_speckle = 1;
% the index of ON pixels
% idx_on = randsample(prod(N),round(prod(N)*rho_speckle));
% mask_on = zeros(N);
% mask_on(idx_on) = 1;
% tmp = (rand(N)-0.5)/0.5/(dx_speckle/dx_pixel).*mask_on;
% TMP = F(tmp);
% 
% % number of speckles
% Nspeckle = prod(N)*(dx_pixel/dx_speckle)^2*rho_speckle;
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



switch psd_type
  case 'gaussian'
        
        sigma_x = dx_speckle; 
        sigma_psd = 1/(2*pi*sigma_x);
        % parameter 2: large ampl means large gradient of phase -> smaller g
%         max_ph = pi*10;
        % only the shape factor (max_ph/sigma_x) defines the phase gradient and
        % thus g
        
        % unit in F is dx^2
        % then unit in Ft is (du*N)^2; note N comes from definition of IDFT
        % in matlab
        % another way to say this is that each pixel should accumulate all
        % the phase in the bins as the total phase at that pixel
        ph_kernel = exp(-(xx.^2+yy.^2)/2/sigma_x^2)/(2*pi*sigma_x^2);
        ph_kernel=ph_kernel/sum(ph_kernel(:));

        %PH_psd = abs(F(ph_kernel)*dx_pixel^2).^2;
%         PH_psd2 = (max_ph*N^2/(2*pi*(sigma_psd/du)^2))^2*...
%             exp(-2*rho2/(2*sigma_psd^2));
%         PH_psd = (max_ph/(2*pi*(sigma_psd)^2))^2*...
%             exp(-2*rho2/(2*sigma_psd^2));
        
%     case 'power'
%         % referece: Martin, Flatte, Appl. Opt. 1988 27
%         % parameter: alpha 5/3 for Kolmogorov
%         %            Cn: index strcture constant, relates the 'strength'
%         %            parameter
%         C_factor = 2*gamma(1-alpha/2)*gamma(1+alpha)/alpha/(1+alpha/2)*...
%             cos(pi*alpha/4)*sin((alpha-1)*pi/2);
%         % strength parameter
%         sigma2 = C_factor*Cn^2*k^(2-alpha/2)*R^(1+alpha/2);
%         % required spatial bandwidth (fmax/fmin)
%         R_rho = 2*(sigma2)^(2/alpha)*10^(1/(alpha+2))*...
%             ((1+alpha/2)/(2*alpha*gamma(1+alpha/2)*cos(pi*alpha/4)))^(2/alpha);
%         
%         
%         K_factor = gamma(alpha+1)/4/pi^2*sin((alpha-1)*pi/2);
%         PH_psd = K_factor*Cn^2*ww2.^((-alpha-2)/2);
%     case 'inner'
%         % parameter: alpha: 5/3 for Kolmogorov
%         %            Cn: index strcture constant, relates the 'strength'
%         %            parameter
%         %            beta
%         PH_psd = K_factor*Cn^2*ww2.^((-alpha-2)/2).*exp(-ww2*beta^2/4);
%     case 'out'
%         % parameter: L0: outer scale
%         PH_psd = K_factor*Cn^2*(ww2+1/L0^2).^((-alpha-2)/2);


% %ph_mask = conv2(exp(1i*ph_mask),ph_kernel,'same');

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
% double check


 ph_mean1=mean(ph_mask(:));
 ph_std1=std(ph_mask(:));
% cut off evanecent wave in ph_mask
 ph_mask=real((Ft((F(ph_kernel).*F(ph_mask).*eva))));
 
 ph_mean2=mean(ph_mask(:));
 ph_std2=std(ph_mask(:));
% 


ph_mask=(ph_mask-ph_mean2)*ph_std1/ph_std2+ph_mean1;


% 
itemp=exp(1i*ph_mask);

% ph_kernel = real(Ft(sqrt(PH_psd)));
% figure; plot(x,ph_kernel(N/2+1,:))
% problem here:
% 1) when changing real space pixel size, should only mean the sampling is
% changed, but ph also changes, which is not physical
% 2) the factor 1/(sigma_x/dx) heuristically corrects for the scaling, now
% g and dist of ph/grad_ph does not depends on dx
% oversampling factor dx0/dx is divded out
% ph = real(Ft(F(tmp)*dx^2.*eva.*sqrt(PH_psd)))/(dx0/dx);

% Fourier transform of the seed should satisfy:
% real part satisfy Rayleigh distribution
% phase takes on uniform distribution in [-pi, pi]
% the N factor come from the fact that we are NOT taking an average
% (whereas in Goodman's calculation, the average of N random phasors are
% taken)
% (dx0/dx) factor comes from the consideration that the number of speckle 
%  should be determined by the resolution of the imaging system!
% physical meaning: Fourier spectrum of iid Nspeckle number of speckles.
% within each speckle, they are correlated.
% SEED = raylrnd(sqrt(Nspeckle/6),N(1),N(2)).*exp(1i*(rand(N)-0.5)/0.5*pi);
% ph = real(Ft(SEED.*eva.*sqrt(PH_psd))*(N*du)^2);
% figure; imagesc(ph); colorbar;
% note that ph is complex and both satisfy the statistical requirement
% sqrt(2) might come from this fact, loss of 1/2 energy is not accounted in
% later pure real space modeling.

% figure; hist(ph(:),100);
%% predicting the statistics of phase
% figure; hist(ph(:),100)
% sigma_ph = sqrt(pi/3)*max_ph*sigma_x/dx_speckle*sqrt(rho_speckle);
% % std(ph(:))/sigma_ph
% 
% % figure; imagesc(ph_mask); colorbar; axis image;
% [counts, phs] = hist(ph_mask(:),100);
% prob = counts/sum(counts);
% % figure; plot(phs/pi,prob);
% % xlim([-4*sigma_ph/pi,4*sigma_ph/pi]);
% % xlabel('phase/\pi'); ylabel('probability')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % note:
% % 1) in real space, it is equivalent that a uniform dist random variable is
% % convolved with a gaussian function
% % 2) convolution can be considered as a weighted sum of all the random
% % variables
% % 3) the pdf for sum of random variables is the convolution of pdfs of each
% % random variables
% % 4) the results dist of ph is gaussian
% %    mean = 0, sigma_ph^2 = pi/3*max_ph^2*(sigma_x/dx0)^2
% % physical meaning: var of phase is proportional to max phase^2 and number of
% % speckles within each gaussian envelop
% % note that we need to normalized out the pixel size
% % to get a ph dist independent on # of pixels per sigma area, we can
% % normalize by (sigma_x/dx), this solves the problem that changing dx while
% % keeping sigma_x the same changes the ph dist
% % 5) gradx_ph and grady_ph is also gaussian dist since the only change is
% % the weightings
% % sigma_gx = sigma_gy = sigma_g = sqrt(pi/6)*max_ph/dx0
% 
% % 6) grad_ph is taking the abs of two dimensional gaussian dist random
% % variables, so it is Rayleigh dist
% % f(grad_ph) = grad_ph/sigma_grad^2 * exp(-grad_ph^2/2/sigma_grad^2)
% 
% % 7) calculating g is taking the expecation of cos_theta,
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ph_mask2=unwrap(angle(itemp));

[gradx_ph,grady_ph] = gradient(ph_mask2,dx_pixel);

% sigma_grad = sqrt(pi/6)*max_ph/dx_speckle*sqrt(rho_speckle);
% [counts, grads] = hist(gradx_ph(:),100);
% prob = counts/sum(counts);
% figure; plot(atand(grads/k),prob);
% xlim([-3*atand(sigma_grad/k),3*atand(sigma_grad/k)])


% std(gradx_ph(:))/sigma_grad

grad_ph = sqrt(gradx_ph.^2+grady_ph.^2);
% figure; imagesc(grad_ph); colorbar;
% grad_ph follows Rayleigh distribution, mean = sqrt(pi/2)*sigma_grad
%mean_grad = sigma_grad*sqrt(pi/2);
% mean2(grad_ph(:))/ mean_grad
% figure; hist(grad_ph(:),100)


%%%%
% probability distribution of theta can be calculated accoding to chane
% rule:
% p(theta) =
% tan(theta)/cos(theta)^2*k^2/sigma^2*exp(-tan(theta)^2*k^2/2/sigma^2)
%%%%
%theta_scat = atan(grad_ph/k);

sin_theta=grad_ph/k;
sin_theta=min(sin_theta,1);
cos_theta=abs(sqrt(1-sin_theta));
% check if statistics matches with predictions
% [counts, thetas] = hist(theta_scat(:),100);
% prob = counts/sum(counts);
% figure; plot(thetas/pi*180,prob,'--b');
% hold on;
% p_theta = tan(thetas)./cos(thetas).^2*k^2/sigma_grad^2.*...
%     exp(-tan(thetas).^2*k^2/2/sigma_grad^2)*(thetas(2)-thetas(1));
% % plot(thetas/pi*180,p_theta,'r');
% xlabel('angle (degree)');
% ylabel('probability');
% legend('simulated','theoretical')
% hold off;



%%
%cos_theta = cos(theta_scat);
% figure; hist(cos_theta(:),100)
% estimate anisotropic parameter g
g_sim = mean2(cos_theta);

%%%%%%%%%%%%%%
% the derivation shows that for a gaussian phase screen. Given the desired
% speckle size as defined by the imaging system. The Gaussian blurring
% kernel size is no longer relevant!
% The parameter that controls the g-factor is 
% sigma0 = sqrt(pi/6)*ph_max/dx0*lambda/2/pi
% dx0 is the speckle size
% g = sqrt(pi/2)/sigma0*exp(1/2/sigma0^2)*(1-erf(1/sqrt(2)/sigma0))
%%%%%%%%%%%%%%

% g_theory1= sum(cos(thetas).*p_theta);
% 
% %sigma_grad = sqrt(pi/6)*max_ph/dx0;
% sigma_tan = sigma_grad/k;
% g_theory2 = sqrt(pi/2)/sigma_tan*exp(1/2/sigma_tan^2)*erfc(1/sqrt(2)/sigma_tan);
% g_theory2(isnan(g_theory2)) = 1;


% % compare spectrum between phase and object
% sp_ph = F(ph);
% sp_obj = F(exp(1i*ph)-1);

% %% numerical calculation from the single Gaussian bump
% [gradx_gs,grady_gs] = gradient(ph_kernel,dx);
% 
% grad_gs = sqrt(gradx_gs.^2+grady_gs.^2);
% 
% grad_gs0 = max_ph/sigma_x^2*sqrt(r2).*exp(-r2/2/sigma_x^2);
% 
% ths = atan(grad_gs0/k);
% cos_ths = cos(ths);
% g_gs = mean2(cos_ths);

g_anisotropy = [g_sim];
end

% %% calculator for g
% 
% sigma_g = @(phi_max, sigma_s, lambda) sqrt(1/24/pi)*phi_max.*(lambda./sigma_s)*sqrt(rho_speckle);
% g = @(phi_max, sigma_s, lambda)...
%     sqrt(pi/2)./sigma_g(phi_max, sigma_s, lambda).*...
%     (exp(1./2./sigma_g(phi_max, sigma_s, lambda).^2).*...
%     erfc(1/sqrt(2)./sigma_g(phi_max, sigma_s, lambda)));
% 
% 
% phi_max = [0.1:0.01:20*pi];
% sigma_s = [0.1:0.01:20]*lambda;
% [pp,ss] = meshgrid(phi_max,sigma_s);
% gg = g(pp, ss, lambda);
% gg(isnan(gg)) = 1;
% 
% figure; imagesc(phi_max/pi,sigma_s/lambda,gg);
% ylabel('\sigma_{scat}/\lambda');
% xlabel('maximum phase/\pi')
% set(gcf, 'Color', 'w');
% 

