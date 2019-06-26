%%
clc
clear

ModelSettings;

d=ls;
Ndiffuser=1;

dz=[ones(1,N_diffuser)]*d;

% % Define Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

[x,y] = meshgrid([-N_obj(2)/2:N_obj(2)/2-1]*dx_pixel,[-N_obj(1)/2:N_obj(1)/2-1]*dx_pixel);


 
% NA is too large
%N = 500;


% maximum spatial frequency set by NA
um_m = NA/lambda;
% system resolution based on the NA

sigma_p=sqrt(d/ls/10)*pi;
%max_ph=pi/5;
%max_ph=0;



%N_illum = 1;

ph_mask = zeros(N_obj(1),N_obj(2),N_diffuser);
for m = 1:N_diffuser
    
    [ph_mask(:,:,m)] = RandPhaseScreen_RealSpaceNM(sigma_x, ...
        seed_density, sigma_p, dx_pixel, N_obj, lambda);
% figure; imagesc(ph_mask(:,:,m)); axis image; colormap jet; colorbar;dx_speckle/100
end
% drawnow;

o_slice=exp(1i*ph_mask);


k = 2*pi/lambda;
z1=sum(dz(1:N_diffuser-N_residue));%um

N = min(NA*z1/dx_pixel,N_obj(1)); 
% % NA is too large
% N = 500;

% incident plane wave
%i0 = exp(1i*k*sqrt(x.^2+y.^2+z1.^2))./sqrt(x.^2+y.^2+z1.^2);
sourceType='Plane';
switch sourceType
    case 'Plane'
i0=zeros(N_obj(1),N_obj(2));
i0(:,:)=1;
    case 'Pencil'
waist = 3;% um
i0 = exp(-(x.^2+y.^2)./waist.^2);
    case 'Spherical'
        i0 = exp(-1i*k*sqrt(x.^2+y.^2+z1.^2))./sqrt(x.^2+y.^2+z1.^2);
        [xx,yy] = meshgrid([-N_obj(2)/2:N_obj(2)/2-1],[-N_obj(1)/2:N_obj(1)/2-1]);
        mask = zeros(N_obj(2),N_obj(1));
         mask(xx.^2+yy.^2 < N^2)=1;
         i0=i0.*mask;            
    case'Gaussian'
        waist0 = lambda/NA/pi;% um
        zR=waist0^2/lambda;
        waist=waist0*sqrt(1+z1^2/zR^2);
        z2=0;
        i0 = exp(-(x.^2+y.^2)./waist.^2).*exp(-1i*k*(x.^2+y.^2)./2/(z2+zR^2/z2));
        
end



phi = zeros(N_obj(1),N_obj(2),N_diffuser);
phi(:,:,1) = i0; % incident field of 1st slice is illumination
% initialize output field at each slice
psi = zeros(N_obj(1),N_obj(2),N_diffuser);
psi(:,:,1) = i0;

% psi(:,:,1)=zeros(N_obj(1),N_obj(2)); % point source
% psi(round(N_obj(1)/2),round(N_obj(2)/2),1)=1;
um_obj= 1/dx_pixel/2;

Xsize = (N_obj(2)-1)*dx_pixel;% total span of x
du = 2*pi/(Xsize);
umax = pi/(dx_pixel);
u = -umax:du:umax;
[U,V]=meshgrid(u,u);
% propogator
%k2 = pi*lambda*(U.^2+V.^2);
k2 = (U.^2+V.^2);

eva = double(k2<(2*pi/lambda)^2);  


for m = 2:(N_diffuser+1)
    
    %H = exp(1i*2*pi*(-z2)*real(sqrt((1/lambda^2-k2/pi/lambda).*eva)));
     
    
    H = exp(1i*(dz(m-1))*real(sqrt(((2*pi)^2/lambda^2-k2).*eva))*1.33);
    
    
    % propagate from neiboring slices
    phi(:,:,m) = Ft((F(psi(:,:,m-1))).*H); % propagation term
    % output field = incidence * object
    psi(:,:,m) = phi(:,:,m).*o_slice(:,:,m-1); % phase mask operation
    
    %imagesc(abs(psi(:,:,m)).^2/max(max(abs(psi(:,:,m)).^2)));pause(1)
    
end


outputWavefront=squeeze(psi(:,:,end));


%
% %% calculate g
k = 2*pi/lambda;
o_slice=exp(1i*ph_mask);
itemp=squeeze(o_slice(:,:,1));
% for ii=2:10
%     itemp=itemp.*squeeze(o_slice(:,:,ii));
%     
% end
% 
 [nx,ny]=size(itemp);

for mm=1:nx
    itemp2=itemp(mm,:);
    phase=(angle(itemp2));
    kx(mm,:)=gradient(phase,dx_pixel);
    
    
end
kx=kx(2:end-1,2:end-1);

for mm=1:ny
    itemp2=itemp(:,mm);
    phase=(angle(itemp2));
    ky(:,mm)=gradient(phase,dx_pixel);
end
ky=ky(2:end-1,2:end-1);
grad_ph = sqrt(kx.^2+ky.^2);

sin_theta=grad_ph/k;
%sin_theta=min(sin_theta,1);
cos_theta=real(sqrt(1-sin_theta.^2));
g = mean2(cos_theta);


clear d ls

fprintf('The anisotropy factor g is %.2f \n',g)








