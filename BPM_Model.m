%%
%addpath('/projectnb/npbvan/ns/WavefrontShaping/PhaseMaskRepo/Generate_TM/')

% clc
% clear
 rng(jobid);


g=[0.8];
factor_x_all=[0.2];

ls_all=[20:20:100,150:50:500]/5;

 
 for ii=1:length(factor_x_all)
     
     for jj=1:length(ls_all)
         
         factor_x=factor_x_all(ii);
         ls=ls_all(jj);
     




   
psd_type='gaussian';

lambda = 0.5; % um


dx_speckle = lambda/2; 
rho_speckle = lambda/4; % seed density

dx_pixel=rho_speckle; %um
N_obj = [6000,6000];

N_diffuser=70;

N_residue=20;% number of layers used to measure axial resolution


d=2;% layer distance
dz=[ones(1,N_diffuser)]*d;



Prop_mode = 1;

% % Define Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

NA = 0.5;

dx_obj=dx_pixel;
[x,y] = meshgrid([-N_obj(2)/2:N_obj(2)/2-1]*dx_obj,[-N_obj(1)/2:N_obj(1)/2-1]*dx_obj);


 
% NA is too large
%N = 500;


% maximum spatial frequency set by NA
um_m = NA/lambda;
% system resolution based on the NA

max_ph=sqrt(d/ls/10)*pi;

%max_ph=0;



%N_illum = 1;

ph_mask = zeros(N_obj(1),N_obj(2),N_diffuser);
for m = 1:N_diffuser
    
    [ph_mask(:,:,m), g_anisotropy] = RandPhaseScreen_RealSpaceNM(psd_type, dx_speckle*factor_x, ...
        rho_speckle, max_ph, dx_pixel, N_obj, lambda);
% figure; imagesc(ph_mask(:,:,m)); axis image; colormap jet; colorbar;dx_speckle/100
end
% drawnow;

o_slice=exp(1i*ph_mask);


%%



%k2 = pi*lambda*(u.^2+v.^2);




%phi = zeros(N1,N2,Nslice);
%i0=
%% illumination beam

% angle_max = [0.25,0.25];
% 
% angle_v=linspace(-angle_max(1),angle_max(1),N_obj(1));
% angle_h=linspace(-angle_max(2),angle_max(2),N_obj(2));
% 
% [angle_hh, angle_vv] = meshgrid(angle_h, angle_v);
% 
% sin_thetav = sin(angle_vv);
% sin_thetah = sin(angle_hh);
% 
% 
% % corresponding spatial freq for each illumination angle
% v_illum = sin_thetav/lambda;
% u_illum = sin_thetah/lambda;
% 
% 
% 
% i0=exp(1i*2*pi*(u_illum(m)*x+v_illum(m)*y));


k = 2*pi/lambda;
z1=sum(dz(1:N_diffuser-N_residue));%um

N = min(NA*z1/dx_pixel,N_obj(1)); 
% % NA is too large
% N = 500;

% incident plane wave
%i0 = exp(1i*k*sqrt(x.^2+y.^2+z1.^2))./sqrt(x.^2+y.^2+z1.^2);
sourceType='Gaussian';
switch sourceType
    case 'Point'
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
        z2=z1;
        i0 = exp(-(x.^2+y.^2)./waist.^2).*exp(-1i*k*(x.^2+y.^2)./2/(z2+zR^2/z2));
        
end



%N=500;
% 




%%
phi = zeros(N_obj(1),N_obj(2),N_diffuser);
phi(:,:,1) = i0; % incident field of 1st slice is illumination
% initialize output field at each slice
psi = zeros(N_obj(1),N_obj(2),N_diffuser);
psi(:,:,1) = i0;

% psi(:,:,1)=zeros(N_obj(1),N_obj(2)); % point source
% psi(round(N_obj(1)/2),round(N_obj(2)/2),1)=1;
um_obj= 1/dx_obj /2;

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
     
    
    H = exp(1i*(dz(m-1))*real(sqrt(((2*pi)^2/lambda^2-k2).*eva)));
    
    
    % propagate from neiboring slices
    phi(:,:,m) = Ft((F(psi(:,:,m-1))).*H); % propagation term
    % output field = incidence * object
    psi(:,:,m) = phi(:,:,m).*o_slice(:,:,m-1); % phase mask operation
    
    %imagesc(abs(psi(:,:,m)).^2/max(max(abs(psi(:,:,m)).^2)));pause(1)
    
end


outputWavefront=squeeze(psi(:,:,end-N_residue));
outputXZ=squeeze(psi(:,round(end/2),:));
% 
% figure
% imagesc((abs(psi(500:end-500,500:end-500,m)).^2/max(max(abs(psi(500:end-500,500:end-500,m)).^2))));
% %%
% figure
% imagesc(abs(psi(:,:,m)).^2)




    
    %Etemp=squeeze((psi(:,:,mm)));
%     Etemp=squeeze((psi(round(end/2),round(end/2),:)));
%     
%     
%     Iall(nn,:)=abs(Etemp(:)).^2;
    
save(['Output_g_',num2str(g(ii)),'_ls_',num2str(ls),'config',num2str(jobid),'.mat'],'outputWavefront','outputXZ')

end


end


% %%
% figure
% I=mean(Iall);
%plot(d*(1:length(I(2:20))),log(I(2:20)));xlabel('z (um)');ylabel('log(I)')
figure;imagesc((abs(i0)));title('incident field')

figure
imagesc(abs(outputWavefront));title('Output field')
%
% %% calculate g
% k = 2*pi/lambda;
% o_slice=exp(1i*ph_mask);
% itemp=squeeze(o_slice(:,:,1));
% % for ii=2:10
% %     itemp=itemp.*squeeze(o_slice(:,:,ii));
% %     
% % end
% 
% [nx,ny]=size(itemp);
% 
% for mm=1:nx
%     itemp2=itemp(mm,:);
%     phase=(angle(itemp2));
%     kx(mm,:)=gradient(phase,dx_pixel);
%     
%     
% end
% kx=kx(2:end-1,2:end-1);
% 
% for mm=1:ny
%     itemp2=itemp(:,mm);
%     phase=(angle(itemp2));
%     ky(:,mm)=gradient(phase,dx_pixel);
% end
% ky=ky(2:end-1,2:end-1);
% grad_ph = sqrt(kx.^2+ky.^2);
% 
% sin_theta=grad_ph/k;
% %sin_theta=min(sin_theta,1);
% cos_theta=real(sqrt(1-sin_theta.^2));
% g = mean2(cos_theta)
% % % calculating kx





% dx_obj=1; % um
% 
% [x,y] = meshgrid([-N_obj(2)/2:N_obj(2)/2-1]*dx_obj,[-N_obj(1)/2:N_obj(1)/2-1]*dx_obj);
% 










