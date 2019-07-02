%%

%d=20;% layer distance 100
dz=[ones(1,N_diffuser)]*d;



%Prop_mode = 1;

% % Define Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));



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
    
    [ph_mask(:,:,m)] = RandPhaseScreen_RealSpaceNM(sigma_x, ...
        seed_density, max_ph, dx_pixel, N_obj, lambda);

end

o_slice=exp(1i*ph_mask);

k = 2*pi/lambda;
z1=sum(dz(1:N_diffuser-N_residue));%um

N = min([NA*z1/dx_pixel,N_obj(1),N_obj(2)]); 
% % NA is too large
% N = 500;

% incident plane wave
%i0 = exp(1i*k*sqrt(x.^2+y.^2+z1.^2))./sqrt(x.^2+y.^2+z1.^2);

switch sourceType
    case 'Plane'
i0=zeros(N_obj(1),N_obj(2));
i0(:,:)=1;
    case 'Spherical'
        i0 = exp(-1i*k*sqrt(x.^2+y.^2+z1.^2))./sqrt(x.^2+y.^2+z1.^2);
        [xx,yy] = meshgrid([-N_obj(2)/2:N_obj(2)/2-1],[-N_obj(1)/2:N_obj(1)/2-1]);
        mask = zeros(N_obj(2),N_obj(1));
         mask(xx.^2+yy.^2 < N^2)=1;
         i0=i0.*mask;            
    case'Gaussian'
        waist0 = lambda/NA/pi;% um
        zR=waist0^2/lambda;
        z2=z1; % location of the focus
        waist=waist0*sqrt(1+z2^2/zR^2);
        i0 = exp(-(x.^2+y.^2)./waist.^2).*exp(-1i*k*(x.^2+y.^2)./(2*(z2+zR^2/z2)));
        
        
end


%%
phi = zeros(N_obj(1),N_obj(2),N_diffuser);
phi(:,:,1) = i0; % incident field of 1st slice is illumination
% initialize output field at each slice
psi = zeros(N_obj(1),N_obj(2),N_diffuser);
psi(:,:,1) = i0;

% psi(:,:,1)=zeros(N_obj(1),N_obj(2)); % point source
% psi(round(N_obj(1)/2),round(N_obj(2)/2),1)=1;
% um_obj= 1/dx_obj/2;
% 
% Xsize = (N_obj(2)-1)*dx_pixel;% total span of x
% du = 2*pi/(Xsize);
% umax = pi/(dx_pixel);
% u = -umax:du:umax; 

% **********  should be ********
% Xsize = N_obj(2)*dx_pixel;% total span of x
% du = 2*pi/(Xsize);
% umax = pi/(dx_pixel);
% u = -umax:du:umax-du; 


% **************  or new  ***********
x=-N_obj(1)/2:N_obj(1)/2-1;
y=-N_obj(2)/2:N_obj(2)/2-1;
LX = N_obj(1)*dx_pixel;
LY = N_obj(2)*dx_pixel;
u=lambda*x/LX;
v=lambda*y/LY;
[uu,vv] = meshgrid(u,v);
% propogator
%k2 = pi*lambda*(U.^2+V.^2);
k2 = (uu.^2+vv.^2);
eva = double(k2/lambda^2<(1/lambda)^2);  

for m = 2:(N_diffuser+1)
    
    H =exp(1i*k*dz(m-1)*sqrt((1-uu.^2-vv.^2).*eva));     
    
   % H = exp(1i*(dz(m-1))*real(sqrt(((2*pi)^2/lambda^2-k2).*eva)));
    % second:   H = exp(1i*k*(dz(m-1))*real(sqrt(1-k2).*eva));
   
        
    % propagate from neiboring slices
    phi(:,:,m) = Ft((F(psi(:,:,m-1))).*H); % propagation term
    % output field = incidence * object
    psi(:,:,m) = phi(:,:,m).*o_slice(:,:,m-1); % phase mask operation
    
    %imagesc(abs(psi(:,:,m)).^2/max(max(abs(psi(:,:,m)).^2)));pause(1)
    
end


outputWavefront=squeeze(psi(:,:,end-N_residue));
outputXZ=squeeze(psi(:,round(end/2)+1,:));
    
%save(['Output_g_',num2str(g(ii)),'_ls_',num2str(ls),'config',num2str(jobid),'.mat'],'outputWavefront')


figure; imagesc(abs(squeeze(psi(:,round(end/2),:)))); colormap hot; title('amplitude of XZ profile')
figure;imagesc(x,y,(abs(i0)));title('amplitude of incident field'); xlabel('x (\mum)');ylabel('y (\mum)');


figure
imagesc(x,y,abs(outputWavefront));title('amplitude of output field');xlabel('x (\mum)');ylabel('y (\mum)');









