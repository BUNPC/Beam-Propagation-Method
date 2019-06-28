
sourceType='Gaussian';
NA = 0.1; % numerical aperture
ls=100; % scattering mean free path
lambda = 0.5; % wavelength um
sigma_x=0.6*lambda/2; % this determines g. to obtain the relation between g and sigma_x, use BPM_Estimating_g

%dx_speckle = lambda/2; %
seed_density = lambda/4; % phase seed density, 
dx_pixel=seed_density; % pixel size, normally set to be equal to seed density
N_obj = [1024,1024];

N_diffuser=20;% number of layers
d=10;% layer distance
N_residue=0;% number of layers used to measure axial profile beyond the focal plane, set it to be zero if you don't need this 



if d>ls
    
   warning('layer distance d is larger than scattering mean free path ls')
    
end

dz=[ones(1,N_diffuser)]*d;
z1=sum(dz(1:N_diffuser-N_residue));%um

if NA*z1/dx_pixel>N_obj(1)
    warning('The area of the input wavefront is not large enough to cover the desired NA. Try decrease NA, or increase dx_pixel, N_obj.')
    
elseif NA*z1/dx_pixel>N_obj(2)
    warning('The area of the input wavefront is not large enough to cover the desired NA. Try decrease NA, or increase dx_pixel, N_obj.')
    
end




