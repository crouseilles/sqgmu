function d_fft_b_adv = deriv_fft_advection(model, fft_b, w)
% Compute the Fourier transform of the partial derivative  
% of the advection term (adv1) and hyperviscosity (adv2) term

%% Grid of wave vectors
PX=model.grid.MX/2;%frequency of aliasing
kx=model.grid.k.kx;%modes en x
ky=model.grid.k.ky;%modes en y
k2=model.grid.k.k2;%|k|^2

%% Advection term
%%%%%%%%%%%%%%%%%
gradb(:,:,1)=real(ifft2(- 1i * kx .* fft_b));
gradb(:,:,2)=real(ifft2(- 1i * ky .* fft_b));
wgradT=sum(bsxfun(@times,w,gradb),3);
adv1=fft2(wgradT);
adv1(PX(1)+1,:)=0;%Remove aliasing
adv1(:,PX(2)+1)=0;%Remove aliasing

%% Hyperviscosity
%%%%%%%%%%%%%%%%%
k2=model.grid.k_HV.k2;%|k|^2 (different de k2)
adv2 = - model.advection.HV.val * k2 .^ (model.advection.HV.order/2) .* fft_b;

%% Summing terms
d_fft_b_adv=adv1+adv2; 

clear adv1 adv4 gradb wgradT

