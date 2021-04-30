function [fft_u, fft_psi ] = SQG_large_UQ (model, fft_b)
% Compute the streamfunction and the velocity from the
% buoyancy according to usual SQG if a_H = 0 i.e. k_c = inf


%% Fourier grid
kx=model.grid.k.kx;
ky=model.grid.k.ky;
on_k=model.grid.k.on_k;

%% Streamfunction
fft_psi=bsxfun(@times,on_k/model.physical_constant.buoyancy_freq_N,fft_b);

%% Velocity
fft_u(:,:,1,:) = bsxfun( @times, 1i*(-ky) , fft_psi ) ;
fft_u(:,:,2,:) = bsxfun( @times, 1i*(kx) , fft_psi ) ;
