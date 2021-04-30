function [fft_buoy_part, model] = fct_fft_advection_sto(model,  fft_buoy_part)
% Advection of buoyancy using SQG or SQG_MU model
%

%% Folder to save plots and files
model.folder.folder_simu = [ 'images/usual_SQG/' model.type_data ];
% Create the folders
fct_create_folder_plots(model)

% Colormap
load('BuYlRd.mat');
model.folder.colormap = BuYlRd; clear BuYlRd

% Version of matlab
vers = version;
year = str2double(vers(end-5:end-2));
subvers = vers(end-1);
model.folder.colormap_freeze = ...
    (  year < 2014 || ( year == 2014 && strcmp(subvers,'a') ) );

%% Grid

% Spatial grid
model.grid.x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
model.grid.y = model.grid.dX(2)*(0:model.grid.MX(2)-1);

% Grid in Fourier space
model = init_grid_k (model);

%% Initialisation of the spatial fields

% Initial large-scale velocity  (NC : solve "Poisson" equation to get w from buoy)
fft_w = SQG_large_UQ(model, fft_buoy_part);
w=real(ifft2(fft_w));

% Create several identical realizations of the intial buoyancy
fft_buoy_part = repmat(fft_buoy_part(:,:,1),[1 1 1 model.advection.N_ech]);

%% Hyperviscosity

% NC : compute grandient of w
[dxUx,dyUx] = gradient(w(:,:,1)',model.grid.x_ref,model.grid.y_ref);
[dxUy,dyUy] = gradient(w(:,:,2)',model.grid.x_ref,model.grid.y_ref);
s=(dxUx+dyUy)/2;
d=(dyUx+dxUy)/2;
lambda = sqrt(s.^2+d.^2);
model.advection.lambda_RMS = sqrt(mean(lambda(:).^2));
clear s d dxUx dxUy dyUx dyUy lambda

% Order of the hyperviscosity
model.advection.HV.order=8;

% Hyperviscosity coefficient
model.advection.HV.val= ...
    40 * model.advection.lambda_RMS * ...
    (mean(model.grid.dX)/pi)^model.advection.HV.order;

%% Choice of time step : CFL
% CFL of the (large-scale) advection
dX=permute(model.grid.dX,[1 3 2]);
bound1=sum(bsxfun(@times,abs(w),pi./dX),3);
bound1=max(bound1(:));
bound1=1/bound1/4;

% CFL of the hyperviscosity
dX2=(model.grid.dX /pi).^2;
bound2=1/model.advection.HV.val*(prod(dX2)/sum(dX2)) ^ ...
    (model.advection.HV.order/2);
clear dX dX2

% Minimum of the CFL
dt = min([bound1 bound2]); 
clear  bound1 bound2

model.advection.dt_adv = dt;

sigma_on_sq_dt = 0;

%% Loop on time

% Used for the first plot
tt_last = -inf;

% Number of time step
N_t = ceil(model.advection.advection_duration/dt);

% Printing some information
fprintf(['The initial condition is ' model.type_data ' \n'])
fprintf(['1/k_c is equal to ' num2str(1/model.sigma.k_c) ' m \n'])
fprintf(['Time step : ' num2str(dt) ' seconds \n']);
fprintf(['Time of advection : ' num2str(N_t*dt/3600/24) ' days \n']);

for t=1:N_t
    %% Time-uncorrelated velocity (isotropic and homogeneous in space)
    sigma_dBt_dt = 0;    
    
    %% Time-correlated velocity : solve "Poisson" to get w from buoy 
    fft_w = SQG_large_UQ(model, fft_buoy_part);
    w=real(ifft2(fft_w));
    clear fft_w
        
    %% Transport of tracer
    % Runge-Kutta 4 scheme
    fft_buoy_part = RK4_fft_advection(model,fft_buoy_part, w);
    
    %% Plots and save
    tt = floor(t*dt/(3600*24)); % Number of days
    if tt > tt_last
        tt_last = tt;
        fprintf([ num2str(t*dt/(24*3600)) ' days of advection \n'])
        day = num2str(floor(t*dt/24/3600));
        
        % Plots
        fct_plot(model,fft_buoy_part,day)
        
        % Save files
        %save( [model.folder.folder_simu '/files/' day '.mat'], ...
        %    'model','t','fft_buoy_part','w','sigma_dBt_dt', ...
        %    'sigma_on_sq_dt');
    end
end