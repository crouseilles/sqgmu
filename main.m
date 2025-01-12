%%%%%%%%%%%%%%%%%%%%
%%% Main 
%%%%%%%%%%%%%%%%%%%%
init;

%% Main parameters to choose

% Deterministic model
stochastic_simulation = false ;
% Usual SQG model (stochastic_simulation=false)


% Duration of the simulation (in seconds)
advection_duration = 3600*24*20; % 20 days

% Number of realizations in the ensemble
N_ech=1;
% ( N_ech=200 enables moments to converge when the parameter resolution is
%   set to 128 )
% ( N_ech is automatically set to 1 in deterministic simulations )

% Type of initial condtions 
type_data = 'Vortices';
% 'Vortices' : 2 large anticyclones and 2 large cyclones
%   (used in "Geophysical flow under location uncertainty", Resseguier V.,
%    Memin E., Chapron B.)
% 'Perturbed_vortices' : Same flow with slight small-scale modifications
%   (used in "Chaotic transitions and location uncertainty in geophysical 
%    fluid flows", Resseguier V., Memin E., Chapron B.)
% 'Spectrum' : Gaussian random field with a spectrum slope deined by 
%   the variable slop_b_ini (default value  = -5/3)

% Resolution
resolution = 32;
% The number of grid point is resolution^2
% It has to be an even integer

%% Optional parameters

% deterministic simulation ==> a_H = 0 and only one simulation is performed
k_c = inf; % And then a_H = 0
N_ech=1;
plot_moments = false;

% Spectrum slope of sigma dBt
slope_sigma = - 5/3; 

% Rate between the smallest and the largest wave number of sigma dBt
kappamin_on_kappamax = 1/2;

% Spectrum slope of the initial condition (if type_data = 'Spectrum' )
slope_b_ini = - 5/3; 

% Physical parameters
model = fct_physical_param();

% Gather parameters in the structure model
model.sigma.slope_sigma = slope_sigma;
model.sigma.kappamin_on_kappamax = kappamin_on_kappamax;
model.type_data=type_data;
model.advection.N_ech=N_ech;
model.sigma.k_c = k_c;
model.advection.advection_duration=advection_duration;
model.advection.plot_moments = plot_moments;

%% Generating initial buoyancy
[fft_buoy,model] = fct_buoyancy_init(model,resolution);

%% Advection
[fft_buoy_final, model] = fct_fft_advection_sto(model, fft_buoy);
