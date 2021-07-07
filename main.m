%%%%%%%%%%%%%%%%%%%%
%%% Main 
%%%%%%%%%%%%%%%%%%%%
% Cleaning
clear all;close all;clc;
dbstop if error;% mode debug
% Paths
addpath(pwd)
% Cleaning
home;

%% Main parameters to choose
% Duration of the simulation (in seconds)
advection_duration = 3600*24*20; % 20 days

% Resolution (even integer)
resolution = 64;

% Gather parameters in the structure model
angle_grid = 45/360*2*pi; % rad
OMEGA = 2*pi/(24*60*60);
model.physical_constant.f0 = 2* OMEGA * sin( angle_grid );
model.physical_constant.rho=1e3;% Background density
model.physical_constant.g=9.81;% Gravity
model.physical_constant.buoyancy_freq_N = 3 * model.physical_constant.f0;% Background stratification
model.odg_b = 1e-3;% Amplitude of the buoyancy

%%%
model.advection.advection_duration=advection_duration;

%% Initial condition for the buoyancy
[fft_buoy,model] = fct_buoyancy_init(model,resolution);

%% Advection
[fft_buoy_final, model] = fct_fft_advection_sto(model, fft_buoy);
