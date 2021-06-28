function [fft_sst,model] = fct_buoyancy_init(model,resolution)
% Initial condition for the buoyancy 

%% Grid
n=resolution;m=n;
Lx=1e6;
dx=Lx/n;dy=dx;
x= dx*(0:n-1);y= dy*(0:m-1);
model.grid.origin=[0 0];
model.grid.x_ref=x;model.grid.y_ref=y;
[X,Y]=ndgrid(x,y);
model.grid.dX=[dx dy];
MX=[n m];model.grid.MX=MX;

ee=4;
odg_b = model.odg_b;
x=X(:,1);y=Y(1,:);
sigma= 2 * Lx/15; % Length scale close to the Rossby radius

%% Spatial buoyancy field 
%Warm anticyclones    
center1x=x(1/4*model.grid.MX(1)+1);center1y=y(1/4*model.grid.MX(2)+1);
b_S = + exp( -1/(2*sigma^2)* (ee*(X-center1x).^2+(Y-center1y).^2) );
center2x=x(3/4*model.grid.MX(1)+1);center2y=y(1/4*model.grid.MX(2)+1);
b_S = b_S + exp( -1/(2*sigma^2)* (ee*(X-center2x).^2+(Y-center2y).^2) );
% Cold cyclones     
center1x=x(1/4*model.grid.MX(1)+1);center1y=y(3/4*model.grid.MX(2)+1);
b_S = b_S - exp( -1/(2*sigma^2)* (ee*(X-center1x).^2+(Y-center1y).^2) );
center2x=x(3/4*model.grid.MX(1)+1);center2y=y(3/4*model.grid.MX(2)+1);
b_S = b_S - exp( -1/(2*sigma^2)* (ee*(X-center2x).^2+(Y-center2y).^2) );
% Specify the amplitude of the buoyancy
b_S = odg_b* b_S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fft_sst = fft2(b_S);
        


