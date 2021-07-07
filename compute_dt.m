function [ dt ] = compute_dt(model, w)
%compute dt 

% CFL of the (large-scale) advection
dX=model.grid.dX;
bound1=sum(abs(w)*pi/dX(1),3);
bound1=max(bound1(:));
bound1=1/bound1/4;
% CFL of the hyperviscosity
dX2=(model.grid.dX /pi).^2;
bound2=1/model.advection.HV.val*(prod(dX2)/sum(dX2))^(model.advection.HV.order/2);

% Minimum of the CFL
dt = min([bound1 bound2]); 
clear  bound1 bound2 dX dX2

end

