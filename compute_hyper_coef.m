function [ coef ] = compute_hyper_coef( model, w )
%compute hyperviscosity coefficient 

% compute gradient of w using 'gradient' matlab function
[dxUx,dyUx] = gradient(w(:,:,1)',model.grid.x_ref,model.grid.y_ref);
[dxUy,dyUy] = gradient(w(:,:,2)',model.grid.x_ref,model.grid.y_ref);
s=(dxUx+dyUy)/2;
d=(dyUx+dxUy)/2;
lambda = sqrt(s.^2+d.^2);
model.advection.lambda_RMS = sqrt(mean(lambda(:).^2));

% Hyperviscosity coefficient
coef = 40 * sqrt(mean(lambda(:).^2)) * ...
    (mean(model.grid.dX)/pi)^model.advection.HV.order;

end

