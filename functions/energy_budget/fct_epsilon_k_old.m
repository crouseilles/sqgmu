function int_epsilon = fct_epsilon_k_old(model,ft,int_epsilon_dt_m_1,color)
% Compute the spectrum of a function and superimposed a slope -5/3
%

% Color by default
if nargin < 3
    color='b';
end

% Square modulus of the Fourier transform
ft=abs(ft).^2;
ft= sum(ft,3);

switch model.dynamics
    case 'SQG'
        slope_ref = -5/3;
    case '2D'
        slope_ref = -3;
    otherwise
        error('Unknown type of dynamics');
end

% Get parameters
MX=model.grid.MX;
PX=MX/2;
dX=model.grid.dX;
if any(size(ft)~=MX)
    error('wrong size');
end
if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end

persistent idxref MXref dXref kidx

if ~exist('MXref','var') ||  isempty (MXref) || any(MXref ~= MX) ...
        || any(dXref ~= dX)
    MXref=MX;
    dXref=dX;
    
    % Remove aliasing
    ft(PX(1),:)=0;
    ft(:,PX(2))=0;
    
    %% Wave vector
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
    kx=2*pi/model.grid.dX(1)*kx;
    ky=2*pi/model.grid.dX(2)*ky;
    [kx,ky]=ndgrid(kx,ky);
    k=sqrt(kx.^2+ky.^2);
    k(PX(1)+1,:)=inf;
    k(:,PX(2)+1)=inf;
    k=k(:);
    
    %% Wave number
    M_kappa=min(MX);
    P_kappa= M_kappa/2;
    d_kappa = max(1./dX);
    kidx=1/(M_kappa)* (0:(P_kappa-1)) ;
    kidx=2*pi*d_kappa*kidx;
    
    %% Masks associated with the rings of iso wave number
    d_kappa = kidx(2) - kidx(1);
    idx = sparse( bsxfun(@le,kidx, k ) );
    idx = idx & sparse( bsxfun(@lt,k, [ kidx(2:end) kidx(end)+d_kappa ] ) );
    idxref=idx;
end

%% Spectrum
% Integration over the rings of iso wave number
spectrum = idxref' * ft(:);

if strcmp(model.dynamics,'2D')
    spectrum = spectrum ./(kidx').^2;
end

% Division by prod(model.grid.MX) because of the Parseval theorem for
% discrete Fourier transform
% Division by prod(model.grid.MX) again in order to the integration 
% of the spectrum over the wave number yields the energy of the
% buoyancy averaged (not just integrated) over the space
spectrum = 1/prod(model.grid.MX)^2 * spectrum;

% Division by the wave number step
d_kappa = kidx(2)-kidx(1);
spectrum = spectrum / d_kappa;

% Time integral of the dissipation per scale epsilon(k)
int_epsilon = - cumsum(spectrum) * d_kappa;

epsilon = (int_epsilon - int_epsilon_dt_m_1)/model.advection.dt_adv;

%% Plot
slope_ref = 0;
idx_not_inf=~(isinf(log10(epsilon(2:end))) ...
    | epsilon(2:end)<1e-4*max(epsilon(2:end)) | isinf(kidx(2:end)'));
idx_not_inf = [false; idx_not_inf];
% line1= slope_ref * log10(kidx(2:end))  ;
line1= slope_ref * log10(kidx(2:end))  ;
offset = -1 + mean(  log10(epsilon(idx_not_inf)')  ...
    - line1(idx_not_inf(2:end)));
line1 = line1 + offset;
ref=10.^line1;
semilogx(kidx(2:end),ref,'--k');
hold on;
name_plot =semilogx(kidx(2:end) , epsilon(2:end) ,color);
ax=axis;
ax(4)=max([epsilon(2:end); ref']);
% min_ax= 10 ^(slope_ref * log10(kidx(2)*512/2) + offset) ;
% ax(3) = (model.odg_b/(1e-3))^2 * ...
%     6e-2*(kidx(2)/kidx(end))*min([max(epsilon); max(ref)']);
% % ax(3)=6e-2*(kidx(2)/kidx(end))*min([max(epsilon); max(ref)']);
% ax(3) = min( [ax(3) min(ref) min_ax]);
ax(1:2)=kidx(2)*[1 min(model.grid.MX)/2];
if ax(4)>0
    axis(ax)
end
