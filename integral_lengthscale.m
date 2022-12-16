
function [L,fig] = integral_lengthscale (x,dr,options)

% INTEGRAL_LENGTHSCALE estimates the integral length scale of turbulence
% based on the autocorrelation function of turbulent fluctuations.
%
% L = integral_lengthscale(X,DR) computes the integral lengthscale L for a
% vector X representing turbulent fluctuations and a scalar DR denoting
% spatial distance between consecutive sample points in the record [in
% meters], typically DR = (true air speed)/(sampling rate).
%
% L = integral_lengthscale(...,'Method',METHOD) selects the method for 
% estimating L. METHOD can be one of the following:
%
%       'e-decay'     - (default) evaluates L as the distance at which
%                       autocorrelation function decays e-times
%       'integration' - evaluates L as the integral of the autocorrelation 
%                       function up to its first zero. If there is no
%                       zero-crossing in the available range of scales, the
%                       integral is calculated up to the largest scale
%                       available
%
% L = integral_lengthscale(...,'Sampling',SMPL) controls the
% sampling of the autocorrelation function by defining the distance between
% the evaluation points [in meters]. SMPL is rounded to the closest
% multiple of DR.
%
% L = integral_lengthscale(...,'CrossCorrelatedSignal',Y) allows to
% estimate the length scale based on the crosscorrelation function between
% two sinals X and Y instead of the autocorrelation function of X.
%
% L = integral_lengthscale(...,'Plot',PLT) selects whether to show
% a diagnostics plot. PLT can be true or false (default).
%
% [...,FIG] = edr_sfc(...,'Plot',true) provides the handle to the plot
%
% See also REYNOLDS_DECOMPOSITION, TRAPZ


arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    options.Method (1,1) string {mustBeMember(options.Method,{'e-decay','integration'})} = 'e-decay'
    options.Sampling (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = dr
    options.CrossCorrelatedSignal (:,1) {mustBeReal, mustBeFinite}  = []
    options.Plot (1,1) logical = false
end


% Prepare the list of displacements

Lx = length(x);
iv = (0:floor(options.Sampling/dr):Lx/2-1)';
Li = length(iv);
rv = iv*dr;


% Prepare second vector

if isempty(options.CrossCorrelatedSignal)
    y=x;
else
    y=options.CrossCorrelatedSignal;
end


% Calculate correlation function

Fxy = mean(x.*y);
R = nan(Li,1);
for i=1:Li
    R(i) = mean( x(1:Lx-iv(i)) .* y(iv(i)+1:Lx) );
end
rho = R/Fxy;


% Estimate integral scale

if strcmp(options.Method,'e-decay')
    
    ind0 = find(rho(1:end-1).*rho(2:end)<0,1,'first');
    if isempty(ind0)
        ind0 = Li;
    end
    L = interp1( rho(1:ind0), rv(1:ind0), exp(-1) );
    
elseif strcmp(options.Method,'integration')
    
    ind0 = find(rho(1:end-1).*rho(2:end)<0,1,'first');
    if isempty(ind0)
        fprintf('Warning: Correlation function has no zero. Integrating up to r = %.1f m.\n',rv(end))
        L = trapz( rv, rho );
    else
        r0 = interp1( rho(ind0:ind0+1), rv(ind0:ind0+1), 0);
        L = trapz( [rv(1:ind0);r0], [rho(1:ind0);0] );
    end
    
else
    error('Invalid method selected.)
    
end


% Plot

if options.Plot
    
    [fig,ax,co] = fig16x12;
    set(ax,'XLim',[0 rv(end)]);
    plot(rv,rho,'.','Color',co(1,:))
    plot([0,rv(end)],[exp(-1) exp(-1)],'LineWidth',1,'Color',co(3,:))
    plot(L, interp1(rv,rho,L,'linear','extrap'),'.','MarkerSize',30,'Color',co(2,:))
    legend({'$\rho(r)$','$e^{-1}$',['$L=',sprintf('%5.1f',L),'\,\mathrm{m}$']},'Interpreter','latex')
    xlabel('$r\,\mathrm{[m]}$','Interpreter','latex')
    ylabel('$\rho(r)$','Interpreter','latex')
    
else
    fig = [];
end
    

end