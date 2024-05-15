
function [edr,slp,e,fig] = edr_sfc (x,dr,fit_range,C,options)

% EDR_SFC estimates turbulent kinetic dissipation rate with a standard
% method assuming Kolmogorov scaling of the velocity structure function in
% the inertial range of scales as well as estimates the true scaling exponent.
%
% [EDR,S] = edr_sfc(X,DR,FIT_RANGE,C) computes dissipation rate
% EDR and structure function scaling exponent S for a vector X representing
% wind velocity fluctuations, a scalar DR denoting spatial distance between
% consecutive sample points in the record [in meters] (typically DR = 
% (true air speed)/(sampling rate)), a 2-element vector FIT_RANGE 
% defining the range of scales [in meters] where the fit of the assumed 
% funcional form is performed (it should be a part of the inertial range)
% and a proportionality constant C which is typically 2.0 for longitudinal 
% and 2.6 for transverse velocity components.
% 
% [EDR,S] = edr_sfc(...,'Method',METHOD,'FitPoints',NPTS) selects 
% the method for evaluating the structure function and the number of points 
% involved in the fit. METHOD can be one of the following:
%
%       'direct'    - evaluates structure function at all
%                     possible displacements which belong to the
%                     FIT_RANGE and use those values in the fit.
%       'logmean'   - (default) evaluates structure function at all possible
%                     displacements which belong to the FIT_RANGE but
%                     averages those values in NPTS log-equally-distributed 
%                     bins covering the FIT_RANGE before the fit is performed
%
% The default NPTS is 10.
%
% [EDR,S] = edr_sfc(...,'SLOPE',S) allows to use another fixed slope than 2/3
%
% [EDR,S] = edr_sfc(...,'Plot',PLT) selects whether to show a diagnostics
% plot. PLT can be true or false (default).
%
% [...,E] = edr_sfc(...) reports the estimated error structure with the following
% fields: edr, slp, R2, N, O, Ostar
%
% [...,FIG] = edr_sfc(...,'Plot',true) provides the handle to the plot
%
% See also LOGMEAN, EDR_PSD, POLYFIT



arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    fit_range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty}
    C (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = 2.0
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean'})} = 'logmean'
    options.FitPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 10
    options.Slope (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = 2/3
    options.Detrend (1,1) logical = false
    options.Plot (1,1) logical = false
    options.PlotXLim (1,2) {mustBePositive, mustBeFinite, mustBeNonempty} = fit_range
    options.PlotYLim (1,2) {mustBeReal, mustBeNonempty} = [-inf inf];
end


% Detrend

if options.Detrend
    x = detrend(x);
end


% Check if the fit range is valid

Lx = length(x);

if fit_range(1)<dr || fit_range(2)>Lx*dr
    if fit_range(1)<dr
        fit_range(1) = dr;
    else
        fit_range(2) = Lx*dr;
    end
    warning('EDR_SFC:InvalidFitRange','Invalid fit range was changed to [%.2f %.2f].',fit_range(1),fit_range(2))
end


% Prepare the list of displacements

iv = ( ceil(fit_range(1)/dr) : fit_range(2)/dr )';
rv = iv*dr;
Li = length(iv);

if Li<2
    throw(MException('EDR_SFC:TooFewFitPoints','Number of fit points must be at least 2.'))
end


% Calculate structure function values for the displacements from the list

sfc = nan(Li,1);
for i = 1:Li
    sfc(i) = mean( ( x(iv(i)+1:Lx) - x(1:Lx-iv(i)) ).^2 );
end


% Average SFC in log-equal bins (for LOGMEAN method)

if strcmp(options.Method,'logmean')
    [rv_fit,sfc_fit] = logmean(rv,sfc,options.FitPoints);
else
    rv_fit = rv;
    sfc_fit = sfc;
end
e.N = length(sfc_fit);


% Fit (1): fixed slope

slpFixed = options.Slope;

logO = mean(log(sfc_fit)-slpFixed*log(rv_fit));
e.logO = std(log(sfc_fit)-slpFixed*log(rv_fit))/sqrt(length(rv_fit)); % standard error

O = exp(logO);
e.O = O*e.logO;

edr = (O/C)^(1/slpFixed);
e.edr = edr/slpFixed*e.logO;


% Fit (2): free slope

[p,S] = polyfit(log(rv_fit),log(sfc_fit),1);

slp = p(1);
logOstar = p(2);

covarM = (inv(S.R)*inv(S.R)')*S.normr^2/S.df; % covariance matix
e.slp  = sqrt(covarM(1,1));
e.logOstar = sqrt(covarM(2,2));

Ostar = exp(logOstar);
e.Ostar = Ostar*e.logOstar;


% Linear correlation

corrM = corrcoef(log(rv_fit),log(sfc_fit));
e.R2 = corrM(1,2);


% Diagnostic plot

if options.Plot
    
    iv = ( ceil(options.PlotXLim(1)/dr) : options.PlotXLim(2)/dr )';
    rv = iv*dr;
    Li = length(iv);
    Lx = length(x);

    sfc = nan(Li,1);
    for i = 1:Li
        sfc(i) = mean( ( x(iv(i)+1:Lx) - x(1:Lx-iv(i)) ).^2 );
    end

    
    [fig,~,co] = fig16x12('loglog',[1 1],'on','XLim',options.PlotXLim,'YLim',options.PlotYLim);
    
    plot(rv,sfc,'.','Color',co(1,:),'MarkerSize',8)
    plot(rv_fit,sfc_fit,'^','MarkerFaceColor',co(2,:),'MarkerSize',8)
    
    plot(rv_fit,Ostar*rv_fit.^slp,'-','Color',co(4,:),'LineWidth',2)
    plot(rv_fit,C*(rv_fit*edr).^slpFixed,'-','Color',co(5,:),'LineWidth',2)

    xlabel('$r\,[\mathrm{m}]$','Interpreter','latex')
    ylabel('$D\,[\mathrm{m^2\,s^{-2}}]$','Interpreter','latex')
    
    legend({'$D$','fit points',...
        ['$s=$ ',num2str(slp,'%.2f')],...
        '$s=$ 2/3'},...
        'Location','northwest','Interpreter','latex')
    text(0.66,0.10,['$\epsilon = ',sprintf('%.2f',edr/10^floor(log10(edr))),'\cdot10^',...
        sprintf('{%d}',floor(log10(edr))),'\,\mathrm{m^2\,s^{-3}}$',...
        newline,'$R = ',sprintf('%.3f',e.R2),'$'],...
        'FontSize',12,'Units','Normalized',...
        'HorizontalAlignment','left','Interpreter','latex')
    
else
    fig = [];
end


end
