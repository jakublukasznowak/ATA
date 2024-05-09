
function [edrFixed,slpFree,e,fig] = edr_sfc (x,dr,fit_range,C,options)

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
% [EDR,S] = edr_sfc(...,'Plot',PLT) selects whether to show a diagnostics
% plot. PLT can be true or false (default).
%
% [...,R2] = edr_sfc(...) reports the linear correlation coefficient of
% structure function values versus displacement in log-log coordinates for
% the fitting range.
%
% [...,E] = edr_sfc(...) reports the estimated error structure with the following
% fields: offsetFixed, edrFixed, offsetFree, slopeFree, R2
%
% [...,FIG] = edr_sfc(...,'Plot',true) provides the handle to the plot
%
% See also LOGMEAN, EDR_PSD, POLYFIT



arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    fit_range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty, mustBeValidRange(fit_range,x,dr)}
    C (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = 2.0
    options.Slope (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = 2/3
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean'})} = 'logmean'
    options.FitPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 10
    options.Plot (1,1) logical = false
    options.PlotXLim (1,2) {mustBePositive, mustBeFinite, mustBeNonempty, mustBeValidRange(options.PlotXLim,x,dr)} = fit_range
    options.PlotYLim (1,2) {mustBeReal, mustBeNonempty} = [-inf inf];
end


% Prepare the list of displacements

iv = ( ceil(fit_range(1)/dr) : fit_range(2)/dr )';
rv = iv*dr;
Li = length(iv);
Lx = length(x);


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
Li_fit = length(rv_fit);

if Li_fit<options.FitPoints
    if ~strcmp(options.Method,'direct')
        fprintf('Warning in EDR_SFC: Number of fitting points was reduced to %d.\n',Li_fit)
    end
end


% Fit (1): fixed slope

slpFixed = options.Slope;

offsetFixed = mean(log(sfc_fit)-slpFixed*log(rv_fit));
edrFixed = (exp(offsetFixed)/C)^(1/slpFixed);

e.offsetFixed = std(log(sfc_fit)-slpFixed*log(rv_fit))/sqrt(length(rv_fit)); % standard error from LS fit
e.edrFixed = edrFixed/slpFixed*e.offsetFixed; % error propagation


% Fit (2): free slope

[p,S] = polyfit(log(rv_fit),log(sfc_fit),1);

slpFree = p(1);
offsetFree = p(2);
edrFree = (exp(offsetFree)/C)^(1/slpFree);

covarM = (inv(S.R)*inv(S.R)')*S.normr^2/S.df; % covariance matix
e.slopeFree  = sqrt(covarM(1,1));
e.offsetFree = sqrt(covarM(2,2));
e.edrFree    = edrFree/slpFree * sqrt( e.offsetFree^2 + (e.slopeFree*log(edrFree))^2 ); % error propagation


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
    
    plot(rv_fit,C*(rv_fit*edrFree).^slpFree,'-','Color',co(4,:),'LineWidth',2)
    plot(rv_fit,C*(rv_fit*edrFixed).^slpFixed,'-','Color',co(5,:),'LineWidth',2)

    xlabel('$r\,[\mathrm{m}]$','Interpreter','latex')
    ylabel('$D\,[\mathrm{m^2\,s^{-2}}]$','Interpreter','latex')
    
    legend({'$D$','fit points',...
        ['$s=$ ',num2str(slpFree,'%.2f')],...
        '$s=$ 2/3'},...
        'Location','northwest','Interpreter','latex')
%     text(0.66,0.10,['$\epsilon = ',sprintf('%.2f',edrFixed/10^floor(log10(edrFixed))),'\cdot10^',...
%         sprintf('{%d}',floor(log10(edrFixed))),'\,\mathrm{m^2\,s^{-3}}$',...
%         newline,'$R = ',sprintf('%.3f',e.R2),'$'],...
%         'FontSize',12,'Units','Normalized',...
%         'HorizontalAlignment','left','Interpreter','latex')
    
else
    fig = [];
end


end


function mustBeValidRange(a,x,dr)
    if ~ge(a(1),dr)
        eid = 'Range:firstTooLow';
        msg = sprintf('Fitting range must be within [dr dr*length(x)] = [%.2f %.2f].',dr,dr*length(x));
        throwAsCaller(MException(eid,msg))
    end
    if ~le(a(2),length(x)*dr)
        eid = 'Range:lastTooHigh';
        msg = sprintf('Fitting range must be within [dr dr*length(x)] = [%.2f %.2f].',dr,dr*length(x));
        throwAsCaller(MException(eid,msg))
    end
    if ge(a(1),a(2))
        eid = 'Range:notIncreasing';
        msg = 'Fitting range must be of nonzero length.';
        throwAsCaller(MException(eid,msg))
    end
end