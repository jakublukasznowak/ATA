
function [edr,slp,e,fig] = edr_psd (x,dr,fit_range,C,options)

% EDR_PSD estimates turbulent kinetic dissipation rate with a standard
% method assuming Kolmogorov scaling of the velocity spectrum in
% the inertial range of scales as well as estimates the true scaling exponent.
%
% [EDR,S] = edr_psd(X,DR,FIT_RANGE) computes dissipation rate EDR and
% power spectrum scaling exponent S for a vector X representing wind
% velocity fluctuations, a scalar DR denoting spatial distance between
% consecutive sample points in the record [in meters] (typically DR = 
% (true air speed)/(sampling rate)), a 2-element vector FITTING_RANGE 
% defining the range of scales [in meters] where the fit of the assumed 
% funcional form is performed (it should be a part of the inertial range)
% and a proportionality constant C which is typically 0.5 for longitudinal
% and 0.65 for transverse velocity components.
% The power spectrum is evaluated with PWELCH.
% 
% [EDR,S] = edr_psd(...,'Method',METHOD,'FitPoints',NPTS) selects 
% the option for evaluating the power spectrum and the number of points 
% involved in the fit. METHOD can be one of the following:
%
%       'direct'    - evaluates power spectrum within the maximum
%                     available range of normalized frequencies W (where W 
%                     corresponds to the spatial scale R and can be obtained
%                     by 2*PI*DR/R), then selects the points which belong to
%                     the FIT_RANGE and performs a linear fit in
%                     log-log coordinates.
%       'logmean'   - (default) evaluates power spectrum within the maximum available
%                     range of normalized frequencies W, then averages the
%                     values which belong to the FIT_RANGE in NPTS
%                     log-equally-distributed bins covering the FIT_RANGE
%                     and performs a linear fit in log-log coordinates
%                     using the averaged points.
%
% The default NPTS is 20.
%
% [EDR,S] = edr_sfc(...,'SLOPE',S) allows to use another fixed slope than
% 5/3
%
% [EDR,S] = edr_psd(...,'Plot',PLT) selects whether to show a diagnostics
% plot. PLT can be true or false (default).
%
% [EDR,S] = edr_psd(...,'WindowLength',WINDOW,'WindowOverlap',OVERLAP)
% determine the size of the segments (in number of points) which are used
% in the computation of the modified periodograms by the PWELCH function and
% their overlap. Those arguments are simply passed as the inputs to PWELCH.
% The default WINDOW is half the length of the signal X. The default
% OVERLAP is a quarter of the length of the signal X.
%
% [...,E] = edr_sfc(...) reports the estimated error structure with the following
% fields: edr, slp, R2, N, O, Ostar
%
% [...,FIG] = edr_sfc(...,'Plot',true) provides the handle to the plot
%
% See also PWELCH, LOGMEAN, POLYFIT, EDR_SFC


arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    fit_range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty}
    C (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = 0.5
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean'})} = 'logmean'
    options.FitPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 20
    options.WindowLength (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = floor(length(x)/2)
    options.WindowOverlap (1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite, mustBeNonempty} = ceil(length(x)/4)
    options.Slope (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = -5/3
    options.Plot (1,1) logical = false
    options.PlotXLim (1,2) {mustBePositive, mustBeFinite, mustBeNonempty} = fit_range
    options.PlotYLim (1,2) {mustBeReal, mustBeNonempty} = [-inf inf]
end


% Check if the fit range is valid

if fit_range(1)<2*dr || fit_range(2)>dr*options.WindowLength
    if fit_range(1)<2*dr
        fit_range(1) = 2*dr;
    else
        fit_range(2) = dr*options.WindowLength;
    end
    warning('EDR_PSD:InvalidFitRange','Invalid fit range was changed to [%.2f %.2f].',fit_range(1),fit_range(2))
end


% Fit range bounds

r1 = fit_range(1); 
r2 = fit_range(2);
w1 = 2*pi*dr/r2;
w2 = 2*pi*dr/r1; 


% Calculate power spectrum

Nfft = 2^nextpow2(options.WindowLength);

[psd,wv] = pwelch(x,options.WindowLength,options.WindowOverlap,Nfft);
psd = psd(wv>0); wv = wv(wv>0);


% Select fitting range

ind1 = find(wv>=w1,1,'first');
ind2 = find(wv<=w2,1,'last');
wv_fit = wv(ind1:ind2);
psd_fit = psd(ind1:ind2);
Li = length(wv_fit);

if Li<2
    throw(MException('EDR_PSD:TooFewFitPoints','Number of fit points must be at least 2.'))
end


% Average PSD in log-equal bins (for LOGMEAN method)

if strcmp(options.Method,'logmean')
    [wv_fit,psd_fit] = logmean(wv_fit,psd_fit,options.FitPoints);
end
e.N = length(psd_fit);


% Fit (1): fixed slope

slpFixed = options.Slope;

logO = mean(log(psd_fit)-slpFixed*log(wv_fit));
e.logO = std(log(psd_fit)-slpFixed*log(wv_fit))/sqrt(length(wv_fit)); % standard error

O = exp(logO);
e.O = O*e.logO;

edr = (O/C)^(3/2)/dr;
e.edr = 3/2*edr*e.logO;


% Fit (2): free slope

[p,S] = polyfit(log(wv_fit),log(psd_fit),1);

slp = p(1);
logOstar = p(2);

covarM = (inv(S.R)*inv(S.R)')*S.normr^2/S.df; % covariance matix
e.slp  = sqrt(covarM(1,1));
e.logOstar = sqrt(covarM(2,2));

Ostar = exp(logOstar);
e.Ostar = Ostar*e.logOstar;


% Linear correlation

corrM = corrcoef(log(wv_fit),log(psd_fit));
e.R2 = corrM(1,2);


% Diagnostic plot

if options.Plot
    
    rv = 2*pi*dr./wv;
    rv_fit = 2*pi*dr./wv_fit;
  
    [fig,~,co] = fig16x12('loglog',[1 1],'on','XLim',options.PlotXLim,'YLim',options.PlotYLim);
    
    plot(rv,psd,'.','Color',co(1,:),'MarkerSize',8)
    plot(rv_fit,psd_fit,'^','MarkerFaceColor',co(2,:),'MarkerSize',8)
    
    plot(rv_fit,Ostar*wv_fit.^slp,'-','Color',co(4,:),'LineWidth',2)
    plot(rv_fit,C*(dr*edr).^(2/3).*wv_fit.^slpFixed,'-','Color',co(5,:),'LineWidth',2)

    xlabel('$r\,[\mathrm{m}]$','Interpreter','latex')
    ylabel('$P\,[\mathrm{m^2\,s^{-2}\,rad^{-1}}]$','Interpreter','latex')
    
    legend({'$P$','fit points',...
        ['$p=$ ',num2str(abs(slp),'%.2f')],...
        '$p=$ 5/3'},...
        'Location','northwest','Interpreter','latex')
    text(0.66,0.10,['$\epsilon = ',sprintf('%.2f',edr/10^floor(log10(edr))),'\cdot10^',...
        sprintf('{%d}',floor(log10(edr))),'\,\mathrm{m^2\,s^{-3}}$',...
        newline,'$R = ',sprintf('%.3f',abs(e.R2)),'$'],...
        'FontSize',12,'Units','Normalized',...
        'HorizontalAlignment','left','Interpreter','latex')
    
else
    fig = [];
end


end
