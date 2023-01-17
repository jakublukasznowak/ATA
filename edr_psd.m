
function [edrFixed,slopeFree,R2,e_edrFixed,e_slopeFree,fig] = edr_psd (x,dr,fitting_range,component,options)

% EDR_PSD estimates turbulent kinetic dissipation rate with a standard
% method assuming Kolmogorov scaling of the velocity spectrum in
% the inertial range of scales as well as estimates the true scaling exponent.
%
% [EDR,S] = edr_psd(X,DR,FITTING_RANGE) computes dissipation rate EDR and
% power spectrum scaling exponent S for a vector X representing wind
% velocity fluctuations, a scalar DR denoting spatial distance between
% consecutive sample points in the record [in meters] (typically DR = 
% (true air speed)/(sampling rate)), a 2-element vector FITTING_RANGE 
% defining the range of scales [in meters] where the fit of the assumed 
% funcional form is performed (it should be a part of the inertial range)
% and a string COMPONENT selecting between longitudinal ('lon') and
% lateral ('lat') velocity component. The power spectrum is evaluated with
% the PWELCH function.
% 
% [EDR,S] = edr_psd(...,'Method',METHOD,'FittingPoints',NPTS) selects 
% the option for evaluating the power spectrum and the number of points 
% involved in the fit. METHOD can be one of the following:
%
%       'direct'    - evaluates power spectrum within the maximum
%                     available range of normalized frequencies W (where W 
%                     corresponds to the spatial scale R and can be obtained
%                     by 2*PI*DR/R), then selects the points which belong to
%                     the FITTING_RANGE and performs a linear fit in
%                     log-log coordinates.
%       'logmean'   - (default) evaluates power spectrum within the maximum available
%                     range of normalized frequencies W, then averages the
%                     values which belong to the FITTING_RANGE in NPTS
%                     log-equally-distributed bins covering the FITTING_RANGE
%                     and performs a linear fit in log-log coordinates
%                     using the averaged points.
%       'sparse'    - evaluates power spectrum at NPTS log-equally-distributed
%                     normalized frequencies W corresponding to the
%                     FITTING_RANGE only and performs a fit on those points.
%
% The default NPTS is 20.
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
% [...,R2] = edr_psd(...) reports the linear correlation coefficient of
% power spectrum values versus normalized freuencies W in log-log
% coordinates for the fitting range.
%
% [...,E_EDR,E_S] = edr_psd(...) reports the estimated error of the results
% computed by standard error propagation of the fitting errors.
%
% [...,FIG] = edr_sfc(...,'Plot',true) provides the handle to the plot
%
% See also PWELCH, LOGMEAN, POLYFIT, EDR_SFC


arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    fitting_range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty, mustBeValidRange(fitting_range,x,dr)}
    component (1,1) string {mustBeMember(component,{'lon','lat'})}
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean','sparse'})} = 'logmean'
    options.FittingPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 20
    options.Plot (1,1) logical = false
    options.WindowLength (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = floor(length(x)/2)
    options.WindowOverlap (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = ceil(length(x)/4)
end


r1 = fitting_range(1); 
r2 = fitting_range(2);
if r2>dr*options.WindowLength/2
    r2 = dr*options.WindowLength/2;
    fprintf('Warning in EDR_PSD: Fitting range modified to [%.2f %.2f] to comply with the pwelch window length.\n',r1,r2)
end
w1 = 2*pi*dr/r2;
w2 = 2*pi*dr/r1; 

Np = options.FittingPoints; 
Nfft = 2^nextpow2(options.WindowLength);


% Calculate PSD
if strcmp(options.Method,'sparse')
    wv_fit = exp( linspace(log(w1),log(w2),Np) )';
    psd_fit = 2*pwelch(x,options.WindowLength,options.WindowOverlap,wv_fit); % 2*onesided spectrum
else
    [psd,wv] = pwelch(x,options.WindowLength,options.WindowOverlap,Nfft);
    psd = psd(wv>0); wv = wv(wv>0);
    ind1 = find(wv>=w1,1,'first');
    ind2 = find(wv<=w2,1,'last');
    wv_fit = wv(ind1:ind2);
    psd_fit = psd(ind1:ind2);
    if strcmp(options.Method,'logmean')
        [wv_fit,psd_fit] = logmean(wv_fit,psd_fit,Np);
    end
end

if ~options.Plot
    clear wv psd
end


% Select a constant
if strcmp(component,'lat')
    C = 0.65;
else
    C = 0.49;
end


% Fit (1): free slope
[p,S] = polyfit(log(wv_fit),log(psd_fit),1);
covarM = (inv(S.R)*inv(S.R)')*S.normr^2/S.df; % covariance matix

slopeFree = p(1);
e_slopeFree = sqrt(covarM(1,1));
offsetFree = p(2);
% e_offsetFree = sqrt(covarM(2,2));

edrFree = (exp(offsetFree)/C)^(3/2)/dr;
% e_edrFree = 3/2*edrFree*e_offsetFree; % error propagation


% Fit (2): fixed slope
slopeFixed = -5/3;

offsetFixed = mean(log(psd_fit)-slopeFixed*log(wv_fit));
e_offsetFixed = std(log(psd_fit)-slopeFixed*log(wv_fit))/sqrt(length(wv_fit)); % standard error from LS fit

edrFixed = (exp(offsetFixed)/C)^(3/2)/dr;
e_edrFixed = 3/2*edrFixed*e_offsetFixed; % error propagation


% Linear correlation
corrM = corrcoef(log(wv_fit),log(psd_fit));
R2 = corrM(1,2);


% Diagnostic plot
if options.Plot
    
    rv_fit = 2*pi*dr./wv_fit;
    if strcmp(options.Method,'sparse')
        rv = rv_fit; psd = psd_fit;
    else
        rv = 2*pi*dr./wv;
    end
    
    [fig,ax,co]=fig16x12('loglog',[1 1]);
%     set(ax,'XLim',[min(rv) max(rv)],'YLim',[min(psd) max(psd)],'XDir','reverse');
    set(ax,'XLim',[min(rv) fitting_range(2)],'YLim',[min(psd) max(psd)]);
    
    plot(rv,psd,'.','Color',co(1,:))
    plot(rv_fit,psd_fit,'^','MarkerFaceColor',co(2,:))
    
    plot(rv_fit,C*(dr*edrFree).^(2/3).*wv_fit.^slopeFree,'-','Color',co(4,:),'LineWidth',1)
    plot(rv_fit,C*(dr*edrFixed).^(2/3).*wv_fit.^slopeFixed,'-','Color',co(5,:),'LineWidth',1)

    xlabel('$r\,[\mathrm{m}]$','Interpreter','latex')
    ylabel('$P\,[\mathrm{m^2\,s^{-2}\,rad^{-1}}]$','Interpreter','latex')
    
    legend({'PSD','fit points',...
        ['free slope: ',num2str(slopeFree,'%.2f')],...
        ['fixed slope: ','-5/3']},...
        'Location','northwest')
%     text(0.05,0.10,['$\epsilon = ',sprintf('%.2f',edrFixed/10^floor(log10(edrFixed))),'\cdot10^',...
%         sprintf('{%d}',floor(log10(edrFixed))),'\,\mathrm{m^2\,s^{-3}}$',...
%         newline,'$R = ',sprintf('%.3f',R2),'$'],...
%         'FontSize',12,'Units','Normalized',...
%         'HorizontalAlignment','left','Interpreter','latex')
    text(0.66,0.10,['$\epsilon = ',sprintf('%.2f',edrFixed/10^floor(log10(edrFixed))),'\cdot10^',...
        sprintf('{%d}',floor(log10(edrFixed))),'\,\mathrm{m^2\,s^{-3}}$',...
        newline,'$R = ',sprintf('%.3f',R2),'$'],...
        'FontSize',12,'Units','Normalized',...
        'HorizontalAlignment','left','Interpreter','latex')
    
else
    fig = [];
end


end


function mustBeValidRange(a,x,dr)
    if ~ge(a(1),dr*2)
        eid = 'Range:firstTooLow';
        msg = sprintf('Fitting range must be within [dr*2 dr*length(x)/2] = [%.2f %.2f].',dr*2,dr*length(x)/2);
        throwAsCaller(MException(eid,msg))
    end
    if ~le(a(2),length(x)*dr/2)
        eid = 'Range:lastTooHigh';
        msg = sprintf('Fitting range must be within [dr*2 dr*length(x)/2] = [%.2f %.2f].',dr*2,dr*length(x)/2);
        throwAsCaller(MException(eid,msg))
    end
    if ge(a(1),a(2))
        eid = 'Range:notIncreasing';
        msg = 'Fitting range must be of nonzero length.';
        throwAsCaller(MException(eid,msg))
    end
end