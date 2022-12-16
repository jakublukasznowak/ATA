
function [edrFixed,slopeFree,R2,e_edrFixed,e_slopeFree,fig] = edr_sfc (x,dr,fitting_range,component,options)

% EDR_SFC estimates turbulent kinetic dissipation rate with a standard
% method assuming Kolmogorov scaling of the velocity structure function in
% the inertial range of scales as well as estimates the true scaling exponent.
%
% [EDR,S] = edr_sfc(X,DR,FITTING_RANGE,COMPONENT) computes dissipation rate
% EDR and structure function scaling exponent S for a vector X representing
% wind velocity fluctuations, a scalar DR denoting spatial distance between
% consecutive sample points in the record [in meters] (typically DR = 
% (true air speed)/(sampling rate)), a 2-element vector FITTING_RANGE 
% defining the range of scales [in meters] where the fit of the assumed 
% funcional form is performed (it should be a part of the inertial range)
% and a string COMPONENT selecting between longitudinal ('lon') and
% lateral ('lat') velocity component.
% 
% [EDR,S] = edr_sfc(...,'Method',METHOD,'FittingPoints',NPTS) selects 
% the method for evaluating the structure function and the number of points 
% involved in the fit. METHOD can be one of the following:
%
%       'direct'    - evaluates structure function at all
%                     possible displacements which belong to the
%                     FITTING_RANGE and use those values in the fit.
%       'logmean'   - (default) evaluates structure function at all possible
%                     displacements which belong to the FITTING_RANGE but
%                     averages those values in NPTS log-equally-distributed 
%                     bins covering the FITTING_RANGE before the fit is performed
%       'sparse'    - evaluates structure function only at NPTS 
%                     log-equally-distributed displacements within the
%                     FITTING_RANGE.
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
% [...,E_EDR,E_S] = edr_sfc(...) reports the estimated error of the results
% computed by standard error propagation of the fitting errors.
%
% [...,FIG] = edr_sfc(...,'Plot',true) provides the handle to the plot
%
% See also LOGMEAN, EDR_PSD, POLYFIT


arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    fitting_range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty, mustBeValidRange(fitting_range,x,dr)}
    component (1,1) string {mustBeMember(component,{'lon','lat'})}
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean','sparse'})} = 'logmean'
    options.FittingPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 10
    options.Plot (1,1) logical = false
end


r1 = fitting_range(1);
r2 = fitting_range(2);
Np = options.FittingPoints;
Lx = length(x);


% Prepare the list of displacements
if strcmp(options.Method,'sparse')
    iv = unique( round( exp( linspace(log(r1),log(r2),Np) )/dr ) )';
else
    iv = (ceil(r1/dr):r2/dr)';
end
Li = length(iv);


% Calculate structure function values for the displacements from the list
sfc=nan(Li,1);
for i=1:Li
    sfc(i) = mean( ( x(iv(i)+1:Lx) - x(1:Lx-iv(i)) ).^2 );
end


% Average the SFC in log-equal bins
rv = iv*dr;
if strcmp(options.Method,'logmean')
    [rv_fit,sfc_fit] = logmean(rv,sfc,Np);
else
    rv_fit = rv;
    sfc_fit = sfc;
end
Li_fit = length(rv_fit);

if Li_fit<Np
    fprintf('Warning in EDR_SFC: Number of fitting points was reduced to %d.\n',Li_fit)
end

if ~options.Plot
    clear iv rv sfc
end


% Select a constant
if strcmp(component,'lat')
    C = 2.6;
else
    C = 2.0;
end


% Fit (1): free slope
[p,S] = polyfit(log(rv_fit),log(sfc_fit),1);
covarM = (inv(S.R)*inv(S.R)')*S.normr^2/S.df; % covariance matix

slopeFree = p(1);
e_slopeFree = sqrt(covarM(1,1));
offsetFree = p(2);
% e_offsetFree = sqrt(covarM(2,2));

edrFree = (exp(offsetFree)/C)^(1/slopeFree);
% e_edrFree = edr_Free/slope * sqrt( e_offset^2 + (e_slope*log(edrFree))^2 ); % error propagation


% Fit (2): fixed slope
slopeFixed = 2/3;

offsetFixed = mean(log(sfc_fit)-slopeFixed*log(rv_fit));
e_offsetFixed = std(log(sfc_fit)-slopeFixed*log(rv_fit))/sqrt(length(rv_fit)); % standard error from LS fit

edrFixed = (exp(offsetFixed)/C)^(1/slopeFixed);
e_edrFixed = edrFixed/slopeFixed*e_offsetFixed; % error propagation


% Linear correlation
corrM = corrcoef(log(rv_fit),log(sfc_fit));
R2 = corrM(1,2);


% Diagnostic plot
if options.Plot
    
    [fig,ax,co]=fig16x12('loglog',[1 1]);
    set(ax,'XLim',[min(rv) max(rv)],'YLim',[min(sfc) max(sfc)]);
    
    plot(rv,sfc,'.','Color',co(1,:))
    plot(rv_fit,sfc_fit,'^','MarkerFaceColor',co(2,:))
    
    plot(rv_fit,C*(rv_fit*edrFree).^slopeFree,'-','Color',co(4,:),'LineWidth',1)
    plot(rv_fit,C*(rv_fit*edrFixed).^slopeFixed,'-','Color',co(5,:),'LineWidth',1)

    xlabel('$r\,[\mathrm{m}]$','Interpreter','latex')
    ylabel('$D\,[\mathrm{m^2\,s^{-2}}]$','Interpreter','latex')
    
    legend({'SFC','fit points',...
        ['free slope: ',num2str(slopeFree,'%.2f')],...
        ['fixed slope: ','2/3']},...
        'Location','northwest')
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