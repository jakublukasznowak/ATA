
function [edr,cf,fig1,fig2] = edr_var (x,dr,cutoff_scale,component,options)

% EDR_VAR estimates turbulent kinetic dissipation rate with an iterative
% method based on correcting the measured variance of velocity derivatives
% by assuming specific functional form of the energy spectrum in the
% dissipative range of scales. The method is described in sec. 2.c. of the 
% paper:
%
%   Akinlabi E.O., Wacławczyk M., Mellado J.P., and Malinowski S.P., 2019:
%   Estimating Turbulence Kinetic Energy Dissipation Rates in the Numerically
%   Simulated Stratocumulus Cloud-Top Mixing Layer: Evaluation of Different 
%   Methods, Journal of the Atmospheric Sciences, vol. 76(5), pp. 1471–1488,
%   10.1175/JAS-D-18-0146.1
%
% [EDR,CF] = edr_var(X,DR,CUTOFF_SCALE) computes dissipation rate EDR and
% related correction factor CF for a vector X representing wind velocity
% fluctuations, a scalar DR denoting spatial distance between consecutive
% sample points in the record [in meters], typically DR = (true air speed)/
% /(sampling rate), a scalar CUTOFF_SCALE corresponding to the smallest
% spatial scale resolvable by the instrument used [in meters], and a string
% COMPONENT selecting between longitudinal ('lon') and lateral ('lat') 
% velocity component. Typically CUTOFF_SCALE = (true air speed)/(Nyquist freq.)
% but can be larger depending on the measurement and signal acquisition methods.
%
% [EDR,CF] = edr_var(...,'Viscosity',VIS) specifies the kinematic viscosity
% of air [m^2/s]. The default value is 1.506e-5.
%
% [EDR,CF] = edr_var(...,'AbsTol',ABSTOL) controls the absolute tolerance,
% i.e. the absolute difference between the current and the previous
% estimate, at which the iteration is terminated [in m^2/s^3]. The default
% ABSTOL is 1e-7.
%
% [EDR,CF] = edr_var(...,'ModelFunction',MF) specifies the model function
% describing the shape of the energy spectrum in the dissipative range of
% scales. MF needs to be a valid Matlab function handle of single argument.
% The default is the Pope's spectrum given by eq. 6.248 in the textbook 
% "Turbulent flows" by Pope (2000).
%
% [EDR,CF] = edr_var(...,'DissipationSpectrumRange',DSRANGE) controls the
% range of non-dimensionalised wavenumbers k_1* (= k_1*ETA, where ETA is the 
% Kolmogorov microscale) at which the 1-d dissipation spectrum is evaluated
% by integrating the 3-d spectrum function. DSRANGE needs to be 2-element
% vector specifing the upper and lower bounds. The default DSRANGE is 
% [5e-8 1.0] which corresponds to the values of 1-d Pope's dissipation spectrum
% 100-times lower than its maximum. In general, the wider range increases
% the accuracy of CF estimation but slows down the computations.
% 
% [EDR,CF] = edr_var(...,'DissipationSpectrumPoints',DSPTS) controls the
% number of non-dimensionalised wavenumbers k_1* at which the dissipation
% spectrum is evaluated. The default DSPTS is 1000. In general, the higher
% number increases the accuracy at the cost of computation speed.
%
% [EDR,CF] = edr_var(...,'DissipationSpectrum',DSSPECTRUM) allows to
% provide the precomputed dissipation spectrum which can be obtained by
% using DISSIPATION_SPECTRUM function. It is recommended to precompute the 
% dissipation spectrum only once before repeatedly calling the EDR_VAR 
% with the same MF, DSRANGE and DSPTS. DSSPECTRUM should be 2-column matrix
% where the first column contains the non-dimensiolised wavenumbers k_1*
% and the second column contains the corresponding values of the 1-d
% dissipation spectrum.
%
%   Example for longitudinal component:
%       [k1,D,~] = dissipation_spectrum;
%       [edr,cf] = process_timeseries(x,ind,@edr_var,dr,cutoff_scale,'lon',...
%            'DissipationSpectrum',[k1,D]);
%
%   Example for lateral component:
%       [k1,~,D] = dissipation_spectrum;
%       [edr,cf] = process_timeseries(x,ind,@edr_var,dr,cutoff_scale,'lat',...
%            'DissipationSpectrum',[k1,D]);
% 
% [EDR,S] = edr_var(...,'Plot',PLT) selects whether to show diagnostics
% plots. PLT can be true or false (default).
% 
% [...,FIG1,FIG2] = edr_var(...,'Plot',true) provides the handles to the
% diagnostics plots
%
% See also EDR_SFC, EDR_PSD, DISSIPATION_SPECTRUM, REYNOLDS_DECOMPOSITION


arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    cutoff_scale (1,1) {mustBePositive, mustBeFinite, mustBeNonempty, mustBeValidScale(cutoff_scale,x,dr)} % ~ TAS/Nyquist
    component (1,1) string {mustBeMember(component,{'lon','lat'})}
    options.Viscosity (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 1.506e-5
    options.AbsTol (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 1e-7
    options.ModelFunction (1,1) {mustBeA(options.ModelFunction,'function_handle')} = ...
        @(a) exp( -5.2 * ( (a.^4+0.4^4).^0.25 - 0.4 ) )  % Pope (6.248)
    options.DissipationSpectrumRange (1,2) {mustBePositive, mustBeFinite, mustBeNonempty} = [5e-8 1.0] % non-dim k1* (k1*eta)
    % Pope's spectrum decays w.r.t. maximum 100-times at k1* = 5.7e-8 and 1.0, 1000-times at k1* = 5.7e-11 and 1.4
    options.DissipationSpectrumPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 1e3
    options.DissipationSpectrum (:,2) {mustBePositive, mustBeFinite} = [] % [non-dim k1* (k1*eta), non-dim D11* (D11/u_eta^3)]
    options.Plot (1,1) logical = false
end

vis = options.Viscosity;


% Calculate non-dim 1-d dissipation spectrum D11* or D22*

if isempty(options.DissipationSpectrum)
    if strcmp(component,'lat')
        [k1,~,D] = dissipation_spectrum( options.DissipationSpectrumRange,...
            options.DissipationSpectrumPoints, options.ModelFunction );
    else
        [k1,D,~] = dissipation_spectrum( options.DissipationSpectrumRange,...
            options.DissipationSpectrumPoints, options.ModelFunction );
    end
else
    k1 = options.DissipationSpectrum(:,1);
    D  = options.DissipationSpectrum(:,2);
end


% Iterate to estimate the correction factor by integrating dissipation
% spectrum

if strcmp(component,'lat')
    e0 = 15/2*vis*mean(diff(x).^2)/dr^2;
else
    e0 =   15*vis*mean(diff(x).^2)/dr^2;
end

e1 = 1e-5;
e2 = e0;
ev = [e1; e2];

while abs(e2-e1) > options.AbsTol
    
    eta = (vis^3/e2)^0.25;
    k1_cut = 2*pi/cutoff_scale*eta;
    cf = correction_factor(k1,D,k1_cut);
    
    e1 = e2;
    e2 = e0*cf;
    
    if options.Plot
        ev = [ev; e2];   
    end
    
end
edr = e2;


% Diagnostic plots

if options.Plot
    
    % Plot 1: dissipation spectrum
    
    eta = (vis^3/edr)^0.25;
    k1_cut = 2*pi/cutoff_scale*eta;
    k1_cut_v = 2*pi/cutoff_scale*(vis^3./ev).^0.25;
    
    [fig1,~,co] = fig16x12('loglog',[1 1]);
    plot(k1,D,'Color',co(1,:))
    plot(k1_cut_v,interp1(k1,D,k1_cut_v),'o','Color',co(2,:),'MarkerFaceColor',co(2,:))
    plot(k1_cut,interp1(k1,D,k1_cut),'o','Color',co(3,:),'MarkerFaceColor',co(3,:))
    xlabel('$k_1\eta$','Interpreter','latex')
    if strcmp(component,'lat')
        ylabel('$(\epsilon\nu)^\frac{3}{4}D_{22}(k_1)$','Interpreter','latex')
    else
        ylabel('$(\epsilon\nu)^\frac{3}{4}D_{11}(k_1)$','Interpreter','latex')
    end
    text(0.05,0.90,['$\eta = ',sprintf('%.1f',eta*1000),'\,\mathrm{mm}$',newline,...
        '$k_1^{cut}\eta = ', sprintf('%.2f',k1_cut/10^floor(log10(k1_cut))),'\cdot10^',...
        sprintf('{%d}',floor(log10(k1_cut))),'$'],...
        'FontSize',12,'Units','Normalized',...
        'HorizontalAlignment','left','VerticalAlignment','top',...
        'Interpreter','latex')
    
    % Plot 2: correction factor iterations
    
    [fig2,~,co] = fig16x12('linlog',[1 1]);
    plot(1:length(ev),ev,'o','Color',co(2,:),'MarkerFaceColor',co(2,:))
    plot(length(ev),edr,'o','Color',co(3,:),'MarkerFaceColor',co(3,:))
    xlabel('$n$','Interpreter','latex')
    ylabel('$\epsilon\,\mathrm{[m^2 s^{-3}]}$','Interpreter','latex')
    text(0.66,0.10,['$\epsilon = ',sprintf('%.2f',edr/10^floor(log10(edr))),'\cdot10^',...
        sprintf('{%d}',floor(log10(edr))),'\,\mathrm{m^2\,s^{-3}}$',...
        newline, '$C_F = ', sprintf('%5.2f',cf),'$'],...
        'FontSize',12,'Units','Normalized',...
        'HorizontalAlignment','left','Interpreter','latex')
    
else
    fig1 = [];
    fig2 = [];
end

    
end


function CF = correction_factor (x,y,xcut)

ycut = interp1(x,y,xcut,'linear');

indR = (x>xcut);
indL = (x<xcut);

C1 = trapz( [xcut; x(indR)], [ycut; y(indR)] );
C2 = trapz( [x(indL); xcut], [y(indL); ycut] );

CF = 1+C1/C2;

end


function mustBeValidScale(a,x,dr)
    if ~ge(a,dr*2)
        eid = 'Cutoff:tooLow';
        msg = sprintf('Cutoff scale must be within [dr*2 dr*length(x)/2] = [%.2f %.2f].',dr*2,dr*length(x)/2);
        throwAsCaller(MException(eid,msg))
    end
    if ~le(a,length(x)*dr/2)
        eid = 'Cutoff:tooHigh';
        msg = sprintf('Cutoff scale must be within [dr*2 dr*length(x)/2] = [%.2f %.2f].',dr*2,dr*length(x)/2);
        throwAsCaller(MException(eid,msg))
    end
end


function mustBeA(a,cls)
    if ~isa(a,cls)
        eid = 'Class:wrong';
        msg = sprintf('Input must be valid %s.',cls);
        throwAsCaller(MException(eid,msg))
    end
end