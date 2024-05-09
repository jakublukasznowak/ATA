%% Turbulence in trade wind cumulus clouds
%% Introduction
% Shallow cumulus clouds in the trade-wind regions cool the planet by reflecting 
% solar radiation. The response of trade cumulus clouds to climate change is a 
% key uncertainty in climate projections. The properties of those clouds and their 
% coupling with mesoscale circulations in the atmosphere were investigated by 
% the EUREC4A measurement campaign in winter 2020 (<https://essd.copernicus.org/articles/13/4067/2021/ 
% Stevens et al., 2021>). The French SAFIRE ATR42 research aircraft (ATR) sampled 
% the atmosphere around the cloud base level following a fixed flight pattern 
% of 120 x 20 km rectangles (<https://essd.copernicus.org/articles/14/2021/2022/ 
% Bony et al., 2022>). Among others, the measurements involved turbulence (<https://essd.copernicus.org/articles/13/3379/2021/ 
% Brilouet et al., 2021>) and microphysical parameters. The flight tracks were 
% divided into segments and the respective data was published in a number of datasets. 
% 
% This example shows how to use the functions of the Atmospheric Turbulence 
% Analysis toolbox in a consistent workflow to analyse measurements performed 
% in the atmosphere with a research aircraft. It is achieved using the data collected 
% in the segment R2B of ATR flight RF12 during EUREC4A.
%% Import data
% Add the toolbox functions to MATLAB path.

addpath(genpath(['..',filesep]))
% Download data files
% Two Net-CDF data files serve as the input:
%% 
% * |EUREC4A_ATR_turbulent_fluctuations_20200205_RF12_R2B_L3_v1.9.nc| containing 
% the timeseries of turbulent wind velocity, temperature and humidity measured 
% in the selected segment from the dataset <https://doi.org/10.25326/128 Brilouet, 
% P. & Lothon, M. (2020)>,
% * |EUREC4A_ATR_PMA_Composite-CDP-2DS_20200205_F12_v1.nc| containg microphysical 
% parameters and derived products (e.g. cloud mask) measured in the whole flight 
% from the dataset <https://doi.org/10.25326/237 Coutris, P. (2021)>.

file_turb = 'EUREC4A_ATR_turbulent_fluctuations_20200205_RF12_R2B_L3_v1.9.nc';
http_turb = 'https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/TURB_FLUCTUATIONS/L3/v1.9/longlegs/RF12/';

file_mcph = 'EUREC4A_ATR_PMA_Composite-CDP-2DS_20200205_F12_v1.nc';
http_mcph = 'https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/PMA/PROCESSED/CloudComposite/';

path_turb = websave(file_turb,[http_turb,file_turb]);
path_mcph = websave(file_mcph,[http_mcph,file_mcph]);
% Load data from files
% Read start and end time of the segment. Load detrended 25 Hz timeseries of 
% turbulent temperature, humidity and wind velocity components: longitudinal, 
% transversal, vertical. For details about the data, see <https://essd.copernicus.org/articles/13/3379/2021/ 
% Brilouet et al. (2021)>.

fsamp_turb = 25; % [Hz]

time_turb_start = ncreadatt(path_turb,'/','time_leg_start');
time_turb_end   = ncreadatt(path_turb,'/','time_leg_end');

time_turb = ncread(path_turb,'time'); % time [seconds from 2020-01-01]
T = ncread(path_turb,'T_DET');  % temperature
R = ncread(path_turb,'MR_DET'); % water vapor mixing ratio
U = ncread(path_turb,'UL_DET'); % longitudinal wind velocity
V = ncread(path_turb,'VT_DET'); % transverse   wind velocity
W = ncread(path_turb,'W_DET');  % vertical     wind velocity
%% 
% Load 1 Hz timeseries of liquid water content and cloud mask for the whole 
% flight.

fsamp_mcph = 1; % [Hz]

time_mcph = ncread(path_mcph,'time'); % time [seconds from 2020-01-01]
LWC = ncread(path_mcph,'LWC');  % liquid water content
CLD = ncread(path_mcph,'CLOUD_mask'); % cloud mask
%% 
% Specify fixed true air speed of the aircraft. The true measured TAS is indeed 
% nearly constant in horizontal segments. The respective data on exact TAS can 
% be found in the dataset <https://doi.org/10.25326/298 CNRM/TRAMM, SAFIRE, Laboratoire 
% d'Aérologie (2021)>.

TAS = 100; % [m/s]
% Arrange the imported data
% Convert time stamps and time vectors into native Matlab |datetime|.

epoch = datetime('2020-01-01 00:00:00.000');

time_turb = datetime(time_turb,'ConvertFrom','epochtime','Epoch',epoch,...
    'Format','yyyy-MM-dd HH:mm:ss.SS','TimeZone','UTC');
time_mcph = datetime(time_mcph,'ConvertFrom','epochtime','Epoch',epoch,...
    'Format','yyyy-MM-dd HH:mm:ss.SS','TimeZone','UTC');

time_turb_start = datetime(dateshift(time_turb(1),'start','day')+duration(time_turb_start),...
    'Format','yyyy-MM-dd HH:mm:ss.SS','TimeZone','UTC');
time_turb_end  = datetime(dateshift(time_turb(1),'start','day')+duration(time_turb_end),...
    'Format','yyyy-MM-dd HH:mm:ss.SS','TimeZone','UTC');
%% 
% Cut out only the selected segment from the whole-flight microphysics data.

ind1 = find(time_mcph >= time_turb_start,1,'first');
ind2 = find(time_mcph <= time_turb_end,  1,'last');

time_mcph = time_mcph(ind1:ind2);
LWC = LWC(ind1:ind2);
CLD = CLD(ind1:ind2);
%% 
% Convert cloud mask from logical vector into a list of time stamps denoting 
% the start and end of cloud penetrations.

CLDind = time_mcph(mask2ind(CLD))+duration(0,0,0.5/fsamp_mcph);
% Plot the imported data

figure('Units','normalized','Position',[0 0 0.6 0.3])
hold on, grid on
co = get(gca,'ColorOrder');
plot(datenum(time_turb),R/3+2,'Color',co(4,:))
plot(datenum(time_turb),T+1,'Color',co(2,:))
plot(datenum(time_mcph),LWC*3,'Color',co(1,:),'LineWidth',1)
axis tight
yl = get(gca,'YLim');
for k=1:size(CLDind,1)
    p=patch(datenum(CLDind(k,[1 2 2 1])),[yl(1) yl(1) yl(2) yl(2)],...
        co(1,:),'FaceAlpha',0.4,'EdgeColor','none');
end
datetick('x','HH:MM')
axis tight
legend({'r/4+2','T+1','LWC*4','clouds'})
xlabel('Time')
title('ATR RF12 R2B thermodynamics')
figure('Units','normalized','Position',[0 0 0.6 0.3])
hold on, grid on
plot(datenum(time_turb),[U V W])
axis tight
yl = get(gca,'YLim');
for k=1:size(CLDind,1)
    p=patch(datenum(CLDind(k,[1 2 2 1])),[yl(1) yl(1) yl(2) yl(2)],...
        co(1,:),'FaceAlpha',0.4,'EdgeColor','none');
end
datetick('x','HH:MM')
axis tight
legend({'u','v','w','clouds'})
xlabel('Time')
title('ATR RF12 R2B wind velocity')
%% Reynolds decomposition
% Reynolds averaing/decomposition is the partition of any signal $x\left(t\right)$into 
% slowly varying mean component $\left\langle x\left(t\right)\right\rangle$ and 
% rapidly fluctuating component $x^{\prime } \left(t\right)\ldotp$
% 
% $$x\left(t\right)=\left\langle x\left(t\right)\right\rangle +x^{\prime } \left(t\right)$$
% 
% Decompose vertical velocity signal using |reynolds_decomposition| function 
% for a cutoff scale of 500 m with two available methods: |movmean| (sliding average) 
% and |butter| (Butterworth filter) selecting the filter order of 6.

cutoff_scale = 500; % m
window_Re = cutoff_scale/TAS*fsamp_turb; % # points

[Wp_movmean,Wm_movmean] = reynolds_decomposition(W,window_Re,'Method','movmean');
[Wp_butter, Wm_butter ] = reynolds_decomposition(W,window_Re,'Method','butter','FilterOrder',6);
%% 
% Compare the mean components for the two methods with the original signal.

figure('Units','normalized','Position',[0 0 0.6 0.3])
hold on, grid on
plot(datenum(time_turb),[W Wm_movmean Wm_butter])
datetick('x','HH:MM')
axis tight
xlabel('Time [s]'), ylabel('w [m/s]')
title('Reynolds decomposition of vertical wind velocity')
legend({'orginal','movmean','butter'})
%% 
% Compare the spectra of the signal components obtained with the two methods.

[psd_W,fv] = pwelch(W,[],[],[],fsamp_turb);
psd_Wm_movmean = pwelch(Wm_movmean,[],[],[],fsamp_turb);
psd_Wm_butter  = pwelch(Wm_butter, [],[],[],fsamp_turb);

figure
hold on, grid on
plot(fv,[psd_W,psd_Wm_movmean,psd_Wm_butter])
plot(ones(2,1)*TAS/cutoff_scale,[min(psd_W) max(psd_W)],'LineWidth',2)
set(gca,'XScale','log','YScale','log')
axis tight
xlabel('Frequency [Hz]'), ylabel('PSD')
legend({'orginal','movmean','butter','cutoff'},'Location','southwest')
title('Spectral properties of decomposed signals')
%% 
% Decompose all turbulent signals using the |butter| method.

[Tp,Tm] = reynolds_decomposition(T,window_Re,'Method','butter','FilterOrder',6);
[Rp,Rm] = reynolds_decomposition(R,window_Re,'Method','butter','FilterOrder',6);
[Up,Um] = reynolds_decomposition(U,window_Re,'Method','butter','FilterOrder',6);
[Vp,Vm] = reynolds_decomposition(V,window_Re,'Method','butter','FilterOrder',6);
[Wp,Wm] = reynolds_decomposition(W,window_Re,'Method','butter','FilterOrder',6);
%% Integral length scales
% Integral length scale of a turbulent signal$x\left(r\right)$is defined as 
% the integral of the autocorrelation function.
% 
% $$\rho \left(r\right)=\frac{\left\langle x\left(r^{\prime } \right)x\left(r^{\prime 
% } +r\right)\right\rangle }{\left\langle {x\left(r^{\prime } \right)}^2 \right\rangle 
% \;}\;$$
% 
% $$L=\int_0^{\infty } \rho \left(r\right)\;\textrm{dr}$$
% 
% In practice, indefinite integration cannot be achieved in the case of measurement 
% data. Therefore, the |integral_lengthscale| function implements two standard 
% practical methods used in such a situation:
%% 
% * |e-decay:| evaluates a distance at which autocorrelation function $\rho 
% \left(r\right)$decays e-times,
% * |integration:| numerically integrate the autocorrelation function up to 
% its first zero.
%% 
% Compute spatial distance between sample points.

dr = TAS/fsamp_turb;
%% 
% Compute the integral length scale for vertical wind velocity obtained with 
% two available methods. Show the autocorrelation function in diagnostics plots.

[L_W_edecay,~] = int_ls_short(W,'Method','e-decay','dr',dr,'Plot',true);
L_W_integration = int_ls_short(W,'Method','integrate','dr',dr,'Plot',true);
%% 
% Compute integral lengthscales for other variables and compare the results 
% of the two methods. Note that the choice of cutoff scale in Reynolds decomposition 
% strongly affects the computed integral length scale if the fluctuating part 
% $x^{\prime } \left(t\right)$is used as the input (<https://journals.ametsoc.org/view/journals/atsc/79/10/JAS-D-22-0028.1.xml 
% Waclawczyk et al. 2022>). For this reason, it is advisable to asses the integral 
% length scale using the full (not-decomposed) signal or perform the decomposition 
% with considerably large cutoff scale (subject to the length of the record and 
% expected large-scale effects).

ILS = table('Size',[5 2],'VariableTypes',{'double','double'},...
    'VariableNames',{'e-decay','integrate'},'RowNames',{'T','R','U','V','W'});

for i=1:size(ILS,1)                     % iterate over turbulent variables (=rows of the table)
    var = ILS.Properties.RowNames{i};   % current variable name
    ILS{var,'e-decay'}     = int_ls_short(eval(var),'Method','e-decay','dr',dr);
    ILS{var,'integrate'} = int_ls_short(eval(var),'Method','integrate','dr',dr);
end

ILS
%% Turbulent moments
% Compute the variance (2nd moment) of turbulent wind velocity using |turb_moment| 
% function.

u2 = turb_moment(Up,2);
v2 = turb_moment(Vp,2);
w2 = turb_moment(Wp,2);
%% 
% Calculate turbulence kinetic energy (TKE).

TKE = 0.5*(u2+v2+w2)
%% 
% Compute the 2nd, 3rd and 4th moment of turbulent signals with |turb_moment| 
% function. Estimate systematic and random sampling errors of the outcome according 
% to <https://journals.ametsoc.org/view/journals/atot/11/3/1520-0426_1994_011_0661_hlilew_2_0_co_2.xml 
% Lenschow et al. (1994)>. Use the integral length scales calculated before as 
% the input for error estimation scheme.

MOM = table('Size',[5 3*3],'VariableTypes',repmat({'double'},1,3*3),...
    'VariableNames',{'second','second_err_sys','second_err_ran',...
                     'third', 'third_err_sys', 'third_err_ran',...
                     'fourth','fourth_err_sys','fourth_err_ran'},...
    'RowNames',{'T','R','U','V','W'});

for i=1:size(MOM,1)                     % iterate over turbulent variables (=rows of the table)
    var = MOM.Properties.RowNames{i};   % current variable name
    for k=2:4                           % iterate over moments 2nd to 4th
        c = 3*(k-2)+1;                  % column index of the kth moment value
        [MOM{var,c},MOM{var,c+1},MOM{var,c+2}] = turb_moment(eval([var,'p']),k,'Lx',ILS{var,'integrate'});
    end
end

MOM
%% 
% The error estimations are fractions relative to the flux value.
% 
% Compute standard deviation, skewness and kurtosis.

MOM.std = sqrt(MOM.second);
MOM.skewness = MOM.third ./ MOM.second.^1.5;
MOM.kurtosis = MOM.fourth ./ MOM.second.^2;

MOM(:,{'std','skewness','kurtosis'})
%% Turbulent fluxes with errors
% Calculate turbulent fluxes of heat $\left\langle w^{\prime } \;T^{\prime } 
% \right\rangle \;$and moisture $\left\langle w^{\prime } r^{\prime } \right\rangle$ 
% with |turb_moment| function. Estimate systematic and random sampling errors 
% of the outcome according to <https://journals.ametsoc.org/view/journals/atot/11/3/1520-0426_1994_011_0661_hlilew_2_0_co_2.xml 
% Lenschow et al. (1994)>.
% 
% First, estimate the parameters which are needed in flux error computation:
%% 
% * integral length scale of the product of the two signals,
% * integral length scale for the cross-correlatation between the two signals,
% * correlation coefficient of the two signals.

ILS2 = table('Size',[2 3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'product','cross','corrcoef'},'RowNames',{'WT','WR'});

ILS2{'WT','product'} = int_ls_short(W.*T-mean(W.*T),'dr',dr);
ILS2{'WR','product'} = int_ls_short(W.*R-mean(W.*R),'dr',dr);

ILS2{'WT','cross'} = int_ls_short(W,T,'dr',dr);
ILS2{'WR','cross'} = int_ls_short(W,R,'dr',dr);

r1 = corrcoef(W,T); ILS2{'WT','corrcoef'} = r1(1,2);
r2 = corrcoef(W,R); ILS2{'WR','corrcoef'} = r2(1,2);

ILS2
%% 
% Compute vertical fluxes of heat and moisture together with their systematic 
% and random sampling errors.

FLX = table('Size',[2 3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'flux','err_sys','err_ran'},'RowNames',{'WT','WR'});

[FLX{'WT','flux'},FLX{'WT','err_sys'},FLX{'WT','err_ran'}] = ...
    turb_moment(Wp,Tp,'Lff',ILS2{'WT','product'},'Lxy',ILS2{'WT','cross'},'Rxy',ILS2{'WT','corrcoef'});
[FLX{'WR','flux'},FLX{'WR','err_sys'},FLX{'WR','err_ran'}] = ...
    turb_moment(Wp,Rp,'Lff',ILS2{'WR','product'},'Lxy',ILS2{'WR','cross'},'Rxy',ILS2{'WR','corrcoef'});

FLX
%% 
% The error estimations are fractions relative to the flux value.
%% TKE dissipation rate
% Estimate turbulence kinetic energy dissipation rate with the 3 methods implemented 
% in the toolbox:
%% 
% * structure function method (|edr_sfc|) is based on the Kolmogorov scaling 
% (2/3) of the 2nd order structure function in the inertial range of scales,
% * power spectrum method (|edr_psd|) is based on the Kolmogorov scaling (-5/3) 
% of the power spectrum in the inertial range of scales,
% * iterative method of corrected variance of derivatives (|edr_var|) assumes 
% specific functional form of the spectrum in the dissipative range of scales 
% and was outlined in sec. 2.c. of <https://journals.ametsoc.org/view/journals/atsc/76/5/jas-d-18-0146.1.xml 
% Akinlabi et al. (2019)>.
% Find the inertial range
% Prior to calculating the TKE dissipation rate, find out what is the range 
% of scales where one can expect the Kolmogorov scaling characteristic for the 
% inertial casacade. Evaluate the power spectrum with |pwelch| using its default 
% parameters and look for the linear trend in log-log plot. 

[psd_Up,fv] = pwelch(Up,[],[],[],fsamp_turb);
psd_Vp = pwelch(Vp,[],[],[],fsamp_turb);
psd_Wp = pwelch(Wp,[],[],[],fsamp_turb);

figure, hold on, grid on
plot(TAS./fv,[psd_Up,psd_Vp,psd_Wp])
set(gca,'XScale','log','YScale','log','XDir','reverse')
axis tight
xlabel('r [m]'), ylabel('PSD')
legend({'u','v','w'},'Location','northwest')
title('Spectrum of turbulent velocity fluctuations')
% Structure function method
% Select the fitting range. In principle, the lower limit in the case of structure 
% function method can be as small as the distance between sample points |dr|. 
% For more information on the best choice of the fitting range, see <https://www.mdpi.com/2073-4433/11/2/199 
% Waclawczyk et al. (2020)>.

sfc_fit_range = [10 100]; % [m]
%% 
% Estimate TKE dissipation rate with |edr_sfc| function using its default method 
% |logmean| with 12 fitting points. This method evaluates structure function at 
% all possible displacements which belong to the fitting range but averages those 
% values in a given number of log-equally-distributed bins covering the fitting 
% range before the fit is performed. Such averaging prevents from the large number 
% of points concentrated on the side of the larger scales to determine the fit 
% at the cost of small number of points on the side of the smaller scales. Alternatively, 
% the fit can be performed without any prior averaging using |direct| option or 
% the structure function can be evaluated only at the displacements corresponding 
% to the log-equally-distributed fitting points using |sparse| options. For details, 
% see the documentation of |edr_sfc| function.
% 
% In addition to dissipation rate, inspect the slope of the line fitted in log-log 
% coordinates (equivalent to the structure function scaling exponent). Remember 
% to declare whether the input signal corresponds to longitudinal or lateral direction. 
% Show the diagnostics plot. 

EDR = table('Size',[3 3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'SFC','PSD','VAR'},'RowNames',{'U','V','W'});
SLP = table('Size',[3 2],'VariableTypes',{'double','double'},...
    'VariableNames',{'SFC','PSD'},'RowNames',{'U','V','W'});

[EDR{'U','SFC'},SLP{'U','SFC'}] = edr_sfc(Up,dr,sfc_fit_range,2.0,'FitPoints',12,'Plot',true);
[EDR{'V','SFC'},SLP{'V','SFC'}] = edr_sfc(Vp,dr,sfc_fit_range,2.6,'FitPoints',12,'Plot',true);
[EDR{'W','SFC'},SLP{'W','SFC'}] = edr_sfc(Wp,dr,sfc_fit_range,2.6,'FitPoints',12,'Plot',true);
% Power spectrum method
% Select the fitting range. In principle, the lower limit in the case of power 
% spectrum method can be as small as the Nyquist limit, i.e. the smallest resolvable 
% wavelength which is typically |2*dr|. For more information on the best choice 
% of the fitting range, see <https://www.mdpi.com/2073-4433/11/2/199 Waclawczyk 
% et al. (2020)>.

psd_fit_range = [10 100]; % [m]
%% 
% Estimate TKE dissipation rate with |edr_psd| function using its default method 
% |logmean| with 16 fitting points. This method evaluates the power spectrum within 
% the maximum available range of normalized frequencies, then averages the values 
% which belong to the fitting range in a given number of log-equally-distributed 
% bins covering the fitting_range and performs a linear fit in log-log coordinates 
% using the averaged points. Analogously to |edr_sfc|, the fit can be performed 
% without any prior averaging using |direct| option. The |sparse| option is also 
% available but not recommended due to its poor performance. For details, see 
% the documentation of |edr_psd| function.
% 
% In addition to dissipation rate, inspect the slope of the line fitted in log-log 
% coordinates (equivalent to the power spectrum scaling exponent). Declare whether 
% the input signal corresponds to longitudinal or lateral direction. Show the 
% diagnostics plot. 

[EDR{'U','PSD'},SLP{'U','PSD'}] = edr_psd(Up,dr,psd_fit_range,0.5,'FitPoints',16,'Plot',true);
[EDR{'V','PSD'},SLP{'V','PSD'}] = edr_psd(Vp,dr,psd_fit_range,0.65,'FitPoints',16,'Plot',true);
[EDR{'W','PSD'},SLP{'W','PSD'}] = edr_psd(Wp,dr,psd_fit_range,0.65,'FitPoints',16,'Plot',true);
% Corrected variance of derivatives method
% Select the cutoff scale up to which the turbulent fluctuations are well resolved 
% in the measured signal. In principle, the lower limit is the smallest resolvable 
% wavelength |2*dr| but can be larger depending on the measurement and signal 
% acquisition methods (<https://www.mdpi.com/2073-4433/11/2/199 Waclawczyk et 
% al., 2020>).

var_cutoff_scale = 2*dr; % [m]
%% 
% Estimate TKE dissipation rate with |edr_var| function using the default Pope 
% dissipation spectrum. Declare whether the input signal corresponds to longitudinal 
% or lateral direction. Show the diagnostics plots: normalized dissipation spectrum 
% and iteration series. The first plot informs what part of the dissipations spectrum 
% is resolved in the signal (left from the dot) and what part is modeled (right 
% from the dot). In both plots, the red dots denote the iterations while the yellow 
% dot denotes the final estimation. For details about the method, see <https://journals.ametsoc.org/view/journals/atsc/76/5/jas-d-18-0146.1.xml 
% Akinlabi et al. (2019)>.

EDR{'U','VAR'} = edr_var(Up,dr,var_cutoff_scale,'lon','Plot',true);
EDR{'V','VAR'} = edr_var(Vp,dr,var_cutoff_scale,'lat','Plot',true);
EDR{'W','VAR'} = edr_var(Wp,dr,var_cutoff_scale,'lat','Plot',true);
% Compare the results
% Compare the results of different methods applied to the three wind velocity 
% components: dissipation rates obtained with the three methods and the fitted 
% slopes for the structure function and power spectrum methods. Remember the expected 
% theoretical values of those slopes are 2/3 and -5/3, respectively. Deviations 
% might signal non-equilibrium or non-stationary conditions (<https://acp.copernicus.org/articles/21/10965/2021/ 
% Nowak et al. 2021>, <https://journals.ametsoc.org/view/journals/atsc/79/10/JAS-D-22-0028.1.xml 
% Waclawczyk et al. 2022>).

EDR
SLP
%% TKE dissipation rate variability
% TKE dissipation rate is a microscopic parameter highly variable in space. 
% With the help of Turbulence Toolbox, one can evaluate the variability of any 
% parameter along the sampling track. The functions |define_av_windows| and |process_timeseries| 
% allow to easily apply the same procedures as used above for the whole segment 
% to the short controlable pieces of the signal called _averaging windows_.
% 
% Define the averaging windows corresponding to the spatial extent of 400 m, 
% displaced with respect to each other by half of their length, with the first 
% one beginnig at the first point of the signal. Select the option |discard| to 
% remove the windows which might extend outside the original signal vector; in 
% this way, ensure each window has the same size.

av_scale = 400;                      % m
av_window = av_scale/TAS*fsamp_turb; % # points
av_step = 0.5*av_window;             % # points
L = length(time_turb);               % # points

window_ind = define_av_windows(L,av_window,1,'Step',av_step,...
    'ReferencePosition','front','EndPoints','discard');
%% 
% Inspect the indices defining first four and last four of the averaging windows.

window_ind([1:4,end-3:end],:)
%% 
% Compute the series of dissipation rates with structure function method for 
% three velocity components.

edr_sfc_U_ts = process_timeseries(Up,window_ind, @edr_sfc,...
    dr,sfc_fit_range,2.0,'Method','logmean','FitPoints',12);
edr_sfc_V_ts = process_timeseries(Vp,window_ind, @edr_sfc,...
    dr,sfc_fit_range,2.6,'Method','logmean','FitPoints',12);
edr_sfc_W_ts = process_timeseries(Wp,window_ind, @edr_sfc,...
    dr,sfc_fit_range,2.6,'Method','logmean','FitPoints',12);
%% 
% Prepare the corresponding time vector with values referring to the middle 
% of the averaging windows.

time_ts = time_turb(window_ind(:,1)) + ...
    0.5*( time_turb(window_ind(:,2)) - time_turb(window_ind(:,1)));
%% 
% Plot the obtained results.

figure('Units','normalized','Position',[0 0 0.6 0.3])
hold on, grid on
plot(datenum(time_ts),[edr_sfc_U_ts edr_sfc_V_ts edr_sfc_W_ts])
axis tight
yl = get(gca,'YLim');
for k=1:size(CLDind,1)
    p=patch(datenum(CLDind(k,[1 2 2 1])),[yl(1) yl(1) yl(2) yl(2)],...
        co(1,:),'FaceAlpha',0.4,'EdgeColor','none');
end
datetick('x','HH:MM')
axis tight
legend({'u','v','w','clouds'})
xlabel('Time')
title('ATR RF12 R2B TKE dissipation rate (SFC method)')
%% 
% The same can be analogously done with the two other methods using |edr_psd| 
% and |edr_var|.

edr_psd_U_ts = process_timeseries(Up,window_ind, @edr_psd,dr,psd_fit_range,0.50,'Method','logmean','FitPoints',12);
edr_psd_V_ts = process_timeseries(Vp,window_ind, @edr_psd,dr,psd_fit_range,0.65,'Method','logmean','FitPoints',12);
edr_psd_W_ts = process_timeseries(Wp,window_ind, @edr_psd,dr,psd_fit_range,0.65,'Method','logmean','FitPoints',12);
%% 
% However, in the case of the repetitive application of the |edr_var| with the 
% same options (e.g. in |process_timeseries|), it is recommended to precompute 
% the normalized dissipation spectrum only once and insert it as the parameter 
% to |edr_var|. Otherwise, this computation would be repeated at each call of 
% |edr_var| which significantly slows down the procedure. The example below illustrates 
% the correct implementation.

[k1,Dlon,Dlat] = dissipation_spectrum;
edr_var_U_ts = process_timeseries(Up,window_ind, @edr_var,dr,var_cutoff_scale,'lon','DissipationSpectrum',[k1,Dlon]);
edr_var_V_ts = process_timeseries(Vp,window_ind, @edr_var,dr,var_cutoff_scale,'lat','DissipationSpectrum',[k1,Dlat]);
edr_var_W_ts = process_timeseries(Wp,window_ind, @edr_var,dr,var_cutoff_scale,'lat','DissipationSpectrum',[k1,Dlat]);
%% 
% Compare the results for vertical component between the three methods.

figure('Units','normalized','Position',[0 0 0.6 0.3])
hold on, grid on
plot(datenum(time_ts),[edr_sfc_W_ts edr_psd_W_ts edr_var_W_ts])
axis tight
yl = get(gca,'YLim');
for k=1:size(CLDind,1)
    p=patch(datenum(CLDind(k,[1 2 2 1])),[yl(1) yl(1) yl(2) yl(2)],...
        co(1,:),'FaceAlpha',0.4,'EdgeColor','none');
end
datetick('x','HH:MM')
axis tight
legend({'sfc','psd','var','clouds'})
xlabel('Time')
title('ATR RF12 R2B TKE dissipation rate (vertical velocity)')
% Clouds versus clear air
% It is clearly visible in the plots above that cloudy regions usually coincide 
% with the largest values of TKE dissipation rate. Compare the mean values between 
% cloudy and clear air pieces of the segment.
% 
% First, resample the cloud mask so that it complies with |time_ts| time vector. 
% Then, average the values corresponding to clouds and clear air separately.

CLD_ts = (interp1(time_mcph,CLD,time_ts,'linear')>0.5);

EDR_ts_mean = table('Size',[3 3*2],'VariableTypes',repmat({'double'},1,3*2),...
    'VariableNames',{'sfc cloud','sfc clear',...
                     'psd cloud','psd clear',...
                     'var cloud','var clear'},...
    'RowNames',{'U','V','W'});

for i = 1:size(EDR_ts_mean,1)                   % iterate over velocity components (=rows of the table)
    cmp = EDR_ts_mean.Properties.RowNames{i};   % current variable name
    for j = 1:size(EDR_ts_mean,2)/2             % iterate over methods (sfc,psd,var)
        c = 2*(j-1)+1;                          % auxiliary column index
        met = EDR_ts_mean.Properties.VariableNames{j}(1:3); % method
        x = eval(['edr_',met,'_',cmp,'_ts']);
        EDR_ts_mean{cmp,c}   = mean(x(CLD_ts));
        EDR_ts_mean{cmp,c+1} = mean(x(~CLD_ts));
    end
end

EDR_ts_mean