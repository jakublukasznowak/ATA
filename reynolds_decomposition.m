
function [xp,xm] = reynolds_decomposition (x,window_size,options)

% REYNOLDS_DECOMPOSITION decomposes a signal into slowly-varying mean and 
% fluctuations. 
% 
% [XP,XM] = reynolds_decomposition(X,WINDOW_SIZE) for a vector X, representing
% a measured timeseries, and a positive integer WINDOW_SIZE, defining the cutoff
% scale as a number of sample points, computes two vectors, XP and XM, 
% respresenting rapidly fluctuating and slowly-varying components of the input
% signal, such that X = XM + XP.
%
% [XP,XM] = reynolds_decomposition(...,'Method',METHOD) specifies the method of
% decomposition. METHOD can be one of the following:
%  
%        'movmean'    - (default) XM is a centered moving average computed 
%                       by sliding a window of length WINDOW_SIZE along X.
%        'butter'     - XM is a lowpass filtered X. The zero-phase digital
%                       filtering is performed with the Butterworth filter
%                       of the normalized cutoff 1/(2*WINDOW_SIZE) by
%                       processing the input data x in both the forward and
%                       reverse directions.
% 
% [XP,XM] = reynolds_decomposition(...,'FilterOrder',ORDER) specifies
% filter order (positive integer) when METHOD is 'butter'. The default is 10.
%
% [XP,XM] = reynolds_decomposition(...,'EndPoints',ENDPT) controls how 
% XM is calculated at the endpoints of X, where there are not enough 
% elements to fill the window when METHOD is 'movmean'. ENDPT can be one of
% the following:
% 
%         'shrink'    - (default) compute XM over the number of elements 
%                       of X that are inside the window, effectively
%                       reducing the window size to fit X at the endpoints.
%         'discard'   - compute XM  only when the window is filled with
%                       elements of X, discarding partial endpoint
%                       calculations and their corresponding elements in XM.
%                       This truncates the output.
% 
% See also DETREND, MOVMEAN, BUTTER, FILTFILT


arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    window_size (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty}
    options.Method (1,1) string {mustBeMember(options.Method,{'movmean','butter'})} = 'movmean'
    options.EndPoints (1,1) string {mustBeMember(options.EndPoints,{'discard','shrink'})} = 'shrink'
    options.FilterOrder (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 10
end

    
if strcmp(options.Method,'movmean')
    xm = movmean(x,window_size,'EndPoints',options.EndPoints);
    margin = (length(x)-length(xm))/2;
    xp = x( ceil(margin)+1:end-floor(margin) ) - xm;
    
elseif strcmp(options.Method,'butter')
    [z,p,k] = butter(options.FilterOrder,2/window_size);
    sos = zp2sos(z,p,k);
    xm = filtfilt(sos,1,x);
    xp = x - xm;
    
else
    error('Invalid method selected.')
    
end
        

end