
function varargout = process_timeseries(x,ind,fun,varargin)

% PROCESS_TIMESERIES processes a signal by aplying the same function to the
% given pieces of the input timeseries.
% 
% [...] = process_timeseries(X,IND,FUN,...) applies the function FUN to each
% piece of the vector X, where those pieces are defined by a Nx2 matrix of 
% window indices IND. IND is typically generated using DEFINE_AV_WINDOWS.
%
% The following input arguments are forwarded to the function FUN and
% should be given in an order required by the specific function of choice.
% The output arguments are the same as the outputs of the function FUN but
% concanetenated over all the processed pieces of the signal X.
%
% [...] = process_timeseries(...,'EstimateTime',ET) selects whether to
% estimate the elapsed and remaining time needed to finish the calculation.
% ET is the number of seconds which defines the frequency of time
% estimation monit. ET equal to 0 disables this option.
% 
% Examples:
%
%   [EDR,S] = process_timeseries(X,IND,@EDR_SFC,DR,FITTING_RANGE) computes
%   a vector of eddy disipation rate values EDR and structure function 
%   scaling exponents inside the windows given by IND using the method 
%   implemented in the function EDR_SFC for the given parameters DR and 
%   FITTING_RANGE.
%
%   TIME_IND = process_timeseries(TIME,IND,@mean) can be used to generate
%   from the orignal time vector TIME the time vector TIME_IND which is 
%   consistent with the other results obtained for the same windows IND.
%
% See also DEFINE_AV_WINDOWS, REYNOLDS_DECOMPOSITION


arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    ind (:,2) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty}
    fun
end

arguments (Repeating)
    varargin
end


Nw = size(ind,1); % number of windows
out = cell(Nw,nargout);


% Check if EstimateTime option is on
argTime = cellfun(@(x) strcmp(x,'EstimateTime'),varargin);
indTime = find(argTime,1,'first');
if ~isempty(indTime)
    everySeconds = varargin{indTime+1};
    if everySeconds>0
        us = etd(clock,0,Nw,everySeconds);
        ifTime = true;
    else
        ifTime = false;
    end
    varargin(indTime:indTime+1) = [];
else
    ifTime = false;
end


% Iterate over windows
for i=1:Nw
    [out{i,1:nargout}] = fun( x(ind(i,1):ind(i,2)), varargin{:} );
    if ifTime, us = etd(us,i); end
end


% Arrange output
varargout = cell(1,nargout);
for j=1:nargout
    varargout{j} = cell2mat(out(:,j));
end


end
