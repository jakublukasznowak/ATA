
function IND = define_av_windows (signal_length,window_size,querry_points,options)

% DEFINE_AV_WINDOWS prepares indices representing averaging windows which
% can be further used to obtain various turbulence parameters from those
% pieces of the signal.
%
% IND = define_av_windows(SIGNAL_LENGTH,WINDOW_SIZE,QUERRY_POINTS) computes
% a Nx2 matrix of first and last indices of the N averaging windows of the 
% given size positioned at querry points across the length of the original
% signal. SIGNAL_LENGTH and WINDOW_SIZE are the length of the signal and the size
% of the averaging window, respectively, given as a number of sample
% points. QUERRY_POINTS is a vector specifying the positions of the windows
% in the orignal timeseries.
%  
% IND = define_av_windows(...,QUERRY_POINTS,'Step',STEP) conctructs a
% vector of querry points when QUERRY_POINTS is a scalar, by adding the lag
% of STEP sample points to the initial scalar value until SIGNAL_LENGTH is
% reached.
%
% IND = define_av_windows(...,'ReferencePosition',REFPOS) specifies the
% reference point with respect to which the averaging windows are defined.
% REFPOS can be one of the two:
% 
%        'front'      - (default) QUERRY_POINTS refer to the first point
%                       of the averaging window
%        'center'     - QUERRY_POINTS refer to the center of
%                       the averaging window
%
% IND = define_av_windows(...,'EndPoints',ENDPT) controls how the windows
% are handled at the endpoints of the signal where there not enough
% elements. ENPT can be one of the following:
%
%         'shrink'    - (default) the size of the windows is reduced on one
%                       side (non-symetrically) to fit inside the valid
%                       signal
%         'discard'   - the windows which shall effectively extend outside
%                       the original signal are removed from the output IND
%
% See also PROCESS_TIMESERIES, REYNOLDS_DECOMPOSITION


arguments
    signal_length (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty}
    window_size (1,1) {mustBePositive, mustBeFinite, mustBeNonempty}
    querry_points (:,1) {mustBeReal, mustBeFinite, mustBeNonempty} = 1;
    options.Step (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 1
    options.ReferencePosition (1,1) string {mustBeMember(options.ReferencePosition,{'front','center'})} = 'front'
    options.EndPoints (1,1) string {mustBeMember(options.EndPoints,{'shrink','discard'})} = 'shrink'
end


% If QUERRY_POINTS is a scalar, assume it defines the first query point and
% the rest is constructed with step increments
if isscalar(querry_points)
    querry_points = (querry_points:options.Step:signal_length)';
end

% If querries corresponds to the center of the window, move them to the front
if strcmp(options.ReferencePosition,'center')
    querry_points = querry_points - window_size/2;
end


% Create index list
IND = [ceil(querry_points), floor(querry_points+window_size)-1];


% Handle indices from outside the available range
if strcmp(options.EndPoints,'shrink')
    IND(IND(:,1)<1,1)=1;
    IND(IND(:,2)>signal_length,2)=signal_length;
elseif strcmp(options.EndPoints,'discard')
    IND( or(IND(:,1)<1,IND(:,2)>signal_length),: ) = [];
end
    

end

