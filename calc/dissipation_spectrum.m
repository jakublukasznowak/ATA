
function [k1,D11,D22] = dissipation_spectrum (k1_range,Np,f_eta)

% DISSIPATION_SPECTRUM computes the 1-d dissipation spectrum by integrating
% given 3-d energy spectrum function. It can be used to precompute the dissipation
% spectrum before repetitive calculations with EDR_VAR.
%
% [k1,D11,D22] = dissipation_spectrum(K1_RANGE,NP,F_ETA) calculates the
% vector of non-dimensionalised wavenumbers K1 together with the
% corresponding longitudinal D11 and transversal D22 dissipation spectra
% evaluated by integrating the model function F_ETA in the range of
% non-dimensionalised wavenumbers K1_RANGE at NP sampling points.
%
% F_ETA needs to be a valid Matlab function handle of single argument.
% The default is the Pope's spectrum given by eq. 6.248 in the textbook 
% "Turbulent flows" by Pope (2000).
% K1_RANGE needs to be 2-element vector specifing the upper and lower bounds.
% The default K1_RANGE is [5e-8 1.0] which corresponds to the values of 
% the Pope's dissipation spectrum 100-times lower than its maximum.
% The default NP is 1000.


arguments
    k1_range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty} = [5e-8 1.0] % non-dim k1*
    % Pope's spectrum decays w.r.t. maximum 100-times at k1* = 5.7e-8 and 1.0, 1000-times at k1* = 5.7e-11 and 1.4
    Np (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 1e3
    f_eta (1,1) {mustBeA(f_eta,'function_handle')} = ...
        @(a) exp( -5.2 * ( (a.^4+0.4^4).^0.25 - 0.4 ) )  % Pope (6.248)
end


k1 = exp( linspace( log(k1_range(1)), log(k1_range(2)), Np )' );
% k1 = linspace( k1_range(1), k1_range(2), Np )';

g1 = @(a) a.^(-8/3).*f_eta(a);
g2 = @(a) a.^(-14/3).*f_eta(a);

G1 = nan(size(k1));
G2 = nan(size(k1));
for i=1:Np
    G1(i) = integral(g1,k1(i),Inf); % default AbsTol=1e-10, RelTol=1e-6
    G2(i) = integral(g2,k1(i),Inf);
end

D11 = 2*1.5*    (k1.^2.*G1 - k1.^4.*G2);
D22 = 2*1.5*0.5*(k1.^2.*G1 + k1.^4.*G2);

end


function mustBeA(a,cls)
    if ~isa(a,cls)
        eid = 'Class:wrong';
        msg = sprintf('Input must be valid %s.',cls);
        throwAsCaller(MException(eid,msg))
    end
end