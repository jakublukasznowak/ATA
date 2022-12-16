
function [F,e_s,e_r,Lx,Ly,Lxy,Lff] = turb_moment (x,y,dr,options)

% TURB_MOMENT calculates statistical moment of turbulent fluctuations of a 
% given order for a single variable or a turbulent flux, i.e. a covariance,
% for two variables. Systematic and random relatice errors of the result are
% estimated with the equations derived in:
%
%   Lenschow, D., Mann, J., and Kristensen, L.: How Long Is Long Enough
%   When Measuring Fluxes and Other Turbulence Statistics?, J. Atmos. Ocean.
%   Tech., 11, 661–673, 1994.
%
% M = turb_moment(X,N,DR) computes the moment M of the order specified by
% an integer N for a vector X representing turbulent fluctuations timeseries
% and a scalar DR denoting spatial distance between consecutive sample points
% in the record [in meters], typically DR = (true air speed)/(sampling rate)
%
% F = turb_moment(X,Y,DR) computes the flux, i.e. the covariance of the
% vectors X and Y which need to be of the same size.
% 
% [...,ES,ER] = turb_moment(...) reports the estimated systematic error ES
% and random error ER of the result according to Lenschow et al. 1994.
%
% [...,ES,ER] = turb_moment(...,'IntegralLengthScaleMethod',METHOD) selects
% the method for the estimation of integral lengthscales used in the error
% estimation which is achieved by calling the function INTEGRAL_LENGTHSCALE. 
% METHOD can be 'e-decay' or 'integration'.
%
% [...,ES,ER] = turb_moment(...,'IntegralLengthSampling',SMPL) controls
% the sampling interval [in meters] for calculating the correlation function
% by calling INTEGRAL_LENGTHSCALE(...,'Sampling',SMPL)
% 
% [...,LX] = turb_moment(X,N,DR) reports the integral lengthscale of X 
% used in error calculation.
%
% [...,LX,LY,LXY,LFF] = turb_moment(X,N,DR) reports the integral lengthscales
% of X, Y, the integral lengthscale obtained based on crosscorrelation, the
% integral lengthscale of the product timeseries X*Y, respectively, which
% are used in error calculation.
%
% [...,ES,ER] = turb_moment(...,'Lx',LX,'Ly',LY,'Lxy',LXY,'Lff',LFF,'Rxy',RXY)
% accepts precomputed integral lengthscales of X, Y, the lengthscale based on
% crosscorrelation function between X and Y, the integral lengthscale of
% X'Y', the correlation coefficient RXY betwenn X and Y, respectively,
% intead of calculating those inside this function.
%
% See also INTEGRAL_LENGTHSCALE, REYNOLDS_DECOMPOSITION


arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    y (:,1) {mustBeReal, mustBeFinite, mustBeNonempty, mustBeVectorOrInteger(y,x)}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 1
    options.IntegralLengthscaleMethod (1,1) string {mustBeMember(options.IntegralLengthscaleMethod,{'e-decay','integration'})} = 'integration'
    options.IntegralLengthscaleSampling (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = dr
    options.Lx (1,1) {mustBeNonnegative, mustBeFinite, mustBeNonempty} = 0
    options.Ly (1,1) {mustBeNonnegative, mustBeFinite, mustBeNonempty} = 0
    options.Lxy (1,1) {mustBeNonnegative, mustBeFinite, mustBeNonempty} = 0
    options.Lff (1,1) {mustBeNonnegative, mustBeFinite, mustBeNonempty} = 0
    options.Rxy (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = 0
end


if isscalar(y)
    
    % MOMENTS 
    
    order = y;
    F = mean(x.^order);
    
    if nargout>1
        
        L = length(x)*dr;
        if options.Lx>0
            Lx = options.Lx;
        else
            Lx = integral_lengthscale(x,dr,'Method',options.IntegralLengthscaleMethod);
        end
        Ly = [];
        Lxy= [];
        Lff= [];
        
        if order==2
            e_s = Lx/L*2;          % L94 (14)
            e_r = sqrt( Lx/L*2 );  % L94 (36)
            
        elseif order==3
            skw = F/mean(x.^2)^1.5;
            
            fskw=@(a) 2*a.*(3+4*a.^2)./(1+2*a.^2).^(3/2);   % L94 (20)
            % fkrt=@(a) 3*(1+20*a.^2+20*a.^4)./(1+2*a.^2).^2; % L94 (20)
            fun=@(x) fskw(x)-skw;
            a = fzero(fun,0.1);
            
            e_s = Lx/L*3*(2-1/(1+a^2)/(3+4*a^2)); % L94 (21)
            e_r = sqrt( Lx/L * ...
                (1+2*a^2)*(1+147*a^2+1476*a^4+780*a^6)/(a^2*(1+a^2)*(3+4*a^2)^2) ); % LMK93 (B40)
            
        elseif order==4
            e_s = 4*Lx/L;              % L94 (17)
            e_r = sqrt( Lx/L*84/9 );   % L94 (38)
            
        else
            fprintf('Error equations for order %d not implemented :(.\n',order)
        end
        
    end
    
else
    
    % FLUXES
    
    F = mean(x.*y);
    
    if nargout>1
        
        L = length(x)*dr;
        
        % First-choice exact expressions for errors

        if options.Lxy>0
            Lxy = options.Lxy;
        else
            Lxy = integral_lengthscale(x,dr,'Method',options.IntegralLengthscaleMethod,...
                'CrossCorrelatedSignal',y);
        end
        if options.Lff>0
            Lff = options.Lff;
        else
            Lff = integral_lengthscale(x.*y-F,dr,'Method',options.IntegralLengthscaleMethod);
        end
        if abs(options.Rxy)>0
            Rxy = options.Rxy;
        else
            Rxy = abs(F)/sqrt(mean(x.^2)*mean(y.^2));
        end
        
        a = L/Lxy;
        e_s_1 = 2/a-2/a^2+2/exp(a)/a^2;         % L94 (30)
        e_r_1 = sqrt( 2*Lff/L*(1+1/Rxy^2) );    % L94 (48) = LS86

        % Alternative estimation of upper limits of the errors
        % (not used as for now)
        
%         if options.Lx>0
%             Lx = options.Lx;
%         else
%             Lx = integral_lengthscale(x,dr,'Method',options.IntegralLengthscaleMethod);
%         end
%         if options.Ly>0
%             Ly = options.Ly;
%         else
%             Ly = integral_lengthscale(y,dr,'Method',options.IntegralLengthscaleMethod);
%         end
% 
%         e_s_2 = 2/abs(Rxy)*sqrt(Lx*Ly)/L;            % L94 (29)
%         e_r_2 = 2/abs(Rxy)*sqrt( min([Lx,Ly])/L );   % L94 (49)
        
        % Select the smallest of the two
        % (not used as for now)
        
%         e_s = min([e_s_1 e_s_2]);
%         e_r = min([e_r_1 e_r_2]);
        e_s = e_s_1;  
        e_r = e_r_1;
        
        
    end
    
end
        

end


function mustBeVectorOrInteger(a,x)

if isscalar(a)
    if ~(round(a) == a)
        eid = 'Order:NotInteger';
        msg = sprintf('The order of the moment has to be scalar integer.\n');
        throwAsCaller(MException(eid,msg))
    end
else
    if ~(length(a) == length(x))
        eid = 'Vector:DifferentLength';
        msg = sprintf('The input vectors need to be of the same length.\n');
        throwAsCaller(MException(eid,msg))
    end
end

end