function [E, dE, ddE] = ERAMP(varargin)
%ERAMP Compute Young's modulus according to the RAMP
%interpolation scheme
%   EQUATION
%   Computes Young's modulus according to
%
%       E = Emin + x./( 1+ q.*(1-x) ) * (E0-Emin)
%
%   SYNTAX
%   E = ERAMP(x, E0, Emin, q)
%   [E, dE] = ERAMP(...)
%   [E, dE, ddE] = ERAMP(...)
%
%   DESCRIPTION
%   ERAMP computes Young's modulus according to the RAMP
%   interpolation scheme and it's first and second derivative with respect
%   to the design parameters x.
%
%   INPUT ARGUMENTS
%       x      vector of design parameters; 0 <= x <= 1
%       E0     base Young's modulus; for x = 1, the output E = E0.
%       Emin   minimum stiffness; set to > 0 in order to avoid
%              singularities when x = 0 or small.
%       q      interpolation constant
%
%   OUTPUT ARGUMENTS
%       E      vector of Young's modulus corresponding to the parameters
%              given in x
%       dE     vector of first order derivatives of E with respect to x
%       ddE    vector of second order derivatives of E with respect to x
%
% See also: ELin EModSIMP

% LAST MODIFIED: A Sehlstrom    2013-05-23
% Copyright (C)  A Sehlstrom

if nargin < 4
    error('Eramp:argChk', '4 or more inputs needed')
end

x    = varargin{1};
E0   = varargin{2};
Emin = varargin{3};
q    = varargin{4};

E = Emin + x./( 1+ q.*(1-x) ) * (E0-Emin);

dE = (1+q)./((1+q*(1-x)).^2) * (E0-Emin);

ddE = 2*q * (q+1)./(1 + q*(1-x)).^3 *(E0-Emin);

end
