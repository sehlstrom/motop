function [E, dE, ddE] = EModSIMP(x, E0, Emin, p)
%EMODSIMP Compute Young's modulus according to the modified SIMP
%interpolation scheme
%   EQUATION
%   Computes Young's modulus according to
%
%       E = (Emin + x.^p * (E0-Emin))
%
%   SYNTAX
%   E = EMODSIMP(x, E0, Emin, p)
%   [E, dE] = EMODSIMP(...)
%   [E, dE, ddE] = EMODSIMP(...)
%
%   DESCRIPTION
%   EMODSIMP computes Young's modulus according to the modified SIMP
%   interpolation scheme and it's first and second derivative with respect
%   to the design parameters x.
%
%   INPUT ARGUMENTS
%       x      vector of design parameters; 0 <= x <= 1
%       E0     base Young's modulus; for x = 1, the output E = E0.
%       Emin   minimum stiffness; set to > 0 in order to avoid
%              singularities when x = 0 or small.
%       p      interpolation constant
%
%   OUTPUT ARGUMENTS
%       E      vector of Young's modulus corresponding to the parameters
%              given in x
%       dE     vector of first order derivatives of E with respect to x
%       ddE    vector of second order derivatives of E with respect to x
%
% See also: Elin Eramp

% LAST MODIFIED: A Sehlstrom    2013-05-14
% Copyright (C)  A Sehlstrom

if length(varargin) < 4
    error('Emodsimp:argChk', '4 or more inputs needed')
end

x    = varargin{1};
E0   = varargin{2};
Emin = varargin{3};
p    = varargin{4};

E = (Emin + x.^p * (E0-Emin));

dE = p*x.^(p-1) * (E0-Emin);

ddE = (p-1).*p*x.^(p-2) * (E0-Emin);

end
