function [E, dE, ddE] = Elin(x, E0, Emin, q)
%ELIN Compute Young's modulus according to a linear
%interpolation scheme
%   EQUATION
%   Computes Young's modulus according to
%
%       E = E0*x
%
%   SYNTAX
%   E = ELIN(x, E0, Emin, q)
%   [E] = ELIN(...)
%   [E, dE] = ELIN(...)
%   [E, dE, ddE] = ELIN(...)
%
%   DESCRIPTION
%   ELIN computes Young's modulus according to a linear
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
% See also: Emodsimp Eramp

% LAST MODIFIED: A Sehlstrom    2013-05-14
% Copyright (C)  A Sehlstrom

E = x.*E0;

dE = E0;

ddE = 0;

end