function [ V, dV ] = ObjVolume( s )
%OBJVOLUME calculates the volume V and the first order volume
%derivative dV with respect to x for a given design x.
%
%   SYNTAX
%       V = OBJVOLUME( s)
%       [V,dV] = OBJVOLUME( s )
%
%   DESCRIPTION
%       OBJVOLUME computes the colume as:
%
%           V = sum(x*Ve0)
%
%       and the first order derivative of V with respect to x as:
%
%           dV = Ve0
%
%       where Ve0 is the base volume of each element.
%
%   INPUT ARGUMENTS
%       s      a struct with at least the following fields
%              x      vector of design parameters, size(x) = (nelem x 1) 
%                     where nelem is the number of elements; 0 <= x0 <= 1
%              xp     prescribed parameters matrix, size(xp) = (npx x 2) 
%                     where npx is the number of prescribed x parameters; 
%                     each row in xp is [e, pv] where e is the element 
%                     number that is prescribed and 0 <= pv <= 1 is the
%                     prescribed value
%              Ve0    base volume of each elements; either Ve0 is a scalar
%                     or Ve0 is a vector with one entry per element
%
%   OUTPUT ARGUMENTS
%       V      volume
%       dV     sensitivity (first order derivative of V with respect to the
%              parameters x)

% LAST MODIFIED: A Sehlstrom    2013-05-23
% Copyright (C)  A Sehlstrom

% Parse inputs ------------------------------------------------------------
parseo = inputParser;
addRequired(parseo,'s', @isstruct);
parseo.parse(s);

% Extract
x   = s.x;
xp  = s.xp;
Ve0 = s.Ve0;

% Input checks ------------------------------------------------------------
if size(x,2) ~= 1
    error('OCompliance:argChk', '"x" must be a row vector');
end

if isnan(Ve0)
    error('OCompliance:argChk', 'Element base volume "Ve0" in struct "s" has to be set');
end

% Setup -------------------------------------------------------------------
% Make sure prescribed values are prescribed
if size(xp) ~= [0,0]
    x(xp(:,1)) = xp(:,2);
end

% Determin objective ------------------------------------------------------
% Initialize
if size(Ve0, 1) == 1
    V = sum(x*Ve0);
    if nargout > 1
        dV = ones(length(x), 1)*Ve0;
    end
else
    V = sum(x*Ve0);
    if nargout > 1
        dV = Ve0;
    end
end


end

