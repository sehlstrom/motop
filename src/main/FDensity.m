function [x, dF]=FDensity(H, Hs, x, dF)
%FDENSITY applies the density filter on the design parameters x and
%optionally on the objective sensitivities with respect to x dF.
%   EQUATION
%   Each element e is filtered as:
%
%       xeeTilde = 1/sum(Hei)*sum(Hei*xii)
%
%   where Hei is a convolutionary weight factor defined as:
%
%       Hei = max(0, rmin - dist(ee,ii))
%
%   where rmin is the filter radius and dist(ee,ii) gives the
%   centre-to-centre Euclidian distance between element ee and ii.
%
%   SYNTAX
%   x = FDENSITY(H, Hs, x) 
%   [x, dF] = FDENSITY(H, Hs, x, dF)
%
%   DESCRIPTION
%   FDENSITY applies the density filter on the design parameters x and
%   optionally on the objective sensitivities with respect to x dF.
%
%   The convolutionary weight factor is given as arguments H and Hs to the
%   function; thus the factors can be computed separately and used several
%   times. The weight factors can be found using the function FUsetup.
%
%   INPUT ARGUMENTS
%       H      convolutionary operator; one row for each element ee and one
%              colum for each element ii, thus H(ee,ii) is the
%              concolutionary operator for element e with respect to
%              element ii. H is sparse.
%       Hs     row sums of the convolutionary operator; one row for each
%              element ee.
%       x      vector of design parameters; 0 <= x <= 1
%       dF     objective function first order derivatives with respect to x
%
%   OUTPUT ARGUMENTS
%       x      filtered vector of design parameters; 0 <= x <= 1
%       dF     filtered objective function first order derivatives with
%              respect to x
%
% See also: FUsetup FSensitivity

% LAST MODIFIED: A Sehlstrom    2013-05-20
% Copyright (C)  A Sehlstrom
    
x  = (H * x)./Hs; 

if nargout == 2
    dF = H *(dF./Hs);
end
end