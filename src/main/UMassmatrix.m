function [ Me ] = UMassmatrix( lx, ly, t, rho )
%UMASSMATRIX is a utility function that computes the element mass matrix of
%a rectangular element using 2x2 Gauss integration scheme.
%
%    SYNTAX
%    Me=UMASSMATRIX(lx, ly, t, rho)
%
%    DESCRIPTION
%    UMASSMATRIX computes the element mass matrix of a rectangular element
%    with sides lx and ly, thickness t and material density rho. The
%    computation is done using 2x2 Guass integration scheme.
%
%    INPUT ARGUMENTS
%       lx     element length, x-direction
%       ly     element length, y-direction
%       t      element thickness
%       rho    material density
%
%   OUTPUT ARGUMENTS
%       Me     element mass matrix; size(Me) = 8x8
% 

% LAST MODIFIED: A Sehlstrom    2013-08-16
% Copyright (C)  A Sehlstrom

% Half elements side lenthts
a = lx / 2;
b = ly / 2;

% Element artificial coordinates
x1 = -a;   y1 = -b;
x2 =  a;   y2 = -b;
x3 =  a;   y3 =  b;
x4 = -a;   y4 =  b;

% Build integration scheme
ir = 2;
% 2x2 Gauss points
if ir == 2
    % Gauss points
    xgp = [-1  1  1 -1] / sqrt(3);
    ygp = [-1 -1  1  1] / sqrt(3);

    % ... and corresponding weights
    w(:,1) = [ 1  1  1  1];    % x
    w(:,2) = [ 1  1  1  1];    % y
end
wp=w(:,1).*w(:,2);

% Initialize
Me = zeros(8,8);

% Apply integration
for ii=1:length(xgp)
    % Gauss point
    x = xgp(ii) * a;
    y = ygp(ii) * b;
    
    % Shape functions @ Gauss point
    N1 = (x-x2)*(y-y4);
    N2 = -(x-x1)*(y-y3);
    N3 = (x-x4)*(y-y2);
    N4 = -(x-x3)*(y-y1);
    
    % Shape function matrix
    N = [ N1 0  N2 0  N3 0  N4 0  ;
          0  N1 0  N2 0  N3 0  N4 ] / (4*a*b);
    
    % Assemble
    Me = Me + wp(ii)*(N')*N*a*b*rho*t;
end