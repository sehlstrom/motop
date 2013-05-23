%% Example 1 - Topology optimization of an MBB beam using the OC-method
% Optimizes the stiffness of half of a so-called MBB beam, see below, under
% a fixed volume fraction constraint. The beam is loaded with a point load
% P at the midpoint of the beam. The beam is modeled by square plane stress
% elasticity elements (Melsoh) in 2D.
%
% MBB beam: simply supported beam:
%
%                       | P
%                       V
%      o---------------------------------o
%      |                                 |
%      o---------------------------------o
%      A                                 o
%     ¯¯¯                               ¯¯¯
%
% MBB FE-discretization:
%
%
%      | P
%      v
%   |8>o--o--o--o--o--o
%      |  |  |  |  |  |
%   |8>o--o--o--o--o--o
%      |  |  |  |  |  |    nely
%   |8>o--o--o--o--o--o
%      |  |  |  |  |  |
%   |8>o--o--o--o--o--o
%                     o
%                    ¯¯¯
%           nelx 

% Make sure motop main functions are added to the search paths
addpath(genpath(fullfile(fileparts(pwd),'main')), '-BEGIN');

% Check if required CALFEM functions are found on the Matlab path
if exist('hooke','file') ~=2 
    error('CALFEM function "hooke" has to be on the Matlab path')
elseif exist('planre','file') ~=2 
    error('CALFEM function "planre" has to be on the Matlab path')
end

%% Section 1 - Define all parameters
% Define all parameters related to mesh, material, boundary conditions,
% initial conditions etc.

% Generate the mesh -------------------------------------------------------
nelx = 32;                               % number of elements in x-direc.
nely = 20;                               % number of elements in y-direc.
l    = 0.01;                             % element side length

[Ndof, Ncoord, Edof, Ex, Ey] = ...       % make the rectangular mesh
    MURectangle2(nelx, nely, l, l);

ndof  = max(max(Ndof));                  % number of degrees of freedom
nelem = size(Edof,1);                    % number of elements

% Boundary conditions -----------------------------------------------------
pdof = union((nelx+1)*2, ...             % lower right uy dof
             1:((nelx+1)*2):ndof);       % all ux on the right edge
         
bc = zeros(length(pdof), 2);             % set pdfo to 0, i.e. roller
bc(:,1) = pdof;

% Element properties ------------------------------------------------------
E0    = 200e9;                           % base material Young's modulus
Emin  = E0/1e9;                          % minimum Young's modulus
nu0   = 0.3;                             % base material Poisson's ratio
rho0  = 7800;                            % base material density
t     = 0.04*nelx*l;                     % element thickness [m]
ep    = [E0, Emin, nu0, rho0, t, l];     % element properties

ptype = 1;                               % plane stress
D     = hooke(1, ptype, nu0);            % unit constitutive matrix

Ke0   = planre([0 l],[0 l],[ptype t],D); % unit element stiffness matrix
                                         % of a square element

% Loading -----------------------------------------------------------------
P = -1000;                               % Point load of 1 kN
F = zeros(ndof,1);                       % initialize load vector
F(ndof-nelx*2) = P;                      % put P in F

% Design parameters -------------------------------------------------------
x0 = ones(nelem,1)*0.4;                  % initial guess; all elements 40 %
xp = [nelx            1;                 % prescribed parameters (leave
      nelx*(nely-1)+1 1];                % empty if none are to be
                                         % prescribed)

% Objective function ------------------------------------------------------
OFun = @OCompliance;                     % compliance objective function

% Volume fraction constraint ----------------------------------------------
vfrac = 0.3;                             % 30 % of the design domain should
                                         % be used

% Interpolation settings --------------------------------------------------
ip1 = {@EModSIMP, 3};                    % ModSIMP interpolation scheme
                                         % with penalization power 3

ip2 = {@ERAMP, 20};                      % RAMP interpolation scheme with
                                         % interpolation parameter 20

% Filter settings ---------------------------------------------------------
ft1 = {@FDensity, 0.015};                % Density filter with filter
                                         % radius rmin = 0.015 m
                                           
ft2 = {@FSensitivity, 0.015};            % Sensitivity filter with filter
                                         % radius rmin = 0.015 m

% Display settings --------------------------------------------------------
d1  = {1, 1, figure(1)};                 % print all and plott all in
                                         % figure 1
                                           
d2  = {1, 1, figure(2)};                 % print all and plott all in
                                         % figure 2

%% Section 2 - Optimization with different interpolation schemes
% Comparison of the the ModSIMP and the RAMP interpolation schems. The
% density filter is used in both examples

% Optimization using the ModSIMP interpolation scheme
x11 = OptOC2( x0, xp, ep, Edof, bc, Ex, Ey, F, Ke0, OFun, vfrac, ip1, ft1, d1 );

% Optimization using the RAMP interpolation scheme
x12 = OptOC2( x0, xp, ep, Edof, bc, Ex, Ey, F, Ke0, OFun, vfrac, ip2, ft1, d2 );

%% Section 3 - Optimization with different filters
% Comparison of the the sensitivity filter and the density filter. The
% ModSIMP interpolation scheme is used in both examples

% Optimization using the density interpolation scheme
x21 = OptOC2( x0, xp, ep, Edof, bc, Ex, Ey, F, Ke0, OFun, vfrac, ip1, ft1, d1 );

% Optimization using the sensitivity interpolation scheme
x22 = OptOC2( x0, xp, ep, Edof, bc, Ex, Ey, F, Ke0, OFun, vfrac, ip1, ft2, d2 );