function [ x, k ] = OptOC2( x0, xp, ep, Edof, bc, Ex, Ey, F, Ke0, OFun, vfrac, ip, ft, d )
%OPTOC2 Topology optimization accoridng to the optimality cirteria (OC)
%method. Mesh consists of equal size square 2D elements.
%
%   SYNTAX
%   x = OPTOC2( x0, xp, ep, Ke0, Edof, bc, Ex, Ey, F, OFun, vfrac )
%   x = OPTOC2( x0, xp, ep, Ke0, Edof, bc, Ex, Ey, F, OFun, vfrac, ip )
%   x = OPTOC2( x0, xp, ep, Ke0, Edof, bc, Ex, Ey, F, OFun, vfrac, ip, ft )
%   x = OPTOC2( x0, xp, ep, Ke0, Edof, bc, Ex, Ey, F, OFun, vfrac, ip, ft, d )
%   [x, k] = OPTOC2(...)
%
%   DESCRITPTION
%   All elements are square elements with sides of length l, nodes,
%   coordinates, and degrees of freedom defined according to the figure
%   below:
%
%       [u7,u8]       [u5,u6]
%       (x4,y4)       (x3,y3)
%        4 o-------------o 3
%          |             |
%          |             |
%          |             |
%          |      e      |  l
%          |             |
%          |             |
%          |             |
%        1 o-------------o 2
%       (x1,y1)        (x2,y2)
%       [u1,u2]   l    [u3,u4]
%
%   If the function should print the result of each iteration, the output
%   will look like:
%
%       It.: 1 Obj.: 1.212e-02 Vol.frac.: 0.80 Ch.: 2.000e-01 N.: 2.576e-02
%
%   where "It." is the iteration number, "Obj." the objective value,
%   "Vol.frac." the current volume fraction, "Ch." the norm of the change
%   in x from the current iteration and the previous, and "N." is the norm
%   of the first order derivative of the objective.
%
%   INPUT ARGUMENTS
%       x0     vector of initial guesses for the design parameters;
%              0 <= x0 <= 1
%       xp     prescribed parameters; one row [e, xp] for each element e
%              which has a prescribed value xp for it's design variable xe.
%              The values are given within the range 0 <= xp <= 1. If no
%              parameters are to be prescribed
%       ep     element properties; ep = [E0, Emin, t, l] where
%              E0       Young's modulus of the base material [Pa]
%              Emin     minimum Young's modulus [Pa]
%              t        element thickness [m]
%              l        element side length [m]
%       Ke0    element stiffness matrix with unit stiffness (E=1)
%       Me0    element mass matrix with unit mass (m=1)
%       Edof   elemnet degrees of freedom; one row for each element
%              [u1,u2,u3,...,un]
%       Ex     element x-coordinates; one row for each element
%              [ex1,ex2,...,exn]
%       Ey     element y-coordinates; one row for each element
%              [ey1,ey2,...,eyn]
%       F      vector of loads
%       vfrac  volume fraction constraint; 0 < vfrac < 1
%       ip     (optional) cell array with interpolation function and
%              related interpolation parameter; ip = {@EFun, param}.
%              If ip is not set, the ELin interplation function will be
%              used.
%       ft     (optional) cell array with filter function and related
%              filter parameter; ft = {@FFun, rmin}
%       d      (optional) cell array with display settings;
%              d = {print, plot, figure}
%              print    1 - if each iteration should print it's result
%                       0 - if no print should be done
%              plot     1 - if the structure should be plotted in each
%                           iteration
%                       0 - if the structure should not be plotted in each
%                           iteration
%              figure   figure handle that should be used when plotting
%
%   OUTPUT ARGUMENTS
%       x      optimized design parameters
%       k      number of iterations
%
% See also: ELin, EModSIMP, ERAMP, FDensity, FSensitivity

% LAST MODIFIED: A Sehlstrom    2013-05-23
% Copyright (C)  A Sehlstrom

% Input checks ------------------------------------------------------------
% Design parameters x
if size(x0,2) ~= 1
    error('OptOC2:argChk', '"x0" must be a row vector');
else
    nelem = size(x0, 1);
end

% Prescribed design paramters x
if size(xp) ~= [0,0]
    felem = setdiff(1:nelem, xp(:,1));
else
    felem = 1:nelem;
end

% Element properties
ep(1,:) = ep;
if size(ep, 2) ~= 4
    error('OptOC2:argChk', '"ep" has to have 4 elements');
else
    E0   = ep(:,1);
    Emin = ep(:,2);
    t    = ep(:,3);
    l    = ep(:,4);
end

% Unit stiffness matrix
if size(Ke0,1) ~= size(Ke0,2)
    error('OptOC2:argChk', '"Ke0" has to be square');
else
    ndofe = size(Ke0,1);
end

% Element degrees of freedom
if size(Edof, 1) ~= nelem
    error('OptOC2:argChk', '"Edof" has to have nelem rows');
elseif size(Edof, 2) ~= ndofe
    error('OptOC2:argChk', '"Edof" has to have ndofe columns');
else
    ndof = max(max(Edof));
end


% Boundary conditions
if size(bc,1) == 0
    error('OptOC2:argChk','"bc" has to have at least 1 row')
elseif size(bc, 2) ~= 2
    error('OptOC2:argChk','"bc" has to have 2 columns')
else
    fdof = setdiff(1:ndof, bc(:,1));
end

% Element coordinates
if size(Ex) ~= size(Ey)
    error('OptOC2:argChk', '"Ex" and "Ey" has to have the same size');
elseif size(Ex, 1) ~= nelem
    error('OptOC2:argChk', '"Ex" and "Ey" has to have nelem rows');
end

% Load vector
if size(F,1) ~= ndof
    error('OptOC2:argChk', '"F" must have same number of degrees of freedom as specified in "Edof"');
elseif size(F,2) ~= 1
    error('OptOC2:argChk', '"F" must be a row vector');
end

% Volume fraction
if size(vfrac) ~= [1,1]
    error('OptOC2:argChk', '"vfrac" must be a scalar');
elseif vfrac <0
    error('OptOC2:argChk', '"vfrac" must be grater or equal than 0');
elseif vfrac > 1
    error('OptOC2:argChk', '"vfrac" must be smaller than or equal to 1');
end

% Additional functions ----------------------------------------------------
% Young's modulus by interpolation setup
if exist('ip','var')
    interp_fun   = ip{1};
    interp_param = ip{2};
else
    interp_fun   = @ELin;
    interp_param = [];
end

% Filter setup
if exist('ft','var')
    filter_fun  = ft{1};
    filter_name = func2str(filter_fun);
    filter_rmin = ft{2};
    ec = [(Ex(:,1)+Ex(:,2))/2, (Ey(:,1)+Ey(:,3))/2];
    [filter_H, filter_Hs] = FUsetup(ec, filter_rmin);
end

% Display setup
if exist('d','var')
    display_print  = d{1};
    display_plot   = d{2};
    display_figure = d{3};
else
    display_print = 0;
    display_plot  = 0;
end

% Initial guess -----------------------------------------------------------
xk = x0;

% Plot initial guess if required
if display_plot
    display_figure = figure(display_figure);
    display_figure = clf(display_figure);
    hold on;
    colormap(flipud(gray));
    for ii = 1:nelem
        fill(Ex(ii,:), Ey(ii,:), xk(ii));
    end
    caxis([0 1]);
    colorbar;
    axis equal;
    axis tight;
    axis off;
    pause(1e-8);
end

% Constraint --------------------------------------------------------------
V0     = l*l*t;
VOmega = nelem*V0;
dV     = ones(length(x0),1)*V0;

% Filter volume derivatives if needed
if exist('filter_name','var') && strcmp(filter_name,'FDensity')
    [~, dV]    = filter_fun(filter_H, filter_Hs, xk, dV);
end

% Limits ------------------------------------------------------------------
% TODO: Make it possible to change these values via the function call
% Stop criteria limits
stop_maxIter = 50;
stop_absTol  = 1e-3;
stop_gradTol = 1e-6;

% Optimality criteria method settings
oc_ltol    = 1e-4;      % Tolerance for lagrangian multiplier
oc_move    = 0.2;       % Move limit for update of x
oc_eta     = 0.5;       % Numerical damping parameter

% Optimize ----------------------------------------------------------------
k = 0;
while 1
    k = k + 1;
    xold = xk;
    
    % INTERPOLATE MATERIAL PROPERTIES
    [E, dE] = interp_fun(xk, E0, Emin, interp_param);
    
    % ASSEMBLE STIFFNESS MATRIX
    K = zeros(ndof, ndof);
    for ii = 1:length(x0)
        K(Edof(ii,:),Edof(ii,:)) = K(Edof(ii,:),Edof(ii,:)) + E(ii)*Ke0;
    end
    
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [O, dO] = OFun(xk, xp, Edof, bc, F, K, Ke0, E, dE);
    
    % APPLY FILTER
    if exist('filter_name','var')
        [~, dO]    = filter_fun(filter_H, filter_Hs, xk, dO);
    end
    
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    l1 = 0; l2 = 100000;
    fac = max(0, -dO)./dV;
    while (l2-l1 > oc_ltol)
        % Lagrangian midvalue guess
        lmid = 0.5*(l2+l1);
        
        % Do update
        x_new = max(0,max(xk-oc_move,min(1.,min(xk+oc_move,xk.*(fac/lmid).^oc_eta))));
        
        % Set prescribed values
        if size(xp) ~= [0,0]
            x_new(xp(:,1)) = xp(:,2);
        end
        
        % Apply filter on x if needed
        if exist('filter_fun','var')
            x_phys = filter_fun(filter_H, filter_Hs, x_new, dO);
        end
        
        % Set prescribed values
        if size(xp) ~= [0,0]
            x_phys(xp(:,1)) = xp(:,2);
        end
        
        % Test optimality
        if sum(sum(x_phys*V0)) - vfrac*VOmega > 0;
            l1 = lmid;
        else
            l2 = lmid;
        end
    end
    
    xk = x_phys;
    
    % PLOT DENSITIES
    if display_plot
        display_figure = figure(display_figure);
        display_figure = clf(display_figure);
        hold on;
        colormap(flipud(gray));
        for ii = 1:nelem
            fill(Ex(ii,:), Ey(ii,:), xk(ii));
        end
        caxis([0 1]);
        colorbar;
        axis equal;
        axis tight;
        axis off;
        pause(1e-8);
    end
    
    % PRINT RESULTS
    % Find change
    c = max(max(abs(xk-xold)));
    
    % Find norm
    n = norm(dO(felem));
    
    % Print
    if display_print
        V = sum(xk*V0,1);
        fprintf('It.: %4i Obj.: %3.3e Vol.frac.: %1.2f Ch.: %3.3e N.: %3.3e \n', k, O, V/VOmega, c, n);
    end
    
    % STOPPING CRITERIA
    % Check criteiras and return when the first criteria is met
    if k >= stop_maxIter
        x = xk;
        if display_print
            display('Stop criteria: max iterations');
        end
        return;
    elseif c < stop_absTol
        x = xk;
        if display_print
            display('Stop criteria: change tolerance');
        end
        return;
    elseif n < stop_gradTol
        x = xk;
        if display_print
            display('Stop criteria: gradient tolerance');
        end
        return;
    end
end

end
