function [ strucutre ] = UFactoryStrucutre( x, xp, Edof, bc, varargin )
%UFACTORYSTRUCTURE Generates a struct representind a strucutre with
%design parameters, element degrees of freedom and boundary conditions.
%Optionally other parameters can be added.
%   Detailed explanation goes here

% LAST MODIFIED: A Sehlstrom    2013-05-23
% Copyright (C)  A Sehlstrom

parseo = inputParser;
parseo.CaseSensitive = true;

% Input ckecks
issquare    = @(x) (isnumeric(x) && (size(x,1) == size(x,2)) && ~isscalar(x));
isrowvector = @(x) (isnumeric(x) && (size(x,2) == 1));

% Required
addRequired(parseo,'x',    isrowvector);
addRequired(parseo,'xp');
addRequired(parseo,'Edof', @isnumeric);
addRequired(parseo,'bc');

% Optional
addParamValue(parseo,'F',    NaN, isrowvector);
addParamValue(parseo,'K',    NaN, issquare);
addParamValue(parseo,'Ke0',  NaN, issquare);
addParamValue(parseo,'E',    NaN, isrowvector);
addParamValue(parseo,'E0',   NaN, isrowvector);
addParamValue(parseo,'Emin', NaN, @isscalar);
addParamValue(parseo,'dE',   NaN, isrowvector);
addParamValue(parseo,'M',    NaN, issquare);
addParamValue(parseo,'Me0',  NaN, issquare);
addParamValue(parseo,'m',    NaN, isrowvector);
addParamValue(parseo,'dm',   NaN, isrowvector);
addParamValue(parseo,'Ve0',  NaN, isrowvector);
addParamValue(parseo,'V0',   NaN, @isscalar);
addParamValue(parseo,'Ex',   NaN, @isnumeric);
addParamValue(parseo,'Ey',   NaN, @isnumeric);
addParamValue(parseo,'Ez',   NaN, @isnumeric);

% Parse
parseo.parse(x, xp, Edof, bc, varargin{:});

% Create strucutre
strucutre = ...
    struct('x',    parseo.Results.x, ...
           'xp',   parseo.Results.xp, ...
           'Edof', parseo.Results.Edof, ...
           'bc',   parseo.Results.bc, ...
           'F',    parseo.Results.F, ...
           'K',    parseo.Results.K, ...
           'Ke0',  parseo.Results.Ke0, ...
           'E',    parseo.Results.E, ...
           'E0',   parseo.Results.E0, ...
           'Emin', parseo.Results.Emin, ...
           'dE',   parseo.Results.dE, ...
           'M',    parseo.Results.M, ...
           'Me0',  parseo.Results.Me0, ...
           'm',    parseo.Results.m, ...
           'dm',   parseo.Results.dm, ...
           'Ve0',  parseo.Results.Ve0, ...
           'V0',   parseo.Results.V0, ...
           'Ex',   parseo.Results.Ex, ...
           'Ey',   parseo.Results.Ez, ...
           'Ez',   parseo.Results.Ey);

end

