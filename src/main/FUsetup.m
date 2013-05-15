function [H, Hs]=FUsetup(ec, rmin)
%FUSETUP is a filter utility function that computes the convolutionary
%factor H and it's row sums Hs used in mesh independend filters.
%    EQUATION
%    The convolutionary operator is determined according to
%
%        H(ee, ii)  = max(0, rmin-r)
%    
%    where r is the distance between the centroid of element ee and ii. The
%    row sums are computed as
%
%        Hs = sum(H,2)
%
%    SYNTAX
%    [H, Hs]=FUsetup(ec, rmin)
%
%    DESCRIPTION
%    FUSETUP computes the linear convolutionary operator and it's row sums
%    used in a mesh independend filter. For element ee, the filter area
%    (2D) or volume (3D) is bounded by the perimiter of a circle or sphear
%    with radius rmin. Thus, only the properties of elements ii that has
%    their centroid whitin the filter area/volume will be taken into
%    consideration by the filter; outside the filter area, the filter has
%    no effect.
%
%    Within the filter area/volume the effect of the filter is linearly
%    decreasing with the distance between element ee and the surounding
%    elements inside the filter area/volume.
%
%    INPUT ARGUMENTS
%       ec     element centoid coordinates; each row in ec corresponds to
%              the centroid coordinates [x,y] in 2D or [x,y,z] in 3D of an
%              element ee.
%       rmin   filter radius; for each element ee, the filter will operate
%              on the elements ii whose centroid is within the circle with
%              radius rmin in 2D or the sphear with radius rmin in 3D where
%              the centerpoint of the filter area/volume is placed in the
%              centroid of element ee.
%
%   OUTPUT ARGUMENTS
%       H      convolutionary operator; one row for each element ee and one
%              colum for each element ii, thus H(ee,ii) is the
%              concolutionary operator for element e with respect to
%              element ii. H is sparse.
%       Hs     row sums of the convolutionary operator; one row for each
%              element ee.
%
% See also: Fdensity Fsensitvity
% 

% LAST MODIFIED: A Sehlstrom    2013-05-15
% Copyright (C)  A Sehlstrom

[rows,cols] = size(ec);

H  = zeros(rows);

if cols == 2
    for ee = 1:rows
        % Centroid, element e
        xpe  = ec(ee, 1);
        ype  = ec(ee, 2);
        
        for ii = 1:rows
            % Centroid, element i
            xpi = ec(ii, 1);
            ypi = ec(ii, 2);
            
            % Distance between e and i
            r  = sqrt((xpe-xpi)^2+(ype-ypi)^2);
            
            % Convolutionary operator
            H(ee, ii)  = max(0, rmin-r);
        end
    end
elseif cols == 3
    for ee = 1:rows
        % Centroid, element e
        xpe  = ec(ee, 1);
        ype  = ec(ee, 2);
        zpe  = ec(ee, 3);
        
        for ii = 1:rows
            % Centroid, element i
            xpi = ec(ii, 1);
            ypi = ec(ii, 2);
            zpi = ec(ii, 3);
            
            % Distance between e and i
            r  = sqrt((xpe-xpi)^2+(ype-ypi)^2+(zpe-zpi)^2);
            
            % Convolutionary operator
            H(ee, ii)  = max(0, rmin-r);
        end
    end
else
    error('FUsetup:UnsuportedDimension','The centroid coordinate has to be defined in 2D or 3D')
end

% Make H sparse
% most of the entries will be 0 since rmin is usually much smaller than the
% overall size of the considered FE mesh.
H = sparse(H);

% Find the row sums of H
Hs = sum(H,2);
end