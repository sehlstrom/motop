function [Ndof, Ncoord, Edof, Ex, Ey] = MURectangle2(nelx, nely, lx, ly)
%MURECTANGLE2 is a mesh utility function that generates a rectangluar mesh
%of equally sized 2D 4-node elemnets with 2 degrees of freedom per node.
%The size of the mesh is determined by the number of elements in x- and y-
%direction and by the side lengths lx and ly of the elements.
%
%   DESCRIPTION
%       MURECTANGLE2 genereates a rectangular FE-mesh of (nelx x nely)
%       rectangluar elemets of size (lx x ly). The mesh will be generated
%       according to the figure below
%
%           o---o---o---o---o...o---o    -+
%           |   |   |   |   |   |   |     |
%           o---o---o---o---o...o---o     |
%           :   :   :   :   :   :   :     |
%           o---o---o---o---o...o---o    nely
%           |   |   |   |   |   |   |     |
%           o---o---o---o---o...o---o     |
%           |e=1| 2 | 3 | 4 |   |   |     |
%           o---o---o---o---o...o---o    -+
%          n=1  2   3   4   5
%
%           |                       |
%           +--------- nelx --------+
%
%       where each element has the following definitions
%
%       [ux4,uy4]       [ux3,uy3]
%        (x4,y4)         (x3,y3)
%         4 o---------------o 3
%           |               |
%           |               |
%           |               |
%           |       e       |  ly
%           |               |
%           |               |
%           |               |
%         1 o---------------o 2
%        (x1,y1)          (x2,y2)
%       [ux1,ux1]   lx   [ux2,uy2]
%
%   INPUT ARGUMENTS
%       nelx   number of elements in x-direction
%       nely   number of elements in y-direction
%       lx     element lenght in x-direction
%       ly     element lenght in y-direction
%
%   OUTPUT ARGUMENTS
%       Ndof   matrix with nodal degrees of freedom; each row represents
%              the degrees of freedom of the node [ux uy]
%       Ncoord matrix with nodal coordinates; each row represents one node
%              [x y]
%       Edof   matrix with element degrees of freedom; each row represents
%              the degrees of freedom of the element
%              [ux1 uy1 ux2 uy2 ux3 uy3 ux4 uy4]
%       Ex     matrix with element x-coordinates; each row represents the
%              element x-coordinates [x1 x2 x3 x4]
%       Ex     matrix with element y-coordinates; each row represents the
%              element y-coordinates [y1 y2 y3 y4]

% LAST MODIFIED: A Sehlstrom    2013-05-22
% Copyright (C)  A Sehlstrom
    
% Node degrees of freedom and coordinates: Ndof and Ncoord
ni     = 1;
dofi   = 1;
Ndof   = zeros((nely+1)*(nelx+1),2);
Ncoord = zeros((nely+1)*(nelx+1),2);
nodes  = zeros(nely+1, nelx+1);

for ny = 1:nely+1
    for nx = 1:nelx+1
        Ndof(ni,:)      = [dofi  dofi+1];
        Ncoord(ni,:)    = [(nx-1)*lx (ny-1)*ly];
        nodes(ny, nx)   = ni;
        
        ni              = ni + 1;
        dofi            = dofi + 2;
    end
end

% Generate element degrees of freedom: Edof
ei = 1;
Edof = zeros(nelx*nely, 8);
for ely = 1:nely
    for elx = 1:nelx
        edof1         = Ndof(nodes(ely,   elx),   :);
        edof2         = Ndof(nodes(ely,   elx+1), :);
        edof3         = Ndof(nodes(ely+1, elx+1), :);
        edof4         = Ndof(nodes(ely+1, elx),   :);
        
        Edof(ei, :)   = [edof1 edof2 edof3 edof4];
        
        ei            = ei + 1;
    end
end

% Generate element coordinates: Ex and Ey
ei = 1;
Ex = zeros(nelx*nely, 4);
Ey = zeros(nelx*nely, 4);
for ely = 1:nely
    for elx = 1:nelx
        Ex(ei, :) = [elx*lx-lx, elx*lx, elx*lx, elx*lx-lx];
        Ey(ei, :) = [ely*ly-ly, ely*ly-ly, ely*ly, ely*ly];
        ei            = ei + 1;
    end
end

end