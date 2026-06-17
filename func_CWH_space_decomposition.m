function [obsBoxes, freeBoxes] = func_CWH_space_decomposition(doPlot)
%ISS_POLYTOPE_DECOMPOSITION
%
% Creates a simple polytope decomposition of a bounded 3D workspace:
%
%
% The spacecraft is represented as a union of simple axis-aligned boxes.
% The free space is decomposed into non-overlapping axis-aligned boxes.
%
% Usage:
%
%   [obsPolys, freePolys, obsBoxes, freeBoxes] = func_CWH_space_decomposition(true);
%
% Inputs:
%
%   doPlot : optional boolean, default true
%
% Outputs:
%
%   obsBoxes  : Nobs x 6 array, rows are [xmin xmax ymin ymax zmin zmax]
%   freeBoxes : Nfree x 6 array, rows are [xmin xmax ymin ymax zmin zmax]
%   domain        : struct with domain.lb and domain.ub
%
%

    if nargin < 1 || isempty(doPlot)
        doPlot = true;
    end



    % ---------------------------------------------------------------------
    % Bounding domain for the decomposition.
    % This is the region of space that will be decomposed.
    % Units can be interpreted as meters, but this is only a rough model.
    % ---------------------------------------------------------------------
    domain.lb = [-70, -45, -25];
    domain.ub = [ 70,  45,  25];

    % ---------------------------------------------------------------------
    % Build an ISS-like obstacle from simple boxes.
    % Each row is [xmin xmax ymin ymax zmin zmax].
    % ---------------------------------------------------------------------
    obsBoxes = makeISSLikeSpacecraftBoxes();

    % ---------------------------------------------------------------------
    % Generate non-overlapping free-space boxes.
    % ---------------------------------------------------------------------
    freeBoxes = decomposeFreeSpaceIntoBoxes(domain, obsBoxes);

    fprintf('Obstacle polytopes:   %d\n', size(obsBoxes,1));
    fprintf('Free-space polytopes: %d\n', size(freeBoxes,1));

    % ---------------------------------------------------------------------
    % Plot if requested.
    % ---------------------------------------------------------------------
    if doPlot
        func_plotISSDecomposition(obsBoxes, freeBoxes, domain);
    end
end


function obs = makeISSLikeSpacecraftBoxes()
%MAKEISSLIKESPACECRAFTBOXES
%
% Returns a rough ISS-like spacecraft model as a collection of boxes.
% Rows are [xmin xmax ymin ymax zmin zmax].
%
% This is intentionally simple: truss, pressurized modules, solar arrays,
% radiators, and docking vehicles are all approximated by rectangular prisms.

    obs = [
        % Main long truss
        -55   55    -1.2   1.2   -1.2   1.2;

        % Central pressurized modules
         -8    8   -19    19     -2.8   2.8;

        % Central node / hub
         -5    5    -5     5     -5.0   5.0;

        % Docked vehicles or extended modules
         -4    4   -27   -19     -2.4   2.4;
         -4    4    19    27     -2.4   2.4;

        % Upper and lower protrusions
         -6    6    -4     4      3.0   7.0;
        -12   12    -4     4     -6.0  -3.0;

        % Solar array wings, positive z side
        -52  -42   -13    13      2.2   3.0;
        -40  -30   -13    13      2.2   3.0;
        -25  -15   -13    13      2.2   3.0;
         15   25   -13    13      2.2   3.0;
         30   40   -13    13      2.2   3.0;
         42   52   -13    13      2.2   3.0;

        % Small masts connecting truss to arrays
        -48  -46    -1.2   1.2    1.2   2.2;
        -36  -34    -1.2   1.2    1.2   2.2;
        -21  -19    -1.2   1.2    1.2   2.2;
         19   21    -1.2   1.2    1.2   2.2;
         34   36    -1.2   1.2    1.2   2.2;
         46   48    -1.2   1.2    1.2   2.2;

        % Radiator-like panels on negative z side
        -50  -35   -10    10     -4.0  -3.2;
         35   50   -10    10     -4.0  -3.2;

		  % Small masts connecting truss to Radiator-like panels
        -45  -40    -1.2   1.2    -3.2   -1.2;
         40   45    -1.2   1.2    -3.2   -1.2;


        % % Additional side equipment boxes
        % -15  -10     4    10     -2.0   2.0;
        %  10   15     4    10     -2.0   2.0;
        % -15  -10   -10    -4     -2.0   2.0;
        %  10   15   -10    -4     -2.0   2.0;
    ];

	% obsCopy = obs;
	% obs = obsCopy(:,[5 6 3 4 1 2]);
end


function freeBoxes = decomposeFreeSpaceIntoBoxes(domain, obsBoxes)
%DECOMPOSEFREESPACEINTOBOXES
%
% Creates a box decomposition of the domain minus the union of obstacle boxes.
%
% Method:
%
%   1. Use all obstacle x/y/z bounds to create an axis-aligned grid.
%   2. Mark grid cells as occupied or free.
%   3. Greedily merge neighboring free cells into larger boxes.
%
% The resulting free boxes have disjoint interiors and cover the bounded
% free space exactly with respect to the axis-aligned obstacle model.

    tol = 1e-10;

    % Build grid breakpoints from domain and obstacle boundaries.
    xEdges = [domain.lb(1), domain.ub(1)];
    yEdges = [domain.lb(2), domain.ub(2)];
    zEdges = [domain.lb(3), domain.ub(3)];

    for i = 1:size(obsBoxes, 1)
        lo = max(obsBoxes(i, [1 3 5]), domain.lb);
        hi = min(obsBoxes(i, [2 4 6]), domain.ub);

        if all(lo < hi)
            xEdges = [xEdges, lo(1), hi(1)]; %#ok<AGROW>
            yEdges = [yEdges, lo(2), hi(2)]; %#ok<AGROW>
            zEdges = [zEdges, lo(3), hi(3)]; %#ok<AGROW>
        end
    end

    xEdges = uniqueTol(xEdges, tol);
    yEdges = uniqueTol(yEdges, tol);
    zEdges = uniqueTol(zEdges, tol);

    nx = numel(xEdges) - 1;
    ny = numel(yEdges) - 1;
    nz = numel(zEdges) - 1;

    freeMask = false(nx, ny, nz);

    % Classify each grid cell by its center.
    % Because obstacle boundaries are included in the grid, this exactly
    % classifies each open cell as free or occupied.
    for ix = 1:nx
        cx = 0.5 * (xEdges(ix) + xEdges(ix+1));
        for iy = 1:ny
            cy = 0.5 * (yEdges(iy) + yEdges(iy+1));
            for iz = 1:nz
                cz = 0.5 * (zEdges(iz) + zEdges(iz+1));
                p = [cx, cy, cz];

                if ~pointInsideAnyBox(p, obsBoxes, tol)
                    freeMask(ix, iy, iz) = true;
                end
            end
        end
    end

    % Greedily merge free cells into larger axis-aligned boxes.
    freeBoxes = func_hyperrectangles_greedyMerge(freeMask, xEdges, yEdges, zEdges);

    % Additional exact face-merging pass.
    freeBoxes = func__hyperrectangles_mergeAdjacent(freeBoxes, tol);
end




function inside = pointInsideAnyBox(p, boxes, tol)
%POINTINSIDEANYBOX
%
% Returns true if point p lies inside any box in boxes.

    if isempty(boxes)
        inside = false;
        return;
    end

    insideEach = ...
        p(1) >= boxes(:,1) - tol & p(1) <= boxes(:,2) + tol & ...
        p(2) >= boxes(:,3) - tol & p(2) <= boxes(:,4) + tol & ...
        p(3) >= boxes(:,5) - tol & p(3) <= boxes(:,6) + tol;

    inside = any(insideEach);
end


function v = uniqueTol(v, tol)
%UNIQUETOL
%
% Sorts and uniquifies a vector using a tolerance.

    v = sort(v(:).');

    if isempty(v)
        return;
    end

    keep = [true, abs(diff(v)) > tol];
    v = v(keep);
end