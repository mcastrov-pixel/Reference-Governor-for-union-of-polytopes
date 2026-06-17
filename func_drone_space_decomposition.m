function [buildingBoxes, freeBoxes, domain] = ...
	func_drone_space_decomposition(nBuildings, doPlot, rngSeed)
%CITY_POLYTOPE_DECOMPOSITION
%
% Generates a bounded 3D space decomposition with:
%
%   1. buildingBoxes/buildingPolys:
%      Axis-aligned boxes representing buildings.
%
%   2. freeBoxes/freePolys:
%      Non-overlapping axis-aligned boxes covering the free space inside
%      the domain, excluding building interiors.
%
% Buildings are visually plotted with roofs and windows, but geometrically
% each building is represented as a rectangular prism.
%
% Inputs:
%
%   nBuildings : number of buildings to generate
%   doPlot     : true/false, whether to plot
%   rngSeed    : random seed for repeatability
%
% Outputs:
%
%   buildingBoxes : N x 6 array, [xmin xmax ymin ymax zmin zmax]
%   freeBoxes     : M x 6 array, [xmin xmax ymin ymax zmin zmax]
%   domain        : struct with domain.lb and domain.ub
%
% Notes:
%
%   - The free-space decomposition is bounded by domain.
%   - Buildings start at z = 0.
%   - Touching at shared faces/edges/corners is allowed.
%   - Free-space boxes have non-overlapping interiors.

if nargin < 1 || isempty(nBuildings)
	nBuildings = 12;
end

if nargin < 2 || isempty(doPlot)
	doPlot = true;
end

if nargin < 3 || isempty(rngSeed)
	rngSeed = 1;
end


rng(rngSeed,'threefry');

% ---------------------------------------------------------------------
% Define bounded city domain.
% ---------------------------------------------------------------------
domain.lb = [0,   0,   0];
domain.ub = [120, 90,  80];

% ---------------------------------------------------------------------
% Generate random non-overlapping buildings.
% ---------------------------------------------------------------------
buildingBoxes = generateBuildings(nBuildings, domain);

% ---------------------------------------------------------------------
% Decompose free space.
% ---------------------------------------------------------------------
freeBoxes = decomposeFreeSpaceIntoBoxes(domain, buildingBoxes);

fprintf('Generated buildings:      %d\n', size(buildingBoxes, 1));
fprintf('Free-space polytopes:     %d\n', size(freeBoxes, 1));

% ---------------------------------------------------------------------
% Optional plot.
% ---------------------------------------------------------------------
if doPlot
	func_plotCityDecomposition(buildingBoxes, freeBoxes, domain);
end
end

function buildingBoxes = generateBuildings(nBuildings, domain)
%GENERATEBUILDINGS
%
% Generates non-overlapping rectangular building footprints with different
% dimensions and heights.

buildingBoxes = zeros(0, 6);

cityX = domain.ub(1) - domain.lb(1);
cityY = domain.ub(2) - domain.lb(2);
maxZ  = domain.ub(3);

minWidth  = 4;
maxWidth  = 20;
minDepth  = 4;
maxDepth  = 20;
minHeight = 15;
maxHeight = 70;

spacing = 2.0*0;

maxAttempts = 5000;
attempts = 0;

while size(buildingBoxes, 1) < nBuildings && attempts < maxAttempts
	attempts = attempts + 1;

	w = minWidth  + rand() * (maxWidth  - minWidth);
	d = minDepth  + rand() * (maxDepth  - minDepth);
	h = minHeight + rand() * (maxHeight - minHeight);

	xMin = domain.lb(1) + 5 + rand() * (cityX - w - 10);
	yMin = domain.lb(2) + 5 + rand() * (cityY - d - 10);

	xMax = xMin + w;
	yMax = yMin + d;
	zMin = 0;
	zMax = min(h, maxZ);

	candidate = [xMin xMax yMin yMax zMin zMax];

	if ~overlapsAnyBuildingFootprint(candidate, buildingBoxes, spacing)
		buildingBoxes(end+1, :) = candidate; %#ok<AGROW>
	end
end

if size(buildingBoxes, 1) < nBuildings
	warning(['Only generated %d buildings out of requested %d. ', ...
		'Try reducing nBuildings or increasing domain size.'], ...
		size(buildingBoxes, 1), nBuildings);
end
end

function tf = overlapsAnyBuildingFootprint(candidate, boxes, spacing)
%OVERLAPSANYBUILDINGFOOTPRINT
%
% Checks whether a candidate building footprint overlaps existing buildings.
% A small spacing buffer is included.

tf = false;

if isempty(boxes)
	return;
end

c = candidate;

cInflated = [
	c(1)-spacing, c(2)+spacing, ...
	c(3)-spacing, c(4)+spacing, ...
	c(5),         c(6)
	];

for i = 1:size(boxes, 1)
	b = boxes(i, :);

	overlapX = cInflated(1) < b(2) && cInflated(2) > b(1);
	overlapY = cInflated(3) < b(4) && cInflated(4) > b(3);

	if overlapX && overlapY
		tf = true;
		return;
	end
end
end

function freeBoxes = decomposeFreeSpaceIntoBoxes(domain, obstacleBoxes)
%DECOMPOSEFREESPACEINTOBOXES
%
% Creates a non-overlapping box decomposition of:
%
%   domain \ obstacleBoxes
%
% using an axis-aligned grid induced by all obstacle boundaries, followed by
% greedy box merging.

tol = 1e-10;

xEdges = [domain.lb(1), domain.ub(1)];
yEdges = [domain.lb(2), domain.ub(2)];
zEdges = [domain.lb(3), domain.ub(3)];

for i = 1:size(obstacleBoxes, 1)
	lo = max(obstacleBoxes(i, [1 3 5]), domain.lb);
	hi = min(obstacleBoxes(i, [2 4 6]), domain.ub);

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

for ix = 1:nx
	cx = 0.5 * (xEdges(ix) + xEdges(ix+1));

	for iy = 1:ny
		cy = 0.5 * (yEdges(iy) + yEdges(iy+1));

		for iz = 1:nz
			cz = 0.5 * (zEdges(iz) + zEdges(iz+1));
			p = [cx cy cz];

			if ~pointInsideAnyBox(p, obstacleBoxes, tol)
				freeMask(ix, iy, iz) = true;
			end
		end
	end
end

freeBoxes = func_hyperrectangles_greedyMerge(freeMask, xEdges, yEdges, zEdges);

freeBoxes = func__hyperrectangles_mergeAdjacent(freeBoxes, tol);
end

function inside = pointInsideAnyBox(p, boxes, tol)
%POINTINSIDEANYBOX
%
% Returns true if point p lies inside any obstacle box.

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