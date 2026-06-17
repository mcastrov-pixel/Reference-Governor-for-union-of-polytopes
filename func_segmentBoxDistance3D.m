function [dist, tClosest, closestSeg, closestBox] = func_segmentBoxDistance3D(p0, p1, box)
%SEGMENTBOXDISTANCE3D Distance between a 3D axis-aligned box and a segment.
%
%   dist = segmentBoxDistance3D(p0, p1, box)
%
%   [dist, tClosest, closestSeg, closestBox] = segmentBoxDistance3D(...)
%
%   Inputs:
%       p0  - 1x3 or 3x1 vector, first endpoint of segment
%       p1  - 1x3 or 3x1 vector, second endpoint of segment
%       box - 1x6 vector:
%             [xmin xmax ymin ymax zmin zmax]
%
%   Outputs:
%       dist       - Euclidean distance between the segment and the box
%       tClosest   - segment parameter in [0,1] giving closest point
%       closestSeg - closest point on segment
%       closestBox - closest point on box
%
%   The segment is parameterized as:
%       p(t) = p0 + t * (p1 - p0),   0 <= t <= 1
%
%   Distance is zero if the segment intersects the box.

    p0 = p0(:);
    p1 = p1(:);
    box = box(:).';

    if numel(p0) ~= 3 || numel(p1) ~= 3
        error('p0 and p1 must be 3D vectors.');
    end

    if numel(box) ~= 6
        error('box must be a 1x6 vector: [xmin xmax ymin ymax zmin zmax].');
    end

    mn = [box(1); box(3); box(5)];
    mx = [box(2); box(4); box(6)];

    if any(mn > mx)
        error('Invalid box: each min value must be <= corresponding max value.');
    end

    d = p1 - p0;

    % First check whether the segment intersects the box.
    [intersects, tHit] = segmentIntersectsAABB(p0, p1, mn, mx);

    if intersects
        dist = 0;
        tClosest = tHit;
        closestSeg = p0 + tClosest * d;
        closestBox = closestSeg;
        return;
    end

    % The squared distance from p(t) to the box is a convex piecewise
    % quadratic function of t. Breakpoints occur when p(t) crosses one of
    % the six box planes.

    breakpoints = [0; 1];

    for i = 1:3
        if d(i) ~= 0
            tMinPlane = (mn(i) - p0(i)) / d(i);
            tMaxPlane = (mx(i) - p0(i)) / d(i);

            if tMinPlane > 0 && tMinPlane < 1
                breakpoints(end+1,1) = tMinPlane;
            end

            if tMaxPlane > 0 && tMaxPlane < 1
                breakpoints(end+1,1) = tMaxPlane;
            end
        end
    end

    breakpoints = unique(sort(breakpoints));

    bestD2 = inf;
    tClosest = 0;

    % Search each interval between breakpoints.
    for k = 1:numel(breakpoints)-1
        a = breakpoints(k);
        b = breakpoints(k+1);

        candidateT = [a; b];

        if b > a
            tm = 0.5 * (a + b);

            % On this interval, determine which side of the box each
            % coordinate lies on. This determines the quadratic form.
            A = 0;
            B = 0;

            for i = 1:3
                xmid = p0(i) + tm * d(i);

                if xmid < mn(i)
                    target = mn(i);
                elseif xmid > mx(i)
                    target = mx(i);
                else
                    continue;
                end

                % Contribution:
                %   p_i(t) - target = p0_i + d_i t - target
                %
                % Squared term:
                %   d_i^2 t^2 + 2 d_i(p0_i - target)t + ...
                A = A + d(i)^2;
                B = B + d(i) * (p0(i) - target);
            end

            % Minimize A t^2 + 2 B t + C on this interval.
            if A > 0
                tStar = -B / A;

                if tStar > a && tStar < b
                    candidateT(end+1,1) = tStar;
                end
            end
        end

        % Evaluate all candidate parameters.
        for j = 1:numel(candidateT)
            t = candidateT(j);
            q = p0 + t * d;
            c = clampPointToBox(q, mn, mx);

            d2 = sum((q - c).^2);

            if d2 < bestD2
                bestD2 = d2;
                tClosest = t;
            end
        end
    end

    closestSeg = p0 + tClosest * d;
    closestBox = clampPointToBox(closestSeg, mn, mx);
    dist = sqrt(max(bestD2, 0));
end


function c = clampPointToBox(p, mn, mx)
%CLAMPPOINTTOBOX Closest point in an axis-aligned box to point p.
    c = min(max(p, mn), mx);
end


function [intersects, tHit] = segmentIntersectsAABB(p0, p1, mn, mx)
%SEGMENTINTERSECTSAABB Slab test for segment-box intersection.

    d = p1 - p0;

    tEnter = 0;
    tExit = 1;

    tol = 1e-14;

    for i = 1:3
        if abs(d(i)) < tol
            % Segment is parallel to this slab. It must already lie inside.
            if p0(i) < mn(i) || p0(i) > mx(i)
                intersects = false;
                tHit = NaN;
                return;
            end
        else
            t1 = (mn(i) - p0(i)) / d(i);
            t2 = (mx(i) - p0(i)) / d(i);

            tNear = min(t1, t2);
            tFar  = max(t1, t2);

            tEnter = max(tEnter, tNear);
            tExit  = min(tExit, tFar);

            if tEnter > tExit
                intersects = false;
                tHit = NaN;
                return;
            end
        end
    end

    intersects = true;
    tHit = tEnter;
end