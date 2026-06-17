function boxes = func_hyperrectangles_greedyMerge(freeMask, xEdges, yEdges, zEdges)
%GREEDYMERGEGRIDCELLS
%
% Converts a 3D logical free-cell mask into larger non-overlapping boxes.

    [nx, ny, nz] = size(freeMask);
    boxes = zeros(0, 6);

    while any(freeMask(:))
        idx = find(freeMask, 1, 'first');
        [i1, j1, k1] = ind2sub(size(freeMask), idx);

        i2 = i1;
        j2 = j1;
        k2 = k1;

        % Expand in the positive x/y/z directions while possible.
        keepExpanding = true;

        while keepExpanding
            keepExpanding = false;

            gains = [-Inf, -Inf, -Inf];

            % Try expanding in x.
            if i2 < nx
                slab = freeMask(i2+1, j1:j2, k1:k2);
                if all(slab(:))
                    dx = xEdges(i2+2) - xEdges(i2+1);
                    dy = yEdges(j2+1) - yEdges(j1);
                    dz = zEdges(k2+1) - zEdges(k1);
                    gains(1) = dx * dy * dz;
                end
            end

            % Try expanding in y.
            if j2 < ny
                slab = freeMask(i1:i2, j2+1, k1:k2);
                if all(slab(:))
                    dx = xEdges(i2+1) - xEdges(i1);
                    dy = yEdges(j2+2) - yEdges(j2+1);
                    dz = zEdges(k2+1) - zEdges(k1);
                    gains(2) = dx * dy * dz;
                end
            end

            % Try expanding in z.
            if k2 < nz
                slab = freeMask(i1:i2, j1:j2, k2+1);
                if all(slab(:))
                    dx = xEdges(i2+1) - xEdges(i1);
                    dy = yEdges(j2+1) - yEdges(j1);
                    dz = zEdges(k2+2) - zEdges(k2+1);
                    gains(3) = dx * dy * dz;
                end
            end

            [bestGain, bestAxis] = max(gains);

            if isfinite(bestGain)
                keepExpanding = true;

                if bestAxis == 1
                    i2 = i2 + 1;
                elseif bestAxis == 2
                    j2 = j2 + 1;
                elseif bestAxis == 3
                    k2 = k2 + 1;
                end
            end
        end

        % Add box corresponding to merged block.
        newBox = [
            xEdges(i1), xEdges(i2+1), ...
            yEdges(j1), yEdges(j2+1), ...
            zEdges(k1), zEdges(k2+1)
        ];

        boxes = [boxes; newBox]; %#ok<AGROW>

        % Remove these cells from the mask.
        freeMask(i1:i2, j1:j2, k1:k2) = false;
    end
end