
function boxes = func__hyperrectangles_mergeAdjacent(boxes, tol)
%MERGEADJACENTBOXES
%
% Repeatedly merges pairs of boxes that share a full face and whose union
% is again a box.

    if isempty(boxes)
        return;
    end

    changed = true;

    while changed
        changed = false;
        n = size(boxes, 1);

        for i = 1:n
            for j = i+1:n
                [canMerge, mergedBox] = tryMergeTwoBoxes(boxes(i, :), boxes(j, :), tol);

                if canMerge
                    boxes(i, :) = mergedBox;
                    boxes(j, :) = [];
                    changed = true;
                    break;
                end
            end

            if changed
                break;
            end
        end
    end
end


function [canMerge, mergedBox] = tryMergeTwoBoxes(a, b, tol)
%TRYMERGETWOBOXES
%
% Checks whether two boxes can be exactly merged into a single box.

    canMerge = false;
    mergedBox = a;

    axisCols = [
        1 2;
        3 4;
        5 6
    ];

    for ax = 1:3
        cols = axisCols(ax, :);
        otherCols = setdiff(1:6, cols);

        % Other coordinate intervals must match.
        if all(abs(a(otherCols) - b(otherCols)) < tol)

            % The intervals along this axis must touch.
            aMin = a(cols(1));
            aMax = a(cols(2));
            bMin = b(cols(1));
            bMax = b(cols(2));

            touchAB = abs(aMax - bMin) < tol;
            touchBA = abs(bMax - aMin) < tol;

            if touchAB || touchBA
                mergedBox = a;
                mergedBox(cols(1)) = min(aMin, bMin);
                mergedBox(cols(2)) = max(aMax, bMax);
                canMerge = true;
                return;
            end
        end
    end
end