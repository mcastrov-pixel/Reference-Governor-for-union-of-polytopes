function cur_fig = func_plotISSDecomposition(obsBoxes, freeBoxes, domain,figNum)
%func_plotISSDecomposition
%
% Plots obstacle boxes and free-space boxes with transparency.

if nargin<4 || isempty(figNum)
	cur_fig = figure('Color', 'w');
else
	cur_fig = figure(figNum);
	cur_fig.Color = 'W';
end
    hold on;
    grid off;
    axis equal;
    axis vis3d;
    view(35, 25);

    xlabel('x');
    ylabel('y');
    zlabel('z');

    xlim([domain.lb(1), domain.ub(1)]);
    ylim([domain.lb(2), domain.ub(2)]);
    zlim([domain.lb(3), domain.ub(3)]);

    
    % Plot free space first so the spacecraft remains visible.
    hFree = [];

    for i = 1:size(freeBoxes, 1)
        h = drawBoxPatch(freeBoxes(i, :), [0.15 0.45 1.0], 0.1, [0.9,0.9,0.9]);
        if isempty(hFree)
            hFree = h;
        end
    end

    % Plot obstacle boxes.
    hObs = [];

    for i = 1:size(obsBoxes, 1)
        h = drawBoxPatch(obsBoxes(i, :), [0.95,0.95,0.95], 0.75*0+1, [0.25 0.05 0.02]);
        if isempty(hObs)
            hObs = h;
        end
    end

	xticks([domain.lb(1) 0 domain.ub(1)] );xlabel('cross-track','Interpreter','latex','FontSize',18)
	yticks([domain.lb(2) 0 domain.ub(2)] );ylabel('in-track','Interpreter','latex','FontSize',18)
	zticks([domain.lb(3) 0 domain.ub(3)] );zlabel('radial','Interpreter','latex','FontSize',18)
    camlight headlight;
    lighting gouraud;
end


function h = drawBoxPatch(box, faceColor, faceAlpha, edgeColor)
%DRAWBOXPATCH
%
% Draws a rectangular box as a patch object.

    xmin = box(1);
    xmax = box(2);
    ymin = box(3);
    ymax = box(4);
    zmin = box(5);
    zmax = box(6);

    V = [
        xmin ymin zmin;
        xmax ymin zmin;
        xmax ymax zmin;
        xmin ymax zmin;
        xmin ymin zmax;
        xmax ymin zmax;
        xmax ymax zmax;
        xmin ymax zmax
    ];

    F = [
        1 2 3 4;
        5 8 7 6;
        1 5 6 2;
        2 6 7 3;
        3 7 8 4;
        4 8 5 1
    ];

    h = patch( ...
        'Vertices', V, ...
        'Faces', F, ...
        'FaceColor', faceColor, ...
        'FaceAlpha', faceAlpha, ...
        'EdgeColor', edgeColor);
end


function drawBoxWireframe(box, color, lineWidth)
%DRAWBOXWIREFRAME
%
% Draws a wireframe box.

    xmin = box(1);
    xmax = box(2);
    ymin = box(3);
    ymax = box(4);
    zmin = box(5);
    zmax = box(6);

    V = [
        xmin ymin zmin;
        xmax ymin zmin;
        xmax ymax zmin;
        xmin ymax zmin;
        xmin ymin zmax;
        xmax ymin zmax;
        xmax ymax zmax;
        xmin ymax zmax
    ];

    E = [
        1 2;
        2 3;
        3 4;
        4 1;
        5 6;
        6 7;
        7 8;
        8 5;
        1 5;
        2 6;
        3 7;
        4 8
    ];

    for i = 1:size(E, 1)
        plot3(V(E(i,:),1), V(E(i,:),2), V(E(i,:),3), ...
              'Color', color, ...
              'LineWidth', lineWidth);
    end
end

