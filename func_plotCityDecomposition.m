function curFig = func_plotCityDecomposition(buildingBoxes, freeBoxes, domain,curFig_num)
%PLOTCITYDECOMPOSITION
%
% Plots free-space boxes and buildings.
if nargin<4|| isempty(curFig_num)
    curFig = figure('Color', 'w');
else
	curFig = figure(curFig_num);
	curFig.Color='w';
end
    hold on;
    grid off;
    axis equal;
    axis vis3d;

    xlabel('x');
    ylabel('y');
    zlabel('z');

    xlim([domain.lb(1), domain.ub(1)]);
    ylim([domain.lb(2), domain.ub(2)]);
    zlim([domain.lb(3), domain.ub(3)]);

	xlabel('x [m]','FontSize',18,'interpreter','latex');xticks([domain.lb(1) domain.ub(1)/2 domain.ub(1)] );
	ylabel('y [m]','FontSize',18,'interpreter','latex');yticks([domain.lb(2) domain.ub(2)/2 domain.ub(2)] );
	zlabel('z [m]','FontSize',18,'interpreter','latex');zticks([domain.lb(3) domain.ub(3)/2 domain.ub(3)] );

   
    view(38, 25);

    % Plot free space first.
    hFree = [];

    for i = 1:size(freeBoxes, 1)
        h = drawBoxPatch(freeBoxes(i, :), [0.1 0.45 1.0], 0.045, [0.2,0.2,0.2]);
        if isempty(hFree)
            hFree = h;
        end
    end

    % Plot ground plane.
    drawGroundPlane(domain);


    % Plot buildings.
    hBuilding = [];

    for i = 1:size(buildingBoxes, 1)
        h = drawBuilding(buildingBoxes(i, :), i);
        if isempty(hBuilding)
            hBuilding = h;
        end
    end


    camlight headlight;
    lighting gouraud;
end


function drawGroundPlane(domain)
%DRAWGROUNDPLANE
%
% Draws a light gray ground plane.

    x0 = domain.lb(1);
    x1 = domain.ub(1);
    y0 = domain.lb(2);
    y1 = domain.ub(2);
    z0 = domain.lb(3);

    V = [
        x0 y0 z0;
        x1 y0 z0;
        x1 y1 z0;
        x0 y1 z0
    ];

    F = [1 2 3 4];

    patch( ...
        'Vertices', V, ...
        'Faces', F, ...
        'FaceColor', [1 1 1]*0.45, ...
        'FaceAlpha', 0.35, ...
        'EdgeColor', 'none');
end


function h = drawBuilding(box, idx)
%DRAWBUILDING
%
% Draws a building-like rectangular prism with windows and a roof.
% The geometry used for decomposition is still just the input box.

    baseColors = [
        0.45 0.45 0.50;
        0.55 0.50 0.45;
        0.40 0.48 0.55;
        0.50 0.42 0.42;
        0.42 0.50 0.44
    ];

    colorIdx = mod(idx-1, size(baseColors, 1)) + 1;
    buildingColor = baseColors(colorIdx, :);

    h = drawBoxPatch(box, buildingColor, 0.92, [0.15 0.15 0.15]);

    % Slightly darker roof cap.
    drawRoofCap(box);

    % Add windows to two visible side families.
    drawWindowsOnBuilding(box);
end


function drawRoofCap(box)
%DRAWROOFCAP
%
% Draws a simple roof cap.

    xmin = box(1);
    xmax = box(2);
    ymin = box(3);
    ymax = box(4);
    zmax = box(6);

    roofThickness = 0.4;
    overhang = 0.4*0;

    roofBox = [
        xmin-overhang, xmax+overhang, ...
        ymin-overhang, ymax+overhang, ...
        zmax-roofThickness,          zmax
    ];

    drawBoxPatch(roofBox, [0.18 0.18 0.20], 1.0, [0.05 0.05 0.05]);
end


function drawWindowsOnBuilding(box)
%DRAWWINDOWSONBUILDING
%
% Adds simple blue window rectangles to the building sides.
% These are visual only and do not affect the polytope representation.

    xmin = box(1);
    xmax = box(2);
    ymin = box(3);
    ymax = box(4);
    zmin = box(5);
    zmax = box(6);

    width = xmax - xmin;
    depth = ymax - ymin;
    height = zmax - zmin;

    floorSpacing = 4.0;
    windowW = 1.0;
    windowH = 1.5;

    nFloors = max(1, floor(height / floorSpacing) - 1);
    nX = max(1, floor(width / 4));
    nY = max(1, floor(depth / 4));

    % Cap the number of windows so the plot stays responsive.
    nFloors = min(nFloors, 12);
    nX = min(nX, 8);
    nY = min(nY, 8);

    windowColor = [0.3 0.75 1.0];
    alpha = 0.75;

    epsOut = 0.02;

    % Front face, y = ymin.
    for f = 1:nFloors
        zc = zmin + f * floorSpacing;

        for k = 1:nX
            xc = xmin + k * width / (nX + 1);

            drawWindowRectangle( ...
                [xc-windowW/2, ymin-epsOut, zc-windowH/2], ...
                [xc+windowW/2, ymin-epsOut, zc-windowH/2], ...
                [xc+windowW/2, ymin-epsOut, zc+windowH/2], ...
                [xc-windowW/2, ymin-epsOut, zc+windowH/2], ...
                windowColor, alpha);
        end
    end

    % Right face, x = xmax.
    for f = 1:nFloors
        zc = zmin + f * floorSpacing;

        for k = 1:nY
            yc = ymin + k * depth / (nY + 1);

            drawWindowRectangle( ...
                [xmax+epsOut, yc-windowW/2, zc-windowH/2], ...
                [xmax+epsOut, yc+windowW/2, zc-windowH/2], ...
                [xmax+epsOut, yc+windowW/2, zc+windowH/2], ...
                [xmax+epsOut, yc-windowW/2, zc+windowH/2], ...
                windowColor, alpha);
        end
    end
end


function drawWindowRectangle(p1, p2, p3, p4, faceColor, faceAlpha)
%DRAWWINDOWRECTANGLE
%
% Draws a single rectangular window patch.

    V = [p1; p2; p3; p4];
    F = [1 2 3 4];

    patch( ...
        'Vertices', V, ...
        'Faces', F, ...
        'FaceColor', faceColor, ...
        'FaceAlpha', faceAlpha, ...
        'EdgeColor', 'none');
end


function h = drawBoxPatch(box, faceColor, faceAlpha, edgeColor)
%DRAWBOXPATCH
%
% Draws an axis-aligned rectangular box.

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
% Draws a wireframe bounding box.

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

