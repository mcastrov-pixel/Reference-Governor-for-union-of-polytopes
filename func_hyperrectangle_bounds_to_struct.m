function roomsPre = func_hyperrectangle_bounds_to_struct(obstacleBoxes, freeBoxes, domain)
%go from matrix delimiting the free/used space to a structure containing
%the different elements in the connected sequence, obstacle information as
%polytopes,...

nObstacles= size(obstacleBoxes,1);
nFreeSpace = size(freeBoxes,1);
nBounds = size(obstacleBoxes,2);

AMat = kron(eye(3),[-1;1]);

roomsPre.buildingBoxes = obstacleBoxes;
roomsPre.freeBoxes = freeBoxes;
roomsPre.domain = domain;

obstacleBoxes = obstacleBoxes.*[-1 1 -1 1 -1 1];
freeBoxes = freeBoxes.*[-1 1 -1 1 -1 1];


%incorporate the building data
roomsPre.obstacles.HRep.A = cell(1,nObstacles);
roomsPre.obstacles.HRep.A(:) = {AMat};
roomsPre.obstacles.HRep.b = mat2cell(reshape(obstacleBoxes',[],1),ones(1,nObstacles)*nBounds,1)';

%limits of the city
roomsPre.obstacles.limits = domain;
roomsPre.obstacles.limits.lb = roomsPre.obstacles.limits.lb';
roomsPre.obstacles.limits.ub = roomsPre.obstacles.limits.ub';

%incorporate the free space data
roomsPre.n = nFreeSpace;

roomsPre.closedRoom.HRep.A = cell(1,nFreeSpace);
roomsPre.closedRoom.HRep.A(:) = {AMat};
roomsPre.closedRoom.HRep.b = mat2cell(reshape(freeBoxes',[],1),ones(1,nFreeSpace)*nBounds,1)';

end