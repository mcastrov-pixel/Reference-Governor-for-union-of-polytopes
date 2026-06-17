
function rooms = func_hyperrectangles_genWeakExt(cell_HRepA,cell_HRepb,flagVisu)

rooms.n = length(cell_HRepA);
epsSmall = 1e-13;
if nargin<3
	flagVisu = 0;
end
%generate closed room polyhedrons
rooms.closedRoom.HRep.A = cell_HRepA;
rooms.closedRoom.HRep.b = cell_HRepb;
%combineMat is used to quickly generate vertices of hypercube
combineMat = [1 0 1 0 1 0;
	1 0 1 0 0 1;
	1 0 0 1 1 0;
	1 0 0 1 0 1;
	0 1 1 0 1 0;
	0 1 1 0 0 1;
	0 1 0 1 1 0;
	0 1 0 1 0 1;];
boolMat = boolean(combineMat)';
for it = 1:rooms.n
curBounds = diag([-1 1 -1 1 -1 1])*cell_HRepb{it};
	rooms.closedRoom.VRep{it} = [];
	for jt = 1:size(boolMat,2)
		rooms.closedRoom.VRep{it} =[rooms.closedRoom.VRep{it}  curBounds(boolMat(:,jt))] ;
	end
end

%find what room touch each other as well as the corresponding facet
%matrices should be read column wis as in: the index in position (1,2)
%corresponds to the facet of room2 that touches room1.

rooms.graphs.flags = zeros(rooms.n);
rooms.gates = zeros(rooms.n);
rooms.gatesVerify = zeros(rooms.n);
bVals = cell2mat(rooms.closedRoom.HRep.b);
nDim = size(bVals,1)/2;
for it = 1:rooms.n
	curTouch = false(1,rooms.n);
	for kt = 1:nDim
		touchPlusSide = -bVals((kt-1)*2+1,:)==bVals(kt*2,it);
		touchMinusSide = bVals(kt*2,:)== -bVals((kt-1)*2+1,it);
		curTouchCur = touchPlusSide|touchMinusSide;
		for jt = 1:nDim
			if jt~=kt
				curTouchTemp = curTouchCur&(bVals( jt*2,:)<=bVals(jt*2,it))&(bVals( jt*2,:)>-bVals((jt-1)*2+1,it));
				curTouchTemp = [curTouchTemp;
					curTouchCur&(-bVals((jt-1)*2+1,:)<bVals(jt*2,it))&(-bVals((jt-1)*2+1,:)>=-bVals((jt-1)*2+1,it))];
				curTouchTemp = [curTouchTemp;
					curTouchCur&(-bVals((jt-1)*2+1,it)<bVals(jt*2,:))&(-bVals((jt-1)*2+1,it)>=-bVals((jt-1)*2+1,:))];
				curTouchTemp = [curTouchTemp;
					curTouchCur&(bVals(jt*2,it)<=bVals(jt*2,:))&(bVals(jt*2,it)>-bVals((jt-1)*2+1, :))];
				curTouchCur = any(curTouchTemp,1);
			end
		end
		curTouch = curTouch|curTouchCur;
		rooms.gates(:,it) = rooms.gates(:,it) + (kt*2)*(curTouchCur.*touchPlusSide)';
		rooms.gates(:,it) = rooms.gates(:,it) + ((kt-1)*2+1)*(curTouchCur.*touchMinusSide)';
		rooms.gatesVerify(it,:) = rooms.gatesVerify(it,:) + (kt*2)*(curTouchCur.*touchMinusSide);
		rooms.gatesVerify(it,:) = rooms.gatesVerify(it,:) + ((kt-1)*2+1)*(curTouchCur.*touchPlusSide);
	end
	rooms.graphs.flags(:,it) = curTouch';
end

if flagVisu && 0
	curFig = figure;
	blackSheep = find(sum(rooms.graphs.flags,1)==0);
	kt = 2;
	toBeSeen = find(rooms.graphs.flags(:,kt));
for it=1:rooms.n
	popol1 = polytope(rooms.closedRoom.HRep.A{it},rooms.closedRoom.HRep.b{it});
	popol1 = Polyhedron(popol1.vertices');popol1.plot;hold on;
	xlabel('x axis','FontSize',20);ylabel('y axis','FontSize',20);zlabel('z axis','FontSize',20);
	curAx = curFig.Children;
	curObj = curAx.Children(1);
	curObj.FaceColor = 'y';
	curObj.FaceAlpha = .2;
	if any(it == toBeSeen)
		curObj.FaceColor = 'b';
	curObj.FaceAlpha = 1;
	elseif it == kt
		curObj.FaceColor = 'r';
		curObj.FaceAlpha = 1;
	end
	if any(it==blackSheep)
			curObj.FaceColor = 'k';
	curObj.FaceAlpha = 1;
	end
end
keyboard
end


 %% find the corresponding vertices - if need be
%  rooms.verts = cell(rooms.n);
% for it = 1:rooms.n
% 	touches = sum(rooms.graphs.flags(:,it));
% 	touchIdx = find(rooms.graphs.flags(:,it));
% 	for kt = 1:touches
% 		jt = touchIdx(kt);
% 	rooms.verts{jt,it} = find(abs(rooms.closedRoom.HRep.A{it}(rooms.gates(jt,it),:)*rooms.closedRoom.VRep{jt} - rooms.closedRoom.HRep.b{it}(rooms.gates(jt,it)))<=epsSmall);
% 	end
% end



%% generate openings and intersections
rooms.open.HRep.A = cell(rooms.n);
rooms.open.HRep.b = cell(rooms.n);
rooms.intersect.HRep.A = cell(rooms.n);
rooms.intersect.HRep.b = cell(rooms.n);
for it=1:rooms.n;
	for jt = 1:rooms.n
		if rooms.graphs.flags(jt,it)
			%opening is simply by deleting constraint
			rooms.open.HRep.A{jt,it} = rooms.closedRoom.HRep.A{it};
			rooms.open.HRep.b{jt,it} = rooms.closedRoom.HRep.b{it};
			rooms.open.HRep.A{jt,it}(rooms.gates(jt,it),:) = [];
			rooms.open.HRep.b{jt,it}(rooms.gates(jt,it),:) = [];
			% generally intersection is as follow
			% rooms.intersect.A{jt,it} = [rooms.open.HRep.A{jt,it};rooms.closedRoom.HRep.A{jt}];
			% rooms.intersect.b{jt,it} = [rooms.open.HRep.b{jt,it};rooms.closedRoom.HRep.b{jt}];
			%For hyperrectangle we have a simple way
			takeFromNext = rooms.gates(jt,it);
			if mod(rooms.gates(jt,it),2)
				takeFromNext = [takeFromNext;takeFromNext+1];
			else
				takeFromNext = [takeFromNext-1;takeFromNext];
			end
			room_it = rooms.closedRoom.HRep.b{it};
			room_jt = rooms.closedRoom.HRep.b{jt}; % replace indices we do not need to compare
			room_it(takeFromNext) = room_jt(takeFromNext);
			rooms.intersect.HRep.A{jt,it} = rooms.closedRoom.HRep.A{it};
			rooms.intersect.HRep.b{jt,it}= min(room_it,room_jt);
		end
	end
end


%% generate weak Extensions
rooms.weakExtension.HRep.A = cell(rooms.n);
rooms.weakExtension.HRep.b = cell(rooms.n);
for it = 1:rooms.n
	for jt = 1:rooms.n
		if rooms.graphs.flags(jt,it)
			% in general
			% rooms.weakExtension.HRep.A{jt,it} = [rooms.open.HRep.A{jt,it};rooms.open.HRep.A{it,jt}]
			% rooms.weakExtension.HRep.b{jt,it} = [rooms.open.HRep.b{jt,it};rooms.open.HRep.b{it,jt}]
			% for hyperrectangle: 
			rooms.weakExtension.HRep.A{jt,it} = [rooms.closedRoom.HRep.A{it}];
			takeFromNext = rooms.gates(jt,it);
			if mod(rooms.gates(jt,it),2)
				takeFromPrev = takeFromNext+1;
			else
				takeFromPrev =takeFromNext-1;
			end
			room_it = rooms.closedRoom.HRep.b{it};
			room_jt = rooms.closedRoom.HRep.b{jt};
			room_it(takeFromNext) = room_jt(takeFromNext);
			room_jt(takeFromPrev) = room_it(takeFromPrev);
			rooms.weakExtension.HRep.b{jt,it} = min(room_it,room_jt);
		end
	end
end
rooms.graphs.weights = rooms.graphs.flags;
if flagVisu
	curFig = figure;
	kt = 26;
	toBeSeen = find(rooms.graphs.flags(:,kt));
	popol1 = polytope(rooms.closedRoom.HRep.A{kt},rooms.closedRoom.HRep.b{kt});
	popol1 = Polyhedron(popol1.vertices');popol1.plot;hold on;
	xlabel('x axis','FontSize',20);ylabel('y axis','FontSize',20);zlabel('z axis','FontSize',20);
	curAx = curFig.Children;
	curObj = curAx.Children(1);
	curObj.FaceColor = 'r';
	curObj.FaceAlpha = 0;
	curObj.EdgeColor = 'r';
	curObj.LineWidth = 4;
for it=1:length(toBeSeen)
	popol1 = polytope(rooms.closedRoom.HRep.A{toBeSeen(it)},rooms.closedRoom.HRep.b{toBeSeen(it)});
	popol1 = Polyhedron(popol1.vertices');popol1.plot;
	curObj = curAx.Children(1);
	curObj.FaceColor = 'y';
	curObj.FaceAlpha = 0;
	curObj.LineWidth = 3;
	curObj.LineStyle = '-.';
	popol1 = polytope(rooms.weakExtension.HRep.A{toBeSeen(it),kt},rooms.weakExtension.HRep.b{toBeSeen(it),kt});
	popol1 = Polyhedron(popol1.vertices');popol1.plot;
	curObj = curAx.Children(1);
	curObj.FaceColor = 'g';
	curObj.FaceAlpha = .7;
end
% keyboard
end
% generate centers and scaling of unit cube in each dimension

rooms.transformation.center = cell(rooms.n);
rooms.transformation.scaling = cell(rooms.n);
for it = 1:rooms.n
	for jt = 1:rooms.n
		if rooms.graphs.flags(jt,it)
			bounds = (reshape(rooms.weakExtension.HRep.b{jt,it},[],3)');
			rooms.transformation.center{jt,it}  = mean(bounds.*[-1 1],2);
			rooms.transformation.scaling{jt,it}  = mean(bounds,2);
			if any(mean(bounds,2)<=0)
				warning('negative length cube???')
			end
		elseif it==jt
			bounds = (reshape(rooms.closedRoom.HRep.b{it},[],3)');
			rooms.transformation.center{jt,it}  = mean(bounds.*[-1 1],2);
			rooms.transformation.scaling{jt,it}  = mean(bounds,2);
		end
	end
end
