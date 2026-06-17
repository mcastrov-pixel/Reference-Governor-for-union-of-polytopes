
%Main file for on-orbit proximity operation example
%Optimization problem are solved with quadprog using Yalmip %[https://yalmip.github.io/download/]
% uses MPT3 for polytope ploting: https://www.mpt3.org/pmwiki.php/Main/Installation
% Cora toolbox for some light operations on polytopes: https://tumcps.github.io/CORA/


clearvars

folderDirectory = fileparts(mfilename('fullpath'));
cd(folderDirectory);
idNum = 1;


%% 0 -- Model & parameters definition


flg_loadMoas = 0; %compute MOAS or load it
epsVal = 1e-16; % eps value for inner approximation of MOAS (increase if computation takes a lot of time)
sys = func_CWH_dynamics();
name_base = strcat(sys.id,num2str(idNum),'_');



%%  -------- constraints
% we construct the MOAS for different constraints separately and then
% use the intersection of MOAS.


% 1) box constraints: assumes that only the first three outputs are given
box_Ay = kron(eye(3),[-1;1]);
box_by = ones(6,1);
box_C = sys.Cdt(1:3,:);
box_D = sys.Ddt(1:3,:);

% 2) other constraints
% 2.1)velocity constraints
xdotLim = 1.5;
xdot_Ay = kron( eye(sys.nx/2),[-1;1]); % if only velocity output
xdot_Ay = [zeros(size(xdot_Ay,1),3), xdot_Ay,zeros(size(xdot_Ay,1),sys.nv)]; % increase to include the other outputs of the sys.tem ( pos and input)
xdot_by = xdotLim*ones(sys.nx/2*2,1);

%2.2)input constraints
uLim = .1;
u_Ay = kron(eye(sys.nu),[-1;1]);
u_by = ones(sys.nu*2,1)*uLim;
input_Ay = u_Ay*[sys.x2u sys.v2u];
V_input = Polyhedron(u_Ay,u_by).V;

%2.3)assemble it together
other_Ay = [input_Ay;xdot_Ay];
other_by = [u_by;xdot_by];


% generate spacecraft obstacle / free space

[issBoxes, freeBoxes,domain] = func_CWH_space_decomposition(1);
freeBoxes = func_cleanFreeBoxes(freeBoxes,[3 175],1.5);
roomsPre = func_hyperrectangle_bounds_to_struct(issBoxes, freeBoxes, domain);
rooms = func_hyperrectangles_genWeakExt(roomsPre.closedRoom.HRep.A,roomsPre.closedRoom.HRep.b,0);


%% -------- Compute MOAS


if ~flg_loadMoas
	% generate the base MOASES for  the unit cube
	for it = 1:3
		[ MOAS_box.A{it}, MOAS_box.b{it}] = func_genMOAS(sys.Adt,sys.Bdt,box_C(it,:),box_D(it,:),box_Ay((it-1)*2+1:it*2,it),box_by((it-1)*2+1:it*2),epsVal);
	end
	% generate the MOAS for the rest of the constraints
	[ MOAS_other.A, MOAS_other.b] = func_genMOAS(sys.Adt,sys.Bdt,sys.Cdt,sys.Ddt,other_Ay,other_by,epsVal);

	%generate the MOAS for each room (intersection + scaling + recentering)
	rooms = func_genMOASRooms(rooms,MOAS_box,MOAS_other,sys.v2x);

	% Divide the weak extension into each restriction
	rooms = func_genMOASRestrictions(rooms,sys.nv,sys.nx,sys.v2x);

	%save the moas
	save(strcat(name_base,'MOAS'),'rooms','roomsPre',"issBoxes","freeBoxes")
else
	load(strcat(name_base,'MOAS'));
end

%% Simulation

% generate the desired setpoint and find associated element in sequence
rDes = [0,30.2,-8]';
%find box corresponding to desired set point
isInBox = [freeBoxes(:,1:2:6)-rDes', rDes'-freeBoxes(:,2:2:6)];
isInBox = all(isInBox<=0,2);
list_endElem = find(isInBox); list_endElem = list_endElem(1);


% generate the initial conditions and find associated element in sequence
% ground locations
list_startElem = ones(1,0);
p0_list = [35 -35 8 ;-32 6 5;10000 1000 1000; 0 0 0];
lemons =[];

for it =1:size(p0_list,1)
	isInBox = [freeBoxes(:,1:2:6)-p0_list(it,:), p0_list(it,:)-freeBoxes(:,2:2:6)];
	isInBox = all(isInBox<=0,2);
	tempBox = find(isInBox);
	if isempty(tempBox) 
		lemons = [lemons it];
		sprintf("position number %d is not inside the connected sequence. ignored.",it);
	else
		list_startElem(end+1) = tempBox(1);
	end

end
p0_list(lemons,:) = [];
x0_list = sys.v2x*p0_list';
list_startEndPairs = [kron(list_startElem, ones(size(list_endElem))); kron(ones(size(list_startElem)),list_endElem)];


%set up yalmip variables and options

YAL_v = sdpvar(sys.nv,1);
YAL_opts = sdpsettings('solver','quadprog','verbose',0,'cachesolvers',1);


% loop over different initial conditions
for zt  =1:size(list_startEndPairs,2)
	sprintf("simulation %d out of  %d",[zt,size(list_startEndPairs,2)])

	start_elem = list_startEndPairs(1,zt);
	end_elem = list_startEndPairs(2,zt);
	x0 = x0_list(:,zt);

	% ---------------------- generate the weights, graph and path  ----------------------
	rooms.graphs.weightsAlt = func_genDistBasedWeights(rooms,x0(1:3),sys.v2x(1:3,:)*rDes);
	rooms.graphs.graphObject = digraph(rooms.graphs.weightsAlt); % replace with "graph(rooms.graphs.flags);" to get unweighted graph
	path_rooms = shortestpath(rooms.graphs.graphObject,start_elem,end_elem);
	path_length = length(path_rooms);
	path_refs = func_genPathSetPoints(rooms,path_rooms,rDes,sys.v2x,YAL_v,YAL_opts);
	flagInfeas = [];



	% simulation length depends on # of rooms to go through
	NSim = max(round(50*path_length),400);

	% initiate some buffers
	xHist = zeros(sys.nx,NSim + 1);
	vHist = zeros(sys.nv,NSim+1);
	rHist = zeros(sys.nv,NSim+1);
	xHist(:,1) = x0;
	vHist(:,1) = sys.x2v*x0;

	% initiate where we are on the path
	cur_IDX = 1;
	r_IDX = cur_IDX;
	next_IDX = min(cur_IDX + 1,path_length);
	next_room = path_rooms(next_IDX);
	cur_room = path_rooms(cur_IDX);
	prev_room = cur_room;
	tollNew = 1e-7;

	%simulation
	for it = 1:NSim

		% logic to decide in what MOAS we are located
		if cur_IDX<path_length
			%check if we are ready to go from current room to next room
			% --> verify if current state + prev ref in weak extension OR current state + assoc steady state is in weak extension
			test1 = all(rooms.MOAS.A{cur_room,next_room}*[vHist(:,it);xHist(:,it)]-rooms.MOAS.b{cur_room,next_room}<=0);
			test2 = all(rooms.MOAS.A{cur_room,next_room}*[sys.x2v*xHist(:,it);xHist(:,it)]-rooms.MOAS.b{cur_room,next_room}<=0);
			if test1||test2
				cur_IDX = cur_IDX+1;
				sprintf("Progressed to next weak extension at time %d",it)
				prev_room = cur_room;
				next_IDX = min(cur_IDX + 1,path_length);
				next_room = path_rooms(next_IDX);
				cur_room = path_rooms(cur_IDX);
			end
		end

		% hypothesis we are in the MOAS of the current room ("current"  after the update test
		% --> IF NEITHER: current state + prev ref in current room; NOR current state + assoc steady state in current room=
		% we must be in the weak extension and need to transition to the other side of the weak extension (update rIntermediary)

		cur_MOAS.A = rooms.MOAS.A{cur_room,cur_room};
		cur_MOAS.b = rooms.MOAS.b{cur_room,cur_room};
		rHist(:,it) = path_refs(:,2*(cur_IDX-1)+1);
		test1 = ~all(cur_MOAS.A*[vHist(:,it);xHist(:,it)]-cur_MOAS.b<=tollNew);
		test2 = ~all(cur_MOAS.A*[sys.x2v*xHist(:,it);xHist(:,it)]-cur_MOAS.b<=tollNew);
		if test1*test2
			cur_MOAS.A = rooms.MOAS.A{prev_room,cur_room};
			cur_MOAS.b = rooms.MOAS.b{prev_room,cur_room};
			rHist(:,it) = path_refs(:,2*(cur_IDX-1));
		end


		% compute current reference command using command governor
		[vHist(:,it+1),flagInfeas(end+1)] = func_commandGovernor(YAL_v,YAL_opts,cur_MOAS.A,cur_MOAS.b,sys.nv,rHist(:,it),xHist(:,it),vHist(:,it));

		% propagate the dynamics one step forward
		xHist(:,it+1) = sys.Adt*xHist(:,it) + sys.Bdt*vHist(:,it+1);
	end

	%save buffers
	DATA.xHist{zt} = xHist(:,1:end-1);
	DATA.vHist{zt} = vHist(:,2:end);
	DATA.rHist{zt} = rHist(:,1:end-1);
	DATA.pathRefs{zt} = path_refs;
	DATA.pathRooms{zt} = path_rooms;
	DATA.flagInfeas{zt} = flagInfeas;
	DATA.flagRoomConv{zt} = cur_IDX==length(path_rooms);
	DATA.distToDes{zt} = norm(xHist(:,end) - sys.v2x*rDes);
end
save(strcat(name_base,'DATA'),'DATA','sys');


%% trajectory plot

	close (figure(12))
	cur_fig = func_plotISSDecomposition(issBoxes, [], domain,12);

	for it = 1:length(DATA.xHist)

	hold on;
	plot3(DATA.xHist{it}(1,:),DATA.xHist{it}(2,:),DATA.xHist{it}(3,:),'c-','LineWidth',1.5);hold on
	xStart = DATA.xHist{it}(:,1);
	plot3(xStart(1),xStart(2),xStart(3),'g>','LineWidth',6)
	hold off
	func_showRooms(rooms,DATA.pathRooms{it},cur_fig,'none',.41)
	end
	hold on;
	xEnd = sys.v2x*rDes; plot3(xEnd(1,1),xEnd(2,1),xEnd(3),'m<','LineWidth',6)
	hold off;


%% plot velocities


figure(67)
plot((1:length(xHist))*sys.dt,max(abs(xHist(4:6,:)),[],1));hold on;
cellfun(@(y) plot((1:length(y))*sys.dt,max(abs(y(4:6,:)),[],1)),DATA.xHist(1:end),'uniformOutput',false)
hold off;
ylabel('max. velocity component [m/s]','FontSize',18,'interpreter','latex')
xlabel('time [s]','FontSize',18,'interpreter','latex')

%% plot input


	cur_fig = figure(78);
	%plot the input constraint set
	temp_pol = Polyhedron('V',V_input);
	temp_pol.plot();hold on;
	cur_ax = findobj(cur_fig,'Type','axes');
	curObj = cur_ax.Children(1);
	curObj.FaceColor = 'r';
	curObj.FaceAlpha = .2;
	hold on;
	maxCstrVal = {};
	for it = 1:length(DATA.xHist)
		%compute the torque input
		uVals = sys.x2u*DATA.xHist{it}(:,1:end)+sys.v2u*DATA.vHist{it};
		plot3( uVals(1,:),uVals(2,:),uVals(3,:),'LineWidth',1.6);
		maxCstrVal{it} = max(u_Ay*uVals -u_by,[],1);
	end
	xlabel('$u_1$ (x)','FontSize',18,'interpreter','latex')
	ylabel('$u_2$ (y)','FontSize',18,'interpreter','latex')
	zlabel('$u_3$ (z)','FontSize',18,'interpreter','latex')
	hold off

%% check if the trajectory collided with any of the obstacles

bMat = cell2mat(roomsPre.obstacles.HRep.b);
AMat = roomsPre.obstacles.HRep.A{1}*[eye(3) zeros(3)];
DATA.cstrViol = cell(1);
for it = 1:length(DATA.xHist)
	xHistCur = DATA.xHist{it};
	DATA.cstrViol{it} = [];
	if ~isempty(xHistCur)
		for jt = 1:length(xHistCur)
			tempVal = all(AMat*xHistCur(:,jt)<=bMat-1e-6,1); % we allow a tolerance of 1e-6
			DATA.cstrViol{it} = [DATA.cstrViol{it} any(tempVal)];
			if sum(tempVal)>1
				keyboard
			end

		end
	end



end
cstrViolationAtAnyTime = cellfun(@(y) any(y),DATA.cstrViol);
sprintf('number of simulations which collided with an obstacle: %d',sum(any(cstrViolationAtAnyTime,2)))

%% functions

% functions related to the MOAS
function rooms = func_genMOASRooms(rooms,MOAS_box,MOAS_other,v2x)
rooms.MOAS.A = cell(rooms.n);
rooms.MOAS.b = cell(rooms.n);
for it = 1:rooms.n
	for jt = 1:rooms.n
		ACur = MOAS_other.A;bCur = MOAS_other.b;
		if rooms.graphs.flags(jt,it) || it ==jt
			for kt = 1:3
				ABoxCur = MOAS_box.A{kt};
				xv0 = [rooms.transformation.center{jt,it};v2x*rooms.transformation.center{jt,it}];
				bBoxCur = MOAS_box.b{kt}*rooms.transformation.scaling{jt,it}(kt)+ABoxCur*xv0;
				ACur = [ACur;ABoxCur]; bCur = [bCur;bBoxCur];
			end
			rooms.MOAS.A{jt,it} = ACur; rooms.MOAS.b{jt,it} = bCur;
		end
	end
end
end


function rooms = func_genMOASRestrictions(rooms,nv,nx,v2x)
% computes the intersection between the MOAS of a weak extension and the
% next polytope
rooms.MOAS.restrict.A = cell(rooms.n,rooms.n);
rooms.MOAS.restrict.b = cell(rooms.n,rooms.n);

for it = 1:rooms.n
	for jt = 1:rooms.n
		if rooms.graphs.flags(jt,it)
			AyNext = rooms.closedRoom.HRep.A{jt}*[v2x(1:3,1:3) zeros(nv,nx)];
			byNext = rooms.closedRoom.HRep.b{jt};
			rooms.MOAS.restrict.A{jt,it} = [rooms.MOAS.A{jt,it}; AyNext];
			rooms.MOAS.restrict.b{jt,it} = [rooms.MOAS.b{jt,it}; byNext];
		end
	end
end
end


% functions related to the command governor
function [vOut,flagInfeas] = func_commandGovernor(YAL_v,YAL_opts,MOAS_A,MOAS_b,nv,rDes,xCur,vPrev)
if ~isempty(xCur)
	bCur = MOAS_b - MOAS_A(:,nv+1:end)*xCur;
else
	bCur = MOAS_b;
end
YAL_cstr = MOAS_A(:,1:nv)*YAL_v<=bCur;
YAL_obj = (YAL_v - rDes)'*(YAL_v - rDes);
assign(YAL_v,vPrev);
%solve the optimization problem
YAL_diagnosis = optimize(YAL_cstr,YAL_obj,YAL_opts);
vOut = value(YAL_v);
flagInfeas = 0;
% check if we had issues during the solve
if YAL_diagnosis.problem
	%allow for some tolerance (with respect to previous solution
	% before declaring a problem has occured
	if max(MOAS_A(:,1:nv)*vPrev - bCur)>1e-3 || YAL_diagnosis.problem ~=1
		YAL_diagnosis.info;
	end
	vOut = vPrev;
	flagInfeas = 1;
end
end

%functions related to the path
function weights = func_genDistBasedWeights(rooms,pos_0,pos_1)
% computes a weighted graph based on the distance between hyperrectangle and the segment connecting pos 0 and pos 1.
weights = zeros(rooms.n);
for it = 1:rooms.n
	for jt = 1:rooms.n
		if rooms.graphs.flags(jt,it)
				weights(jt,it) = 0.1+func_segmentBoxDistance3D(pos_0, pos_1, rooms.weakExtension.HRep.b{jt,it}'.*[-1 1 -1 1 -1 1]);
		end
	end
end
end


function rPath = func_genPathSetPoints(rooms,path,r_des,v2x,YAL_v,YAL_opts)
eps = 1e-5;
nv = length(r_des);
nPath = length(path);
rPath = zeros(nv,2*(nPath-1)+1);
rPath(:,end) = r_des;
HMat = [eye(nv);v2x];

for it = flip(1:nPath-1)

	MOAS_A = rooms.MOAS.restrict.A{path(it+1),path(it)}*HMat;
	MOAS_b = rooms.MOAS.restrict.b{path(it+1),path(it)}-1e-4; %tightining;
	xCur = [];
	rCur = rPath(:,2*it+1);
	vCur = rooms.transformation.center{path(it+1),path(it)};
	[rPath(:,2*(it)),flagInfeas1] = func_commandGovernor(YAL_v,YAL_opts,MOAS_A,MOAS_b,nv,rCur,xCur,vCur);

	MOAS_A = rooms.MOAS.restrict.A{path(it),path(it+1)}*HMat;
	MOAS_b = rooms.MOAS.restrict.b{path(it),path(it+1)}-1e-4; %tightining;
	xCur = [];
	rCur = rPath(:,2*it);
	vCur = rooms.transformation.center{path(it),path(it+1)};

	[rPath(:,2*(it-1)+1),flagInfeas2] = func_commandGovernor(YAL_v,YAL_opts,MOAS_A,MOAS_b,nv,rCur,xCur,vCur);

	if flagInfeas1||flagInfeas2
		sprintf("problemInOptimization when generating Path")
	end
end
end


function [freeBoxes, bool_smallBoxes] = func_cleanFreeBoxes(freeBoxes,minD_and_V,minD_alone)
dim_freeBoxes = freeBoxes*kron(eye(3),[-1;1]);
minDim = min(dim_freeBoxes,[],2);
boxVol = cumprod(dim_freeBoxes,2);boxVol = boxVol(:,3);
bool_smallBoxes = boolean((minDim<=minD_and_V(1)).*(boxVol<=minD_and_V(2))+(minDim<=minD_alone));
freeBoxes(bool_smallBoxes,:) = [];
end

% plot functions

function cur_fig = func_showRooms(rooms,roomIdx,cur_fig,faceColor,faceAlpha)
if nargin<5 || isempty(faceAlpha)
	faceAlpha = .3;
end
if nargin<4 || isempty(faceColor)
	faceColor = 'y';
end
if nargin<3 || isempty(cur_fig)
	cur_fig = figure;
else
	hold on;
end
if nargin<2 || isempty(roomIdx)
	roomIdx = 1:rooms.n;
end
for it= 1:length(roomIdx)
	jt = roomIdx(it);
	temp_pol = polytope(rooms.closedRoom.HRep.A{jt},rooms.closedRoom.HRep.b{jt});
	temp_pol = Polyhedron(temp_pol.vertices');
	temp_pol.plot;hold on;
	cur_ax = findobj(cur_fig,'Type','axes');
	curObj = cur_ax.Children(1);
	curObj.FaceColor = faceColor;
	curObj.FaceAlpha = faceAlpha;
	curObj.EdgeColor = [0,0,1];
	curObj.EdgeAlpha = faceAlpha;
end
hold off
end

function cur_fig = func_showObstacles(roomsPre,roomIdx,cur_fig,faceColor,faceAlpha)
if nargin<5 || isempty(faceAlpha)
	faceAlpha = .3;
end
if nargin<4 || isempty(faceColor)
	faceColor = 'y';
end
if nargin<3 || isempty(cur_fig)
	cur_fig = figure;
else
	hold on;
end
if nargin<2 || isempty(roomIdx)
	roomIdx = 1:length(roomsPre.obstacles.HRep.A);
end
for it= 1:length(roomIdx)
	jt = roomIdx(it);
	temp_pol = polytope(roomsPre.obstacles.HRep.A{jt},roomsPre.obstacles.HRep.b{jt});
	temp_pol = Polyhedron(temp_pol.vertices');
	temp_pol.plot;hold on;
	cur_ax = findobj(cur_fig,'Type','axes');
	curObj = cur_ax.Children(1);
	curObj.FaceColor = faceColor;
	curObj.FaceAlpha = faceAlpha;

end
hold off
end

function cur_fig = func_showWeakExtensions(rooms,roomIdx,cur_fig,faceColor,faceAlpha)
if nargin<5 || isempty(faceAlpha)
	faceAlpha = .3;
end
if nargin<4 || isempty(faceColor)
	faceColor = 'y';
end
if nargin<3 || isempty(cur_fig)
	cur_fig = figure;
else
	hold on;
end
if nargin<2 || isempty(roomIdx)
	roomIdx = 1:rooms.n;
end
for it= 1:length(roomIdx)-1
	jt = roomIdx(it);
	zt = roomIdx(it+1);
	temp_pol = polytope(rooms.weakExtension.HRep.A{jt,zt},rooms.weakExtension.HRep.b{jt,zt});
	temp_pol = Polyhedron(temp_pol.vertices');
	temp_pol.plot;hold on;
	cur_ax = findobj(cur_fig,'Type','axes');
	curObj = cur_ax.Children(1);
	curObj.FaceColor = faceColor;
	curObj.FaceAlpha = faceAlpha;
	curObj.EdgeColor = [0 0 1];
	curObj.EdgeAlpha = faceAlpha;
end
hold off;
end