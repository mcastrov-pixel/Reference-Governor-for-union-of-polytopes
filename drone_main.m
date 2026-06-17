
%Main file for drone tracking maneuver
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
sys = func_drone_dynamics();
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
xdotLim = 1;
xdot_Ay = kron( eye(sys.nx/2),[-1;1]); % if only velocity output
xdot_Ay = [zeros(size(xdot_Ay,1),3), xdot_Ay,zeros(size(xdot_Ay,1),sys.nv)]; % increase to include the other outputs of the sys.tem ( pos and input)
xdot_by = xdotLim*ones(sys.nx/2*2,1);

%2.2)input constraints
g = 9.81; % gravitational acceleration m/s^2
Tmax = 1.4*g; % max normalized thrust available
epmax = 15*pi/180; % this means that the drone can only lean no more than 10 degree
s1=15;s2 = 10; %for approximation of the constraint set (~ number of vertices)
[V_input,torq_u_Ay,torq_u_by] = func_drone_approximateInputConstraint(s1,s2,Tmax,epmax); % if output are given by the input to the sys.tem
[u_relError, u_V_real, u_V_poly] = coneSpherePolytopeVolumeError(Tmax, epmax, s1, s2)
sprintf("mismatch between polytopic approximation and original constraint set %d percents",u_relError)
torq_Ay = torq_u_Ay*[sys.x2u sys.v2u];

%2.3)assemble it together
other_Ay = [torq_Ay;xdot_Ay];
other_by = [torq_u_by;xdot_by];

% generate urban environment: obstacles / free space

[buildingBoxes, freeBoxes, domain] = func_drone_space_decomposition(35, true, 3);
roomsPre = func_hyperrectangle_bounds_to_struct(buildingBoxes, freeBoxes, domain);
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
	save(strcat(name_base,'MOAS'),'rooms','roomsPre',"buildingBoxes","freeBoxes")
else
	load(strcat(name_base,'MOAS'));
end

%% Simulation

% generate the desired setpoint and find associated element in sequence
rDes = [65,27,25]';
%find box corresponding to desired set point
isInBox = [freeBoxes(:,1:2:6)-rDes', rDes'-freeBoxes(:,2:2:6)];
isInBox = all(isInBox<=0,2);
list_endElem = find(isInBox); list_endElem = list_endElem(1);


% generate the initial conditions: here we make initial conditions that start on the ground and at the top of each building 

% ground locations
list_startElem= 1:rooms.n;
list_startElem(list_endElem) = [];
list_startElem(freeBoxes(list_startElem,5)~=0) = [];
x0_list = zeros(sys.nx,size(list_startElem,2));
for it=1:length(list_startElem)
	roomStart = list_startElem(it);
	x0_list(:,it) = sys.v2x*(rooms.transformation.center{roomStart,roomStart}.*[1;1;1e-1]);
end
%top of building locations
p0_list = buildingBoxes*kron(diag([ 1 1 0]),[1/2;1/2]);
p0_list(:,3) = buildingBoxes(:,6)+1;
x0_list = [x0_list,sys.v2x*p0_list'];
for it =1:size(p0_list,1)
	isInBox = [freeBoxes(:,1:2:6)-p0_list(it,:), p0_list(it,:)-freeBoxes(:,2:2:6)];
	isInBox = all(isInBox<0,2);
	tempBox = find(isInBox);
	if isempty(tempBox)
		keyboard
	end
	list_startElem(end+1) = tempBox;
end
list_startEndPairs = [kron(list_startElem, ones(size(list_endElem))); kron(ones(size(list_startElem)),list_endElem)];

% comment out to run simmulations from all initial conditions
selectedSims = [1 67 size(list_startEndPairs,2)];
list_startEndPairs = list_startEndPairs(:,selectedSims);
x0_list = x0_list(:,selectedSims);






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
	path_refs = func_genPathSetPoints_ditsLine(rooms,path_rooms,x0(1:3,1),sys.v2x(1:3,:)*rDes,sys.v2x);
	flagInfeas = [];



	% simulation length depends on # of rooms to go through
	NSim = max(round(150*path_length),1500);

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
	if mod(zt,10)== 9
		save(strcat(name_base,'temp_DATA'),'DATA','sys');
	end
end
save(strcat(name_base,'DATA'),'DATA','sys');


%% trajectory plot

close (figure(12))
cur_fig = func_plotCityDecomposition(buildingBoxes, [], domain,12);
hold on;
%plot all trajectories
for it = 1:length(DATA.xHist)
	handle_traj = plot3(DATA.xHist{it}(1,:),DATA.xHist{it}(2,:),DATA.xHist{it}(3,:),'Color',[0,1,1,0.5],'LineWidth',1.5);hold on

	xStart = DATA.xHist{it}(:,1);
	handle_IC = plot3(xStart(1),xStart(2),xStart(3),'g>','LineWidth',6);
	if it ==length(DATA.xHist)
		xEnd = sys.v2x*rDes;
		hold on;
		handle_rDes = plot3(xEnd(1,1),xEnd(2,1),xEnd(3),'m<','LineWidth',2,'MarkerSize',10, 'MarkerFaceColor', 'auto');
		hold off;
	end
end
hold off;


%% plot velocities


figure(67);clf;

selectedTraj = [1 2 3];
plot([0 80],[xdotLim xdotLim],'k-.','LineWidth',2.7);hold on;
cellfun(@(y) plot((1:length(y))*sys.dt,max(abs(y(4:6,:)),[],1),'LineWidth',1.7),DATA.xHist(selectedTraj),'uniformOutput',false)
hold off;
ylabel('$\|vel\|_\infty$ [m/s]','FontSize',18,'interpreter','latex')
xlabel('time [s]','FontSize',18,'interpreter','latex')
xlim([0 80])

%% plot input


	cur_fig = figure(78); clf
	%plot the approx input constraint set + input traj

	temp_pol = Polyhedron('V',V_input);
	temp_pol.plot();hold on;
	cur_ax = findobj(cur_fig,'Type','axes');
	handle_cstrSet= cur_ax.Children(1);
	% handle_cstrSet.FaceVertexCData = min(max([0.8 0.06 0.6].*randn(size(temp_pol.V,1),3),[0 0 0]),[1 1 1]);
	handle_cstrSet.FaceColor = [0.3 0.3 0.3];
	handle_cstrSet.FaceAlpha = .3;
	handle_cstrSet.EdgeAlpha = 1;
	handle_cstrSet.LineWidth = 0.5;
	light('Position', [0, 0, 10]);
	camlight headlight
	lighting gouraud;
	hold on;
	for jt = 1:length(selectedTraj)
		%compute the torque input
		it = selectedTraj(jt);
		uVals = sys.x2u*DATA.xHist{it}(:,1:end)+sys.v2u*DATA.vHist{it};
		plot3( uVals(1,:),uVals(2,:),uVals(3,:),'LineWidth',1.6);
	end
	xlabel('x-acc. [m/s$^2$]','FontSize',18,'interpreter','latex'); xticks([-4:2:4] );
	ylabel('y-acc. [m/s$^2$]','FontSize',18,'interpreter','latex');yticks([-4:2:4] );
	zlabel('z-acc. [m/s$^2$]','FontSize',18,'interpreter','latex');zticks([-10:5:5] );
	hold off
	grid on
%% Input plot (how far from polytope boundary

figure(67);clf;

maxCstrVal = {}; %compute max A*u-b
	% winterColor = autumn(length(selectedTraj));
	for jt = 1:length(selectedTraj)
		%compute the torque input
		it = selectedTraj(jt);
		uVals = sys.x2u*DATA.xHist{it}(:,1:end)+sys.v2u*DATA.vHist{it};
		maxCstrVal{jt} = max(torq_u_Ay*uVals -torq_u_by,[],1);

	end
handle_uLim= plot([0 80],[0 0],'k-.','LineWidth',2.7);hold on;
cellfun(@(y) plot((1:length(y))*sys.dt,y,'LineWidth',1.7),maxCstrVal,'uniformOutput',false)
hold off;
ylabel('$\max(A^u u(t) - bu)$','FontSize',18,'interpreter','latex')
xlabel('time [s]','FontSize',18,'interpreter','latex')
xlim([0 80])

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

%used to check how good of an approx of the input constraints we have
function [rel_error, V_real, V_poly] = coneSpherePolytopeVolumeError(Tmax, epmax, s1, s2)

    % Number of unique angular sectors.
    % This assumes alpha = linspace(0,2*pi,s1), so the first and last
    % directions are duplicates.
    n = s1 - 1;

    % Number of radial samples
    m = s2;

    dtheta = 2*pi/n;

    % Normalized radial samples
    a = linspace(0, sin(epmax), m);
    h = sqrt(1 - a.^2);

    radial_sum = 0;

    for j = 1:m-1
        radial_sum = radial_sum + ...
            (a(j) + a(j+1)) * ...
            (a(j+1)*h(j) - a(j)*h(j+1));
    end

    % Real volume
    V_real = (2*pi/3)*Tmax^3*(1 - cos(epmax));

    % Polytope volume
    V_poly = (n*sin(dtheta)/6)*Tmax^3*radial_sum;

    % Relative missing volume
    rel_error = (V_real - V_poly)/V_real;
end


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

function rPath = func_genPathSetPoints_ditsLine(rooms,path,pt1,pt2,v2x);
eps = 1e-5;
nv = length(pt2);
nPath = length(path);
rPath = zeros(nv,2*(nPath-1)+1);
rPath(:,end) = pt2;
YAL_x = sdpvar(nv,1);
YAL_opts = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
HMat = [eye(nv);v2x];
v = pt2-pt1;
Qform = eye(3) - v*v'/(v'*v);
for it = flip(1:nPath-1)
	% set up the line
	YAL_cstr = rooms.MOAS.restrict.A{path(it+1),path(it)}*HMat*YAL_x-rooms.MOAS.restrict.b{path(it+1),path(it)}<=-1e-4;
	YAL_obj = (YAL_x-pt1)'*Qform*(YAL_x-pt1);
	YAL_diagno = optimize(YAL_cstr,YAL_obj,YAL_opts);
	if ~YAL_diagno.problem
		rPath(:,2*(it)) = value(YAL_x);
	else
		problemInOptimization =1
	end
	YAL_cstr = rooms.MOAS.restrict.A{path(it),path(it+1)}*HMat*YAL_x-rooms.MOAS.restrict.b{path(it),path(it+1)}<=-1e-4;
	YAL_obj = (YAL_x - rPath(:,2*it))'*(YAL_x - rPath(:,2*it));
	YAL_diagno = optimize(YAL_cstr,YAL_obj,YAL_opts);
	if ~YAL_diagno.problem
		rPath(:,2*(it-1)+1) = value(YAL_x);
	else
		problemInOptimization =1
	end
end
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