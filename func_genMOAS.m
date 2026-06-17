function [A_Oinf,b_Oinf] = func_genMOAS(A,B,C,D,Ay,by,epsVal)
%System considered is
% x_+ = A x + B v,
% and constraints (Ay y <= by) are on the ouptut defined as
% y = C x + D v

% The resulting moas is given as
% A_Oinf [v' x']' <= b_Oinf


nCstr = size(Ay,1);
nx = size(A,1);

% secondary parameters
toll    = 1e-4; % tolerance parameter
Inx = eye(nx);


% --- Constraints on steady state commands (strictly sterady state!!!), i.e., k = \infty
A_Oinf = Ay*[D+C*((Inx-A)\B), 0*C];
b_Oinf = by*(1-epsVal);
% [A_Oinf, b_Oinf] = elimRedundant(A_Oinf, b_Oinf, toll);

% --- Constraints on predicted state at time k=1,2,..,kfin

flg_done = 0;

YAL_x = [];
YAL_opts = [];
YAL_x = sdpvar(nx+size(B,2),1);
YAL_opts = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
YAL_opts.mosek.MSK_DPAR_ANA_SOL_INFEAS_TOL = 1e-10;

it = -1;
prevDecr = 1;
decrGoal = 0.1;
while ~flg_done
	it = it +1;
	AOi_pre = A_Oinf;
	bOi_pre = b_Oinf;
	% build constraints on current state for O_{\infty,0}
	A_Oinf= [A_Oinf;
		Ay*[D + C*((Inx - A)\(Inx - A^it))*B, C*A^it]];
	b_Oinf =  [b_Oinf;by];

	% delete redundant constraints + check if we have finalized computation of the O_infty
	if it >0 && mod(it,1)==0
		if ~isempty(AOi_pre)
			if issubset(AOi_pre, bOi_pre, A_Oinf, b_Oinf, toll,YAL_x,YAL_opts)
				flg_done = 1;
				popol = polytope(A_Oinf,b_Oinf);
				popol = popol.compact;
				(length(b_Oinf) - length(popol.b))/length(b_Oinf)
				A_Oinf = popol.A;
				b_Oinf = popol.b;
				% [A_Oinf, b_Oinf] = elimRedundant(A_Oinf, b_Oinf, 1e-6,YAL_x,YAL_opts);
			end
			if ~flg_done && mod(it,20) == 10
				% [A_Oinf, b_Oinf] = elimRedundant(A_Oinf, b_Oinf, 1e-6,YAL_x,YAL_opts);
				it
				popol = polytope(A_Oinf,b_Oinf);
				popol = popol.compact;
				(length(b_Oinf) - length(popol.b))/length(b_Oinf)
				A_Oinf = popol.A;
				b_Oinf = popol.b;
			end
		end

	end
	if it >= 1000 % not necessarily a problem but For this simple MOAS I'm not expecting to need more than 1000 it.
		prob = 1
	end
end
end

function [AOut, bOut,IDX] = elimRedundant(AIn, bIn, toll,YAL_x,YAL_opts)
% eliminates redundant constraints from the polyhedron defined as {x | Ain x<= bin} by solving small
% linear problems.

nCstr=size(AIn,1);
options=optimset('Display','off');

IDX = [2:nCstr]; % active constraints
n_del = 0;
for it=1:nCstr
	YAL_cstr = AIn(IDX,:)*YAL_x<=bIn(IDX);
	YAL_obj = -AIn(it,:)*YAL_x;
	YAL_diag = optimize(YAL_cstr,YAL_obj,YAL_opts);
	% [hx,~,EXITFLAG] = linprog(-AIn(it,:),AIn(IDX,:),bIn(IDX),[],[],[],[],options);
	% if EXITFLAG == -3 % meaning problem is unbounded
	%     hxval = Inf;
	% elseif EXITFLAG >=1
	%     hxval= AIn(it,:)*hx;
	% else
	%     prob=2;
	%     hxval = Inf;
	% end
	hxval = AIn(it,:)*value(YAL_x);
	IDX(it-n_del) =  it;
	if hxval<bIn(it)+toll
		IDX(it-n_del) = [];
		it;
		n_del = n_del+1;
	end
	if n_del+1 == nCstr
		IDX = nCstr;
		break;
	end
end

AOut = AIn(IDX,:);
bOut = bIn(IDX);
end

function flg_isSubset = issubset(A1, b1,A2,b2, toll,YAL_x,YAL_opts)
% verifies whether the polyhedron defined by A2 b2 is inside the polyhedron defined with A1, b1. Outputs
% a flag.
% Polyhedrons are defined as {x| Ax \leq b}


nCstr =size(A2,1); % number of inequalitites which define A2
flg_isSubset=true; %
options=optimset('Display','off');
toBeTested = length(b2)-length(b1); %assume most cstr are the same
S1 = polytope(A1,b1);
S2 = polytope(A2,b2);

	try
	flg_isSubset = contains(S2,S1,'approx',1e-2);
	catch
	for it=nCstr-toBeTested+1:nCstr
		f=-A2(it,:);   %select the ith constraint, i=1,...,nCstr, of A1, b1
		YAL_cstr = (A1./b1)*YAL_x<=b1*0+1;
		YAL_obj = f*YAL_x;
		YAL_diag = optimize(YAL_cstr,YAL_obj,YAL_opts);
		% [hx_linprog,~,EXITFLAG] = linprog(f,A1,b1,[],[],[],[],options); %solve an lp problem to determine if redundant
		% if EXITFLAG == -3 % meaning problem is unbounded
		%     remain_linprog = Inf;
		% elseif EXITFLAG >=1
		%     remain_linprog = -f*hx_linprog-b2(it); %
		% else
		%     prob=3
		%     remain_linprog = 3;
		% end
		if YAL_diag.problem == 2
			remain = Inf;
		elseif YAL_diag.problem == 0
			remain = -f*value(YAL_x)-b2(it);
		else
			stopPoint = 1;
		end
		if remain>=toll,
			% it - nCstr
			% remain
			% S1 = polytope(A1,b1);S2 = polytope(A2,b2);
			% a = contains(S1,S2);
			% a
			flg_isSubset=false; break,
		end
		% if (abs(remain)<Inf) && (abs(remain_linprog-remain)>1e-4)
		% 	remain_linprog
		% 	remain
		% 	it
		% else
		% 	it
		% end
	end
end
end
