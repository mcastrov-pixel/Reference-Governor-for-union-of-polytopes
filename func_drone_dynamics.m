function sys = func_drone_dynamics()
% system name
sys.id ='drone_';
% dimensions
sys.nx = 6;
sys.nv = 3;
sys.nu = 3;

% sampling time
sys.dt = .1;

% 3 times double integrator (drone dynamics after exact feedback
% linearization)
sys.A = [zeros(3) eye(3);zeros(3,6)];
sys.B = [zeros(3);eye(3)];
% Outputs are given by each state and input (reference command)
sys.CCstr = [eye(sys.nx);zeros(sys.nv,sys.nx)]; % extracts position states x1 and x2
sys.DCstr = [zeros(sys.nx,sys.nv);eye(sys.nv)];

% discretize
dtSys = c2d(ss(sys.A,sys.B,sys.CCstr,sys.DCstr),sys.dt);
sys.Adt = dtSys.A;
sys.Bdt = dtSys.B;
sys.Cdt = dtSys.C;
sys.Ddt = dtSys.D;

% stabilize the system using LQR
sys.Q = blkdiag(eye(3)*1e1,eye(3)*1e-1);
sys.R = eye(3)*1e-2;

[sys.K,sys.P]= dlqr(sys.Adt,sys.Bdt,sys.Q,sys.R);

% compute closed loop matrices A_cl = A - BK; B_cl = B*K*v2x*p2v where
% sys.v2x is the original mapping from ref to state and p2v is a
% parameterization of the reference command so that it corresponds to
% any position

sys.Adt = sys.Adt-sys.Bdt*sys.K;
sys.Bdt=  sys.Bdt*sys.K*eye(sys.nx,sys.nv);
sys.v2x = (eye(6)-sys.Adt)\sys.Bdt;
sys.v2x(abs(sys.v2x)<1e-15) = 0;
sys.x2v = (sys.v2x'*sys.v2x)\sys.v2x'; % left pseudo inverse
sys.x2v(abs(sys.x2v)<1e-15) = 0;

% parameterize the reference to equal the position
sys.Bdt =sys.Bdt*sys.x2v(1:sys.nv,1:sys.nv);
sys.v2x = eye(sys.nx,sys.nv);
sys.x2v = eye(sys.nv,sys.nx);

% compute the v and x to torq input matrices (used for input
% constraints)
sys.v2u = sys.K*eye(sys.nx,sys.nv)*sys.x2v(1:sys.nv,1:sys.nv);
sys.x2u = -sys.K;

sys.Adt(abs(sys.Adt)<=1e-15) = 0;
sys.Bdt(abs(sys.Bdt)<=1e-15) = 0;
end