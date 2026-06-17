function sys = func_CWH_dynamics()

% system name
sys.id ='CWH_';

%physical parameters
sys.mu = 398600.4; % gravitational constant for earth km^3/s^2
sys.Re = 6378;     % Earth radius (km)
sys.Ra = 400; % altitude km
sys.r0 = sys.Re+sys.Ra;         %km
sys.n = sqrt(sys.mu/sys.r0^3);

% dimensions
sys.nx = 6;
sys.nv = 3;
sys.nu = 3;

% sampling time
sys.dt = 2;

% 3 times double integrator (drone dynamics after exact feedback
% linearization)
sys.A = [zeros(3,3) eye(3);
        3*sys.n^2 zeros(1,3) 2*sys.n 0; %radial
        zeros(1,3) -2*sys.n 0 0; %in track
        zeros(1,2) -sys.n^2 zeros(1,3); %cross track
        ];

% reorder state so that the spacecraft is not oddly aligned with respect to
% orbit
sys.A = blkdiag([0 0 1; 0 1 0; 1 0 0],[0 0 1; 0 1 0; 1 0 0])*sys.A*blkdiag([0 0 1; 0 1 0; 1 0 0],[0 0 1; 0 1 0; 1 0 0]);



sys.B = [zeros(3,3);eye(3)];

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
sys.Q = blkdiag(eye(3)*1e3,eye(3)*1e-1);
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