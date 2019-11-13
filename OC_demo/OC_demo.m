close all;clear;clc
%% Vehicle Parameters:
m = 250;    % vehicle mass (kg)
I_z = 120;  % yaw inertia (kgm²)
lf = 0.8;   % distance COG to front axle (m)
lr = 0.8;   % distance COG to rear axle (m)
% Tire Parameters
A =1800; B =1.5; C =25; D =1; E =20;
% maximum Power
P_max = 80000;      % (W)
%% Track parameters
step_length = 1;
straight_length = 20;
turn_length = 50;
min_radius = 8;
Track = Generate_LinCurv_turn(step_length,straight_length,turn_length,min_radius);

%% Vehicle dynamics
% Define states
n = casadi.SX.sym('n');         % orthogonal path deviation (m)
xi = casadi.SX.sym('xi');       % heading angle deviation (rad)
u = casadi.SX.sym('u');         % vehicle fixed x-velocity (m/s)
v = casadi.SX.sym('v');         % vehicle fixed y-velocity (m/s)
dpsi = casadi.SX.sym('dpsi');   % vehicle yaw rate (rad/s)
x_ir = casadi.SX.sym('x_ir');   % x-Position in global coordinates (m)
y_ir = casadi.SX.sym('y_ir');   % y-Position in global coordinates (m)
psi = casadi.SX.sym('psi');     % yaw angle (rad)

states = [n;xi;u;v;dpsi;x_ir;y_ir;psi];
nx = length(states);
x_init = [0;0;10;0;0;0;0;0];

% Define control inputs
delta = casadi.SX.sym('delta'); % steering angle (rad)
Sxf = casadi.SX.sym('Sxf');     % front long. slip (-)
Sxr = casadi.SX.sym('Sxr');     % rear long. slip (-)

controls = [delta;Sxf;Sxr];
nu = length(controls);
u_init = [0;0;0];

% Curvature variable
Crv = casadi.SX.sym('Crv');

% front/rear slip angles
Saf = -delta + atan((v + lf*dpsi)/u);
Sar = atan((v - lr*dpsi)/u);

% tire forces
% pure slip
Fxpf = A*sin(B*atan(C*Sxf));
Fypf = -A*sin(B*atan(C*tan(Saf)));
Fxpr = A*sin(B*atan(C*Sxr));
Fypr = -A*sin(B*atan(C*tan(Sar)));
% combined slip
Fxf = Fxpf * cos(D*atan(E*tan(Saf)));
Fyf = Fypf * cos(D*atan(E*Sxf));
Fxr = Fxpr * cos(D*atan(E*tan(Sar)));
Fyr = Fypr * cos(D*atan(E*Sxr));

% Dynamics Scaling factor
Sf = (1-n*Crv)/(u*cos(xi)-v*sin(xi));

% Spatial dynamics
rhs = casadi.SX.sym('rhs',nx);
rhs(1) = Sf * (u*sin(xi) + v*cos(xi));
rhs(2) = Sf * dpsi - Crv;
rhs(3) = Sf * (dpsi*v + 1/m *(Fxf*cos(delta) - Fyf*sin(delta) + Fxr));
rhs(4) = Sf * (-dpsi*u + 1/m *(Fyf*cos(delta) + Fxf*sin(delta) + Fyr));
rhs(5) = Sf * 1/I_z*(lf*(Fyf*cos(delta) + Fxf*sin(delta)) - lr*Fyr);
rhs(6) = Sf * (u*cos(psi) - v*sin(psi));
rhs(7) = Sf * (u*sin(psi) + v*cos(psi));
rhs(8) = Sf * dpsi;

f_rhs = casadi.Function('f_rhs',{states,controls,Crv},{rhs});

% approximation of power output
Power =  Fxf*(u*cos(delta)+v*sin(delta))*(1+Sxf) + Fxr*u*(1+Sxr);
f_Power = casadi.Function('f_Power',{states,controls},{Power});

% Runge-Kutta 4 integration
M = 2; % RK4 steps per interval
ds = step_length/M;
X = states;
for i_ = 1:M
    k1 = f_rhs(X,controls,Crv);
    k2 = f_rhs(X+ds/2*k1,controls,Crv);
    k3 = f_rhs(X+ds/2*k2,controls,Crv);
    k4 = f_rhs(X+ds*k3,controls,Crv);
    X = X + ds/6*(k1+2*k2+2*k3+k4);
end
fX = casadi.Function('fX',{states,controls,Crv},{X});

%% build optimization problem
N = Track.N-1;    % number of discretization intervals
% Decision variables (controls)
U = casadi.SX.sym('U',nu,N); 
% Decision variables (states)
X = casadi.SX.sym('X',nx,N+1); 

% dynamics equality constraints
g = X(:,1)-x_init;  
for k = 1:N
    g = [g;X(:,k+1) - fX(X(:,k),U(:,k),Track.curv(k))];
end

% drivetrain power limitation
for k = 1:N
    g = [g;f_Power(X(:,k),U(:,k)) - P_max];
end

obj = 0;
for k = 1:N+1
    obj = obj + ((1-X(1,k)*Track.curv(k))/(X(3,k)*cos(X(2,k))-X(4,k)*sin(X(2,k))))^2;
end

Qu = diag([1e-3,1e-2,1e-2]);
Qdu = diag([1e-1,1e-2,1e-2]);

obj = obj + U(:,1).'*Qu*U(:,1) + (U(:,1)-u_init).'*Qdu*(U(:,1)-u_init);
for k = 2:N
    obj = obj + U(:,k).'*Qu*U(:,k) + (U(:,k)-U(:,k-1)).'*Qdu*(U(:,1)-U(:,k-1));
end

% make the decision variables one column vector
OPT_variables = [reshape(U,nu*N,1);reshape(X,nx*(N+1),1)];
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g);

opts = struct;
% opts.print_time = 0;
% opts.ipopt.max_iter = 100;
% opts.ipopt.print_level =0;%0,3
% opts.ipopt.acceptable_tol =1e-8;
% opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = casadi.nlpsol('solver', 'ipopt', nlp_prob,opts);

% constraints setup
args = struct;
% upper and lower bounds on g
args.lbg = zeros(nx*(N+1) + N,1);  
args.ubg = zeros(nx*(N+1) + N,1);
args.lbg(nx*(N+1)+1:end) = -inf;
% upper and lower bounds on decision variables
args.lbx = -inf * ones(nu*N + nx*(N+1),1);
args.ubx = inf * ones(nu*N + nx*(N+1),1);
args.lbx(1:nu:nu*N) = -30*pi/180;
args.ubx(1:nu:nu*N) = 30*pi/180;
args.lbx(2:nu:nu*N) = -0.5;
args.ubx(2:nu:nu*N) = 0.5;
args.lbx(3:nu:nu*N) = -0.5;
args.ubx(3:nu:nu*N) = 0.5;
args.ubx(1) = 0;
args.lbx(1) = 0;
args.ubx(nu*(N-1)+1) = 0;
args.lbx(nu*(N-1)+1) = 0;
args.lbx(nu*N+1:nx:end) = -1;
args.ubx(nu*N+1:nx:end) = 1;
args.x0 = zeros(nu*N + nx*(N+1),1);
args.x0(nu*N+3:nx:end) = x_init(3);
args.x0(nu*N+6:nx:end) = Track.x;
args.x0(nu*N+7:nx:end) = Track.y;
args.x0(nu*N+8:nx:end) = Track.psi;

%% Solve OCP
sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg);
x_sol = full(sol.x);

% extract result
res.U = x_sol(1:nu*N);
res.U = reshape(res.U,nu,N);
% res.U = [res.U,res.U(:,end)];
res.X = x_sol(nu*N+1:end);
res.X = reshape(res.X,nx,N+1);

% calculate laptime
res.dt = step_length * (1-res.X(1,:).*Track.curv)./(res.X(3,:).*cos(res.X(2,:))-res.X(4,:).*sin(res.X(2,:)));
res.time = [0,cumsum(res.dt(1:N))];
res.laptime = res.time(end);
%% plots

%% path
figure('Color','w');
plot(res.X(6,:),res.X(7,:));hold on;
plot(Track.x,Track.y);grid on;
xlabel('x(m)'),ylabel('y(m)','FontSize',14);
legend({'vehicle path','reference path'},'FontSize',14);
daspect([1,1,1])

%% velocities
figure('Color','w');
plot(Track.S,res.X(3,:));hold on;
plot(Track.S,res.X(4,:));grid on;
xlabel('distance (m)','FontSize',14),ylabel('velocity (m/s)','FontSize',14);
legend({'v_x','v_y'},'FontSize',14);

%% normal path deviation
figure('Color','w');
plot(Track.S,res.X(1,:));grid on;
xlabel('distance (m)','FontSize',14),ylabel('n(m)','FontSize',14);

%% steering angle
figure('Color','w');
plot(Track.S(1:length(res.U(1,:))),res.U(1,:)*180/pi);grid on;
xlabel('distance(m)','FontSize',14),ylabel('steering angle (deg)','FontSize',14);

%% longitudinal slip
figure('Color','w');
plot(Track.S(1:length(res.U(2,:))),res.U(2,:));hold on;
plot(Track.S(1:length(res.U(3,:))),res.U(3,:));grid on;
xlabel('distance (m)','FontSize',14),ylabel('longitudinal slip','FontSize',14);
legend({'front slip','rear slip'},'FontSize',14);

%% yaw angle
figure('Color','w');
plot(Track.S,res.X(8,:)*180/pi);hold on;
plot(Track.S,Track.psi*180/pi);grid on;
xlabel('distance(m)','FontSize',14),ylabel('psi(deg)','FontSize',14);
legend({'vehicle','reference'},'FontSize',14);

%% power output
% approximation of power output
for k = 1:N
    res.Power(k) =  full(f_Power(res.X(:,k),res.U(:,k)));
end

figure('Color','w');
plot(Track.S(1:length(res.Power)),res.Power); grid on
xlabel('distance (m)','FontSize',14); ylabel('Power (W)','FontSize',14);
