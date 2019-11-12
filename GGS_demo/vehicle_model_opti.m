%% initialize Problem
nlp = casadi.Opti();
% decision variables & box constraints
delta = nlp.variable(); nlp.subject_to(-30*pi/180<=delta<=30*pi/180);   % steering angle (rad)
beta = nlp.variable(); nlp.subject_to(-20*pi/180<=beta<=20*pi/180);     % sideslip angle (rad)
dpsi = nlp.variable(); nlp.subject_to(-90*pi/180<=dpsi<=90*pi/180);     % yaw rate (rad/s)
Sxf = nlp.variable(); nlp.subject_to(-0.2<=Sxf<=0.2);                   % front longitudinal slip (-)
Sxr = nlp.variable(); nlp.subject_to(-0.2<=Sxr<=0.2);                   % rear longitudinal slip (-)

%% vehicle dynamics equations (dynamic single track model)
% velocities in vehicle fixed coordinates
dx = Vel*cos(beta);
dy = Vel*sin(beta);

% front/rear slip angles
Saf = -delta + atan((dy + lf*dpsi)/dx);
Sar = atan((dy - lr*dpsi)/dx);

% tire model
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

% sum of forces in vehicle fixed coordinates
Fy = Fyf*cos(delta) + Fxf*sin(delta) + Fyr;
Fx = Fxf*cos(delta) - Fyf*sin(delta) + Fxr;
Mz = lf*(Fyf*cos(delta) + Fxf*sin(delta)) - lr*Fyr;

% accelerations in path tangential coordinates
ay = 1/m * (Fy*cos(beta) - Fx*sin(beta));
ax = 1/m * (Fy*sin(beta) + Fx*cos(beta));

% approximation of power output
Power =  Fxf*(dx*cos(delta)+dy*sin(delta))*(1+Sxf) + Fxr*dx*(1+Sxr);