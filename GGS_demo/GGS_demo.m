%% initialize GGS-data
ax_steps = 20;  % number of discretization points
GGS = struct();
GGS.ax = zeros(ax_steps,1);
GGS.ay = zeros(ax_steps,1);
GGS.delta = zeros(ax_steps,1);
GGS.beta = zeros(ax_steps,1);
GGS.Power = zeros(ax_steps,1);
%% Vehicle Parameters:
m = 250;    % vehicle mass (kg)
lf = 0.8;   % distance COG to front axle (m)
lr = 0.8;   % distance COG to rear axle (m)
% Tire Parameters
A =1800; B =1.5; C =25; D =1; E =20;
% maximum Power
P_max = 80000;      % (W)
% Velocity
Vel = 15;           % (m/s)
%% find max. longitudinal acceleration
nlp = []; vehicle_model_opti;
% objective
nlp.minimize(-ax);
% initialization of decision variables
nlp.set_initial(delta,0);
nlp.set_initial(beta,0);
nlp.set_initial(dpsi,0);
nlp.set_initial(Sxf,0);
nlp.set_initial(Sxr,0);
% power constraint
nlp.subject_to(Power<=P_max);
% solve
nlp.solver('ipopt');
sol = nlp.solve();
% extract results
GGS.ay(1) = sol.value(ay);
GGS.ax(1) = sol.value(ax);
GGS.delta(1) = sol.value(delta);
GGS.beta(1) = sol.value(delta);
GGS.Power(1) = sol.value(Power);
%% find max. longitudinal deceleration
nlp = []; vehicle_model_opti;
% objective
nlp.minimize(ax);
% initialization of decision variables
nlp.set_initial(delta,0);
nlp.set_initial(beta,0);
nlp.set_initial(dpsi,0);
nlp.set_initial(Sxf,0);
nlp.set_initial(Sxr,0);
% solve
nlp.solver('ipopt');
sol = nlp.solve();
% extract results
GGS.ay(end) = sol.value(ay);
GGS.ax(end) = sol.value(ax);
GGS.delta(end) = sol.value(delta);
GGS.beta(end) = sol.value(beta);
GGS.Power(end) = sol.value(Power);
%% find max. steady state lateral acceleration
GGS.ax = linspace(GGS.ax(1),GGS.ax(end),length(GGS.ax));

for i_ = 2:length(GGS.ax)-1
    nlp = []; vehicle_model_opti;
    % objective
    nlp.minimize(-ay);
    % initialization of decision variables
    nlp.set_initial(delta,2*pi/180);
    nlp.set_initial(beta,-1*pi/180);
    nlp.set_initial(dpsi,1);
    nlp.set_initial(Sxf,0);
    nlp.set_initial(Sxr,0);
    % steady state constraints
    nlp.subject_to(Mz == 0);
    nlp.subject_to(ay-Vel*dpsi == 0);
    % longitudinal acceleration constraint
    nlp.subject_to(ax == GGS.ax(i_));
    % solve
    nlp.solver('ipopt');
    sol = nlp.solve();
    % extract results
    GGS.ay(i_) = sol.value(ay);
    GGS.delta(i_) = sol.value(delta);
    GGS.beta(i_) = sol.value(beta);
    GGS.Power(i_) = sol.value(Power);
end

%% plot steering angle
figure('Color','w');
plot(GGS.delta*180/pi,GGS.ax,'k-d','Linewidth',1);hold on;
grid on;
xlabel('steering angle (deg)','FontSize',14);ylabel('A_x (m/s^2)','FontSize',14);
%% plot sideslip angle
figure('Color','w');
plot(GGS.beta*180/pi,GGS.ax,'k-d','Linewidth',1);hold on;
grid on;
xlabel('sideslip angle (deg)','FontSize',14);ylabel('A_x (m/s^2)','FontSize',14);
%% plot Power output
figure('Color','w');
plot(GGS.Power/1000,GGS.ax,'k-d','Linewidth',1);hold on;
grid on;
xlabel('Power output (kW)','FontSize',14);ylabel('A_x','FontSize',14);
%% plot GG - Diagram
figure('Color','w');
plot(GGS.ay,GGS.ax,'k-d','Linewidth',1);hold on;
grid on;
xlabel('A_y (m/s^2)','FontSize',14);ylabel('A_x (m/s^2)','FontSize',14);
daspect([1,1,1]);