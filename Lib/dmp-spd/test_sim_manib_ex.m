clear all
close all
clc
% TOY EXAMPLE
%
% Please acknowledge the authors if this code or part of it is useful for 
% your research.
% @inproceedings{abudakka2020Geometry,
%	title		= {Geometry-aware dynamic movement primitives},
% 	author		= {Abu-Dakka, Fares J. and Kyrki, Ville},
% 	booktitle	= {IEEE International Conference on Robotics and Automation},
% 	pages		= {4421--4426},
% 	address		= {Paris, France},
% 	year		= {2020}
% }
%
% Author: Fares J. Abu-Dakka
% Intelligent Robotics Group, Aalto University
% Website: https://sites.google.com/view/abudakka/
% Sept. 2019; Last update: Feb. 2021

path(path,'libs');
colorss = lines(20);



dt = 1;
nbData = 1/dt;
t = linspace(0,dt*nbData,nbData);

load data/man_gmrCn.mat
for i=1:length(man_gmr)
    M = (triu(ones(2)));
    M(logical(M)) = man_gmr(i).mu(1:3);
    S(:,:,i) = M' * M;
    gmrRes(:,:,i) = S(:,:,i);
    %S(:,:,i) = X(2:3,2:3,i);
end
%tau = length(man_gmr);
tau = 100

dt=.1;tau=10;
[S, ds, dds, t] = generate_spd_data_sim(S, tau, dt);

dt=t(2)-t(1);
N = 40;
alpha_qGoal = 0.5;
alpha_z = 48*2;
alpha_x = 2;
alpha_pPos = 400;
alpha_pOri = 4;
DMP = DMP_init_spd(S(:,:,1), S(:,:,end), dt, N,...
    alpha_qGoal,alpha_z,alpha_x,alpha_pPos,alpha_pOri);
DMP.tau = tau;

%% DMP estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [DMP.w, logspd_orig, f ] = DMP_estimate_spd_sim(S, ds, dds, t, DMP);
 
 %% Initial states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
SS.X = DMP.X0;
SS.X0 = DMP.X0;
SS.X_goal = DMP.X_goal;
 
R = [cosd(90) -sind(90); sind(90) cosd(90)];
SS.X_new_goal = R'*SS.X_goal*R;
% Q = zeros(round((DMP.tau-duration)/DMP.dt), 4);
XX = [];
% W = zeros(round((DMP.tau-duration)/DMP.dt), 3);
dXX = [];
% dW = zeros(round((DMP.tau-duration)/DMP.dt), 3);
ddXX = [];
% t1 = zeros(round((DMP.tau-duration)/DMP.dt), 1);
% x = zeros(round((DMP.tau-duration)/DMP.dt), 1);
SS.dx =  ds(:,1); zeros(3,1);
SS.x = 1;
XX(:,:,1) = SS.X;
dXX(:,1) = SS.dx*DMP.tau;
ddXX(:,1) = dds(:,1);
% Q(1,:) = [S.q.s, S.q.v'];
% W(1,:) = S.o';
dW(1,:) = [0, 0, 0];
t1(1) = 0;
x(1) = 1;
 
i = 2;
 
log_spd = [];
FX=[];
phase=[];
%% Run DMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while SS.x > exp(-DMP.alpha_x * (DMP.tau + DMP.dt) / DMP.tau)
    if i==47
        disp('49')
    end
    [SS, log_s, fx] = DMP_integrate_spd_sim(DMP, SS);
    phase = [phase, SS.x];
    XX(:,:,i) = SS.X;
    dXX(:,i) = SS.dx;
    ddXX(:,i) = SS.ddx;
    log_spd = [log_spd, log_s];
    FX = [FX,fx];

    i = i + 1;
end
 

kk = 100;
nnnnn = XX(:,:,1:kk);

%% Plot distances
manDist=figure('position',[10 10 390 200]);hold on
for i=1:kk
    d_le(i) = distance_Log_Euclidean(nnnnn(:,:,i),S(:,:,i));
end
plot(d_le);
%h = legend('Log-euclidean distance', 'Log Determinant distance');
%set(h, 'Fontsize', 10, 'Interpreter', 'Latex');
ylabel('$distance$', 'Fontsize', 10, 'Interpreter', 'Latex');
xlabel('$t$', 'Fontsize', 10, 'Interpreter', 'Latex');
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 10);
% print(manDist,'-dpng','-r600','results/manDist');
% print(manDist,'-depsc2','-r600','results/manDist');
% print(manDist,'-dsvg','-r600','results/manDist');


%% Plot stiffness over Cartesian
manCart=figure('position',[10 10 390 300]);hold on
for i=1:100
    %xt = [man_gmr(i).t(2),man_gmr(i).t(3)]';
    [Hp1, ~] = plotGMM2(XT(:,i), 5E-2*S(:,:,i), [.6 .6 .6], .4); % Scaled matrix!
end
for i=round(linspace(1,kk,20))
    %xt = [man_gmr(i).t(2),man_gmr(i).t(3)]';
    [Hp2, ~] = plotGMM2(XT(:,i), 5E-2*nnnnn(:,:,i), [0 1 0], .4); % Scaled matrix!
end
axis equal;
h = legend([Hp1(1) Hp2(1)], '$\Upsilon$', '$\hat{\Upsilon}$', 'Location','northwest');
set(h, 'Fontsize', 10, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 10, 'Interpreter', 'Latex');
xlabel('$x$', 'Fontsize', 10, 'Interpreter', 'Latex');
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 10);
% print(manCart,'-dpng','-r600','results/manCart');
% print(manCart,'-depsc2','-r600','results/manCart');
% print(manCart,'-dsvg','-r600','results/manCart');

%% Plot stiffness over Time
manTime=figure('position',[10 10 450 170]);hold on
for i=1:100
    [Hp2, ~] = plotGMM2([i;0], 5E-1*S(:,:,i), [.6 .6 .6], .4); % Scaled matrix!
end
for i=round(linspace(1,kk,20))
    [Hp1, ~] = plotGMM2([i;0], 5E-1*nnnnn(:,:,i), [0 1 0], .4); % Scaled matrix!
end
axis([-10, 105, -10, 10]);
ylabel('$\Upsilon$ and $\hat{\Upsilon}$', 'Fontsize', 10, 'Interpreter', 'Latex');
xlabel('$t$', 'Fontsize', 10, 'Interpreter', 'Latex');
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 10);
% print(manTime,'-dpng','-r600','results/manTime');
% print(manTime,'-depsc2','-r600','results/manTime');
% print(manTime,'-dsvg','-r600','results/manTime');

%% Plot Velocity 
manVel=figure('position',[10 10 450 170]);hold on
v1 = plot(ds','k--');
v2 = plot(dXX'/DMP.tau,'r'); 
set(get(get(v1(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(v1(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(v2(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h = legend('ORIG','DMP');
set(h, 'Fontsize', 10, 'Interpreter', 'Latex', 'Location','northwest');
ylabel('1st time derivative', 'Fontsize', 10, 'Interpreter', 'Latex');
xlabel('$t$', 'Fontsize', 10, 'Interpreter', 'Latex');
axis([0, 100, -320, 310]);
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 10);
% print(manVel,'-dpng','-r600','results/manVel');
% print(manVel,'-depsc2','-r600','results/manVel');
% print(manVel,'-dsvg','-r600','results/manVel');

%% Plot Velocity and acceleration
manDeriv=figure('position',[10 10 450 370]);hold on
subplot(211);hold on
v1 = plot(ds','k--');
v2 = plot(dXX'/DMP.tau,'r'); 
set(get(get(v1(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(v1(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(v2(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h = legend('ORIG','DMP');
set(h, 'Fontsize', 10, 'Interpreter', 'Latex', 'Location','northwest');
ylabel('1st time derivative', 'Fontsize', 12, 'Interpreter', 'Latex');
xlabel('$t$', 'Fontsize', 12, 'Interpreter', 'Latex');
axis([0, 100, -320, 310]);
set(gca,'TickLabelInterpreter','latex');

subplot(212);hold on
a1 = plot(dds','k--');
a2 = plot(ddXX'/DMP.tau,'r'); 
ylabel('2nd time derivative', 'Fontsize', 12, 'Interpreter', 'Latex');
xlabel('$t$', 'Fontsize', 12, 'Interpreter', 'Latex');
axis([0, 100, -900, 1000]);
set(gca,'TickLabelInterpreter','latex');
% print(manDeriv,'-dpng','-r600','results/manDeriv');
% print(manDeriv,'-depsc2','-r600','results/manDeriv');
% print(manDeriv,'-dsvg','-r600','results/manDeriv');



