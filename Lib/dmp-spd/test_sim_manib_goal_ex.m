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
SS.X_new_goal = DMP.X_goal;
 
R = [cosd(30) -sind(30); sind(30) cosd(30)];

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
 
i = 1;
 
log_spd = [];
FX=[];
phase=[];
%% Run DMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while SS.x > exp(-DMP.alpha_x * (DMP.tau+0 + DMP.dt) / DMP.tau)
    if i==50
        SS.X_new_goal = R'*SS.X_goal*R;
    end
    [SS, log_s, fx] = DMP_integrate_spd_sim_goal(DMP, SS);
    phase = [phase, SS.x];
    XX(:,:,i) = SS.X;
    dXX(:,i) = SS.dx;
    ddXX(:,i) = SS.ddx;
    log_spd = [log_spd, log_s];
    FX = [FX,fx];

    i = i + 1;
end;
 

kk = length(XX);
nnnnn = XX(:,:,1:kk);
clrmap = summer(5);

%% Plot distances
stifDistGoal=figure('position',[10 5 390 200]);hold on
for i=1:50
    d_le1(i) = distance_Log_Euclidean(nnnnn(:,:,i),S(:,:,i));
end
for i=50:kk
    d_le2(i-50+1) = distance_Log_Euclidean(nnnnn(:,:,i),SS.X_new_goal);
end
yyaxis left;
plot(1:50,d_le1,'color','b');
ylabel('$distance$', 'Fontsize', 10, 'Interpreter', 'Latex');
yyaxis right
plot(50:kk,d_le2,'color','r');%plot(d_KL);plot(d_ld);plot(d_ri);
plot(50,linspace(0,2.3,50),'.','color','g','markersize',5,'linewidth',1.5);
ylabel('$distance$', 'Fontsize', 10, 'Interpreter', 'Latex');
yyaxis left;
text(50,0.005,'$\leftarrow$ Seperator','interpreter','latex','FontName','Times','Fontsize',12)
h = legend('$d_1$','$d_2$');
set(h, 'Fontsize', 10, 'Interpreter', 'Latex');
xlabel('$t$', 'Fontsize', 10, 'Interpreter', 'Latex');
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 10);
%print(stifDistGoal,'-dpng','-r600','results/stifDistGoal');
%print(stifDistGoal,'-depsc2','-r600','results/stifDistGoal');
%print(stifDistGoal,'-dsvg','-r600','results/stifDistGoal');


%% Plot stiffness over Cartesian
stifCartGoal=figure('position',[10 10 390 300]);hold on
for i=1:kk
    [Hp1, ~] = plotGMM2(XT(:,i), 5E-2*S(:,:,i), [.6 .6 .6], .4); % Scaled matrix!
end
for i=round(linspace(1,kk,20))
    [Hp2, ~] = plotGMM2(XT(:,i), 5E-2*nnnnn(:,:,i), [0 1 0], .4); % Scaled matrix!
end
axis equal;
%[Hp3, ~] = plotGMM2(msd(1).DataP(1:2,50), 5E-5*nnnnn(:,:,50), [0 0 1], .4); % Scaled matrix!
%[Hp4, ~] = plotGMM2(msd(1).DataP(1:2,end), 5E-5*nnnnn(:,:,end), [1 0 0], .4); % Scaled matrix!
%h = legend([Hp1(1) Hp2(1) Hp4(1)], '$\mathbf{{K}}^\mathcal{P}$',...
%    '$\mathbf{\hat{K}}^\mathcal{P}$', '$\mathbf{{K}}_{newGoal}^\mathcal{P}$', 'Location','northwest');
%set(h, 'Fontsize', 10, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 10, 'Interpreter', 'Latex');
xlabel('$x$', 'Fontsize', 10, 'Interpreter', 'Latex');
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 10);
%print(stifCartGoal,'-dpng','-r600','results/stifCartGoal');
%print(stifCartGoal,'-depsc2','-r600','results/stifCartGoal');
%print(stifCartGoal,'-dsvg','-r600','results/stifCartGoal');

%% Plot stiffness over Time
stifTimeGoal=figure('position',[10 10 570 140]);hold on
for i=1:100
    [Hp2, ~] = plotGMM2([i;0], 5E-2*S(:,:,i), [.6 .6 .6], .4); % Scaled matrix!
end
for i=round(linspace(1,kk,20))
    [Hp1, ~] = plotGMM2([i;0], 5E-2*nnnnn(:,:,i), [0 1 0], .4); % Scaled matrix!
end
[Hp3, ~] = plotGMM2([50;0], 5E-2*nnnnn(:,:,50), [0 0 1], .4); % Scaled matrix!
[Hp4, ~] = plotGMM2([100;0], 5E-2*nnnnn(:,:,end), [1 0 0], .4); % Scaled matrix!
axis([-5, 105, -6, 6]);
ylabel('$\mathbf{{K}}^\mathcal{P}$ and $\mathbf{\hat{K}}^\mathcal{P}$', 'Fontsize', 10, 'Interpreter', 'Latex');
xlabel('$t$', 'Fontsize', 10, 'Interpreter', 'Latex');
set(gca,'TickLabelInterpreter','latex', 'Fontsize', 10);
% print(stifTimeGoal,'-dpng','-r600','results/stifTimeGoal');
% print(stifTimeGoal,'-depsc2','-r600','results/stifTimeGoal');
% print(stifTimeGoal,'-dsvg','-r600','results/stifTimeGoal');

%% Plot Velocity 
stifVelGoal=figure('position',[10 10 450 170]);hold on
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
%print(stifVelGoal,'-dpng','-r600','results/stifVelGoal');
%print(stifVelGoal,'-depsc2','-r600','results/stifVelGoal');
%print(stifVelGoal,'-dsvg','-r600','results/stifVelGoal');

%% Plot Velocity and acceleration
stifDerivGoal=figure('position',[10 10 450 370]);hold on
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
axis([0, 100, -300, 100]);
set(gca,'TickLabelInterpreter','latex');
%print(stifDerivGoal,'-dpng','-r600','results/stifDerivGoal');
%print(stifDerivGoal,'-depsc2','-r600','results/stifDerivGoal');
%print(stifDerivGoal,'-dsvg','-r600','results/stifDerivGoal');



