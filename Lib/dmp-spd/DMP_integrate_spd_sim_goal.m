
function [S, log_spd, fx] = DMP_integrate_spd_sim_goal(DMP, S)
% Generates SPD data profile with goal switching
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

Do = symMat2Vec(logmap(S.X_goal,DMP.X0));

% phase variable
dx = -DMP.alpha_x * S.x;
dx = dx / DMP.tau;
S.x = S.x + dx * DMP.dt;

% the weighted sum of the locally weighted regression models
if S.x >= exp(-DMP.alpha_x)
    psi = exp(-(S.x-DMP.c).^2./(2*DMP.sigma2));
    fx = sum( sum((DMP.w(:,:)*S.x).*psi,2)/sum(psi),3);
else
    fx = [0; 0; 0];
end

% integration of angular velocity
log_spd = symMat2Vec(transp_vec(S.X,DMP.X0,logmap(S.X_goal,S.X)));

S.ddx = DMP.alpha_z * (DMP.beta_z * log_spd - S.dx) + fx.*Do;

% temporal scaling
S.ddx = S.ddx / DMP.tau;
S.dx = S.dx + S.ddx * DMP.dt;

%parallel transport dx to S.x pefor expm.
S.X = expmap(transp_vec(DMP.X0,S.X,vec2symMat(S.dx))/ DMP.tau* DMP.dt, S.X);

% Goal Switching: smooth transition to the new goal.
do = DMP.alpha_qGoal * logmap(S.X_new_goal,S.X_goal);
if norm(do) > 1.0e-12
    S.X_goal = expmap(do, S.X_goal);
end

end