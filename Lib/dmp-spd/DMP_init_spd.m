%DMP_INIT_SPD initiate DMP parameters 

% INPUTS:
%   -X0: initial state of DMPs (quaternion position)
%   -X_goal: goal state of DMPs
%   -dt: sample time
%   -N: number of nonlinear radial basis functions per DMP
%   -alpha_qGoal: a gain needed in case of switching the goal
%   -alpha_z: stiffness gain
%   -alpha_x: a decay factor of the discrete canonical system
%   -alpha_pPos: a gain needed in case of position phase stopping
%   -alpha_pOri: a gain needed in case of orientation phase stopping

% OUTPUTS:
%   -DMP: a structure defining the DMP and its parameters

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

function DMP = DMP_init_spd(X0, X_goal,dt,N,...
    alpha_qGoal,alpha_z,alpha_x,alpha_pPos,alpha_pOri)

if nargin<3
   error('At Least Three Input Arguments are Required.')
end

if nargin==3
    N = 20;
    alpha_qGoal = 0.05;
    alpha_z = 48/3;
    alpha_x = 2;
    alpha_pPos = 400;
    alpha_pOri = 4;
elseif nargin==4
    alpha_qGoal = 0.5;
    alpha_z = 48/3;
    alpha_x = 2;
    alpha_pPos = 400;
    alpha_pOri = 4;
elseif nargin==5
    alpha_z = 48/3;
    alpha_x = 2;
    alpha_pPos = 400;
    alpha_pOri = 4;
elseif nargin==6
    alpha_x = 2;
    alpha_pPos = 400;
    alpha_pOri = 4;
elseif nargin==7
    alpha_pPos = 400;
    alpha_pOri = 4;
elseif nargin==8
    alpha_pOri = 4;
end

DMP.N = N;
DMP.dt = dt;
DMP.X0 = X0;
DMP.X_goal = X_goal;
DMP.alpha_qGoal = alpha_qGoal;
DMP.alpha_z = alpha_z;
DMP.beta_z = DMP.alpha_z / 4;
DMP.alpha_x = alpha_x;
DMP.alpha_pPos = alpha_pPos;
DMP.alpha_pOri = alpha_pOri;

%% Gaussian kernel functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
c_lin = linspace(0, 1, DMP.N);
DMP.c = exp(-DMP.alpha_x * c_lin);
DMP.sigma2 = (diff(DMP.c)*0.75).^2;
DMP.sigma2 = [DMP.sigma2, DMP.sigma2(end)];
DMP.w = zeros(3, DMP.N);

end