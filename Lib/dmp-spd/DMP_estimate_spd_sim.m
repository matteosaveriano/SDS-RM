function [w, logspd_out, f] = DMP_estimate_spd_sim(X, dx, ddx, t, DMP)
% Learns SPD data profile
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
%-------------------------------------------------------------------------
Do = symMat2Vec(logmap(DMP.X_goal,DMP.X0));
d = 1./Do;

N = numel(t);
f = zeros(N,3);
A = zeros(N, DMP.N);
logspd_out = [];
for i = 1:N
    
    log_spd = symMat2Vec(transp_vec(X(:,:,i),DMP.X0,...
        logmap(DMP.X_goal,X(:,:,i))));
    
    logspd_out =  [logspd_out, log_spd];
    f(i,:) = d.*(DMP.tau^2*ddx(:,i) + ...
        DMP.alpha_z*DMP.tau*dx(:,i) - ...
        DMP.alpha_z*DMP.beta_z * log_spd);
    
    x = exp(-DMP.alpha_x*t(i)/DMP.tau);
    psi = exp(-(x-DMP.c).^2./(2*DMP.sigma2));
    A(i,:) = x*psi/sum(psi);
end
w = (A \ f)';

end
