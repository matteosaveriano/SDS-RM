function plot_manib_data_cone(S,XX)
% Plots SPD data in SPD cone
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

nbData = length(S);
%% Plot SPD cone space
manCone = plot_spd_cone(150);
colorss = lines(20);

%% Plot geodesics from S1 to S2
nbDrawingSeg = 60; %Number of segments used to draw ellipsoids
for i=1:length(S)-1
    o_geo = plot_geodesic(S(:,:,i),S(:,:,i+1),1);
end
for i=1:length(XX)-1
    n_geo = plot_geodesic(XX(:,:,i),XX(:,:,i+1),1,[0 1 0]);
end

%% Plot datapoints on the manifold
ll = floor(linspace(1,nbData,15));nn = length(ll);

for i=1:length(S)
    odata1 = plot3(reshape(S(1,1,i),[1,1]), reshape(S(2,2,i),[1,1]), reshape(S(2,1,i),[1,1]), '.','markersize',10,'color',[0.1 0.1 0.1]);
end
odata2 = plot3(reshape(S(1,1,1),[1,1]), reshape(S(2,2,1),[1,1]), reshape(S(2,1,1),[1,1]), '^','markersize',10,'MarkerFaceColor',colorss(end,:),'color',colorss(end,:));
odata3 = plot3(reshape(S(1,1,end),[1,1]), reshape(S(2,2,end),[1,1]), reshape(S(2,1,end),[1,1]), '^','markersize',10,'MarkerFaceColor',colorss(floor(length(colorss)/2),:),'color',colorss(end/2,:));


for i=1:length(XX)
    pdata1 = plot3(reshape(XX(1,1,i),[1,1]), reshape(XX(2,2,i),[1,1]), reshape(XX(2,1,i),[1,1]), '*','markersize',10,'color',[0.2 0.4 0.7]);
end
pdata2 = plot3(reshape(XX(1,1,1),[1,1]), reshape(XX(2,2,1),[1,1]), reshape(XX(2,1,1),[1,1]), '^','markersize',15,'MarkerFaceColor',colorss(end,:),'color',colorss(end,:));
pdata3 = plot3(reshape(XX(1,1,end),[1,1]), reshape(XX(2,2,end),[1,1]), reshape(XX(2,1,end),[1,1]), '^','markersize',15,'MarkerFaceColor',colorss(floor(length(colorss)/2),:),'color',colorss(end/2,:));

h_leg = legend([o_geo(1) n_geo(1) odata1(1) odata2(1) odata3(1)],...
   '$\rho^{orig}$', '$\rho^{dmp}$','$\Upsilon_i^{orig}$',...
   '$\Upsilon_1^{orig}$', '$\Upsilon_{100}^{orig}$', 'Location','northeast');
set(h_leg,'Interpreter','latex','FontSize',11);
h_leg = legend([o_geo(1) n_geo(1)],...
   '$\Upsilon$', '$\hat{\Upsilon}$', 'Location','northeast');
set(h_leg,'Interpreter','latex','FontSize',10);

end