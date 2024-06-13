function [X, dx, ddx, tt] = generate_spd_data_sim(X, tau, dt, plot_manifold)
% Generates SPD trajectory of a given points.
% transport all projections to a common tangent space, (e.g. of the first
% point).

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

nbData = round(1 / dt);
t = linspace(dt, 1, nbData);
nbData=100;t = linspace(dt, tau, nbData);

dims = 2;
pdims = dims + dims * (dims - 1)/2;

if(plot_manifold)
    f1 = figure;

    colorss = lines(20);
    for i=1:length(X)
        if (i < length(X)), g12(:,:,:,i) = geodesic(X(:,:,i),X(:,:,i+1),t); end
        [Hp1, ~] = plotGMM2([i; 0], X(:,:,i), colorss(1,:),.1);hold on
        %[Hp2, ~] = plotGMM2([50; 0], dD(j).KO(:,:,i), colorss(2,:),.1);
    end

    %% Plot SPD cone space
    f2 = plot_spd_cone(150);

    %% Plot geodesics from S1 to S2
    nbDrawingSeg = 60; %Number of segments used to draw ellipsoids
    for i=1:length(X)-1
        plot_geodesic(X(:,:,i),X(:,:,i+1),nbData);
    end

    %% Plot datapoints on the manifold
    ll = floor(linspace(1,nbData,15));nn = length(ll);

    for i=1:length(X)
        pdata = plot3(reshape(X(1,1,i),[1,1]), reshape(X(2,2,i),[1,1]), reshape(X(2,1,i),[1,1]), '.','markersize',15,'color',[0.1 0.1 0.1]);
    end
    pdata = plot3(reshape(X(1,1,1),[1,1]), reshape(X(2,2,1),[1,1]), reshape(X(2,1,1),[1,1]), '^','markersize',15,'MarkerFaceColor',colorss(end,:),'color',colorss(end,:));
    pdata = plot3(reshape(X(1,1,end),[1,1]), reshape(X(2,2,end),[1,1]), reshape(X(2,1,end),[1,1]), '^','markersize',15,'MarkerFaceColor',colorss(floor(length(colorss)/2),:),'color',colorss(end/2,:));
end
%%

t_geo = 0:0.05:1;
dg12(:,:,1) = eye(2)*0; clear U;
Zi(:,:,1) = X(:,:,1);
Zii(:,:,1) = X(:,:,1);
spd_ref = eye(2);
for i = 2:length(X)
    
    dX(:,:,i) = (logm(X(:,:,i)) - logm(X(:,:,i-1)))/dt;
    dX(:,:,i) = logmap(X(:,:,i),X(:,:,i-1))/dt;
    
    dX(:,:,i) = transp_vec(X(:,:,i-1),X(:,:,1),dX(:,:,i));
    
    dx(:,i) = symMat2Vec(dX(:,:,i));
    
    Up = (dx(:,i) + symMat2Vec(X(:,:,i-1)));
    
    if(plot_manifold)
        pdata = [pdata plot3(Up(1,:), Up(2,:), Up(3,:), '.','markersize',10,'color',[.8 0 0])];hold on;
        pdata = [pdata plot3([X(1,1,i-1) Up(1)], [X(2,2,i-1) Up(2)], [X(2,1,i-1) Up(3)], '-','linewidth',1,'color',[.8 0 0]) ];
    end
    
	%Zi(:,:,i) = expm(dX(:,:,i)*dt + logm(Zi(:,:,i-1)));
	%Zi(:,:,i) = expmap(dX(:,:,i)*dt,Zi(:,:,i-1));
	Zi(:,:,i) = expmap(transp_vec(X(:,:,1),X(:,:,i-1),dX(:,:,i))*dt,Zi(:,:,i-1));
end



dx(:,1) = dx(:,2);
% Calculate derivatives
tt = linspace(dt*1, tau, length(X));
for j = 1:3
    ddx(j,:) = gradient(dx(j,:), tt);
end
ddx(:,1) = ddx(:,3);
ddx(:,2) = ddx(:,3);

if(plot_manifold)
    for i=1:length(X)
        pdata = plot3(reshape(Zi(1,1,i),[1,1]), reshape(Zi(2,2,i),[1,1]), reshape(Zi(2,1,i),[1,1]), 'o','markersize',10,'color',[0.2 0.4 0.7]);
    end

    for t=50:nbData
        gmr_c = plot3(ddx(1,t),ddx(2,t),ddx(3,t), '.','markersize',25,'color',[0.5 .0 .8]);
    end
end

end
