function p_geo = plot_geodesic(S1,S2,nbData,colors)

% Plot geodesics from S1 to S2

if nargin<4
    colors = [0.6 0.6 0.6];
end

nbDrawingSeg = 60; %Number of segments used to draw ellipsoids
p_geo = [];
U = logmap(S1,S2);
%Up = U + repmat(S2,[1,1,nbData]);
umsh = bsxfun(@times,U,reshape(linspace(0,1,nbDrawingSeg),1,1,60));
msh = expmap(umsh, S2);
p_geo = [p_geo plot3(reshape(msh(1,1,:),[nbDrawingSeg,1]),...
    reshape(msh(2,2,:),[nbDrawingSeg,1]), reshape(msh(2,1,:),...
    [nbDrawingSeg,1]), '-','linewidth',3,'color',colors)];


end