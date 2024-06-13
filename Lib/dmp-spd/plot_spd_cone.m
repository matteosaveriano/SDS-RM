function f2 = plot_spd_cone(r)

if ~nargin
    r = 20;  % or whatever
end
%% Plot SPD cone
f2 = figure('position',[10,10,950,350],'color',[1 1 1]); hold on; axis off; axis equal;

% Plot the SPD convex cone
phi = 0:0.1:2*pi+0.1;
x = [zeros(size(phi)); r.*ones(size(phi))];
y = [zeros(size(phi));r.*sin(phi)];
z = [zeros(size(phi));r/sqrt(2).*cos(phi)];

h = mesh(x,y,z,'linestyle','none','facecolor',[.95 .95 .95],'facealpha',.95);
direction = cross([1 0 0],[1/sqrt(2),1/sqrt(2),0]);
rotate(h,direction,45,[0,0,0])

h = plot3(x(2,:),y(2,:),z(2,:),'linewidth',2,'color',[0 0 0]);
rotate(h,direction,45,[0,0,0])
h = plot3(x(:,63),y(:,63),z(:,63),'linewidth',2,'color',[0 0 0]);
rotate(h,direction,45,[0,0,0])
h = plot3(x(:,40),y(:,40),z(:,40),'linewidth',2,'color',[0 0 0]);
rotate(h,direction,45,[0,0,0])
% Settings
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
axis off
view(70,12);


end