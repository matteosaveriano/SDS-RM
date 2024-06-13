function g = geodesic(S1,S2,t)
g = zeros(size(S2,1),size(S2,1),length(t));


% Shp = S2^.5;
% Smp = S2^-.5;
% for i=1:length(t)
% 	g(:,:,i) = Shp * expm(t(i) .* Smp * S1 * Smp) * Shp;
% end
%g(:,:,1) = S1;
for i=1:length(t)
	
% 	%Interpolation between more than 2 covariances can be computed in an iterative form
% 	for n=1:nbIter
% 		W = zeros(model.nbVar);
% 		for i=1:model.nbStates
% 			W = W + w(i,t) * logmap(model.Sigma(:,:,i), S);
% 		end
% 		S = expmap(W,S);
% 	end
% 	Sigma(:,:,t) = S;

	%Interpolation between two covariances can be computed in closed form
	g(:,:,i) = expmap(t(i)*logmap(S2, S1), S1);
	
end

end

%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = expmap(W,S)
	S = S^.5 * expm(S^-.5 * W * S^-.5) * S^.5;
end

function S = logmap(W,S)
	S = S^.5 * logm(S^-.5 * W * S^-.5) * S^.5;
end
