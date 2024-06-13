function L = logmap(S, D)
% Logarithm map (SPD manifold)
N = size(D,3);
L = 0*D;
for n = 1:N
    L(:,:,n) = S^.5 * logm(S^-.5 * D(:,:,n) * S^-.5) * S^.5;
% 	[v,d] = eig(S\X(:,:,n));
% 	U(:,:,n) = S * v*diag(log(diag(d)))*v^-1;
end
end