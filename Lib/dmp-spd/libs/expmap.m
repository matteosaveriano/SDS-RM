function D = expmap(S, L)
% Exponential map (SPD manifold)
N = size(L,3);
D = 0*L;
for n = 1:N
    D(:,:,n) = S^.5 * expm(S^-.5 * L(:,:,n) * S^-.5) * S^.5;
end
end