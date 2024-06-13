function v = symMat2Vec(S)
% Reduced vectorisation of a symmetric matrix.
% [d,~,N] = size(S);
% 
% v = zeros(d+d*(d-1)/2,N);
% for n = 1:N
%     v(1:d,n) = diag(S(:,:,n));
%     
%     row = d+1;
%     for i = 1:d-1
%         v(row:row + d-1-i,n) = sqrt(2).*S(i+1:end,i,n);
%         row = row + d-i;
%     end
% end

% Vectorization of a tensor of symmetric matrix
[D, ~, N] = size(S);

V = [];
for n = 1:N
	v = diag(S(:,:,n));
	for d = 1:D-1
	 v = [v; sqrt(2).*diag(S(:,:,n),d)]; % Mandel notation
% 	   v = [v; diag(M,n)]; % Voigt notation
%         v = [v; diag(S(:,:,n),d)];
	end
	V = [V v];
end

end