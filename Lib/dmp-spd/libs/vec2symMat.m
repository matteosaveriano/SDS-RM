function S = vec2symMat(V)
% Transforms matrix of vectors to tensor of symmetric matrices

[d, N] = size(V);
D = (-1 + sqrt(1 + 8*d))/2;
for n = 1:N
	v = V(:,n);
	M = diag(v(1:D));
	id = cumsum(fliplr(1:D));

	for i = 1:D-1
%          M = M + diag(v(id(i)+1:id(i+1)),i) + diag(v(id(i)+1:id(i+1)),-i);
	  M = M + diag(v(id(i)+1:id(i+1)),i)./sqrt(2) + diag(v(id(i)+1:id(i+1)),-i)./sqrt(2); % Mandel notation
% 	  M = M + diag(v(id(i+1)+1:id(i+1)),i) + diag(v(id(i+1)+1:id(i+1)),-i); % Voigt notation
	end
	S(:,:,n) = M;
end
end