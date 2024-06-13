function ptv = transp_vec(S, D, T)
% 	if nargin == 3
% 		t = 1;
% 	end
	%U = logmap(S2,S1);
	%ptv = S1^.5 * expm(0.5.*t.*S1^-.5*U*S1^-.5) * S1^-.5 * v * S1^-.5 * expm(0.5.*t.*S1^-.5*U*S1^-.5) * S1^.5;
	% Computationally economic way : ptv = (S2/S1)^.5 * v * ((S2/S1)^.5)';
    ptv = (D/S)^.5 * T * ((D/S)^.5)';
end
