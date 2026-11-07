function [u, h, xi] =  solve_mnid(allA, allB, x0, T, rho, uL)

    if nargin < 6   
        uL = 5;
    end
	n = size(allA, 2);
	p = size(allB, 2);
	N = floor(size(allA, 1)/n);

	assert(size(allA, 1) == N*n);
	assert(size(allB, 1) == N*n);
	assert(length(x0) == n);

	LPA = [zeros(2*n*N, p*T) -ones(2*n*N, 1) -kron(eye(N),ones(2*n,1))];   % LP matrix
	LPb = zeros(2*n*N, 1);       	% LP vector
	LPc = [zeros(p*T, 1); 1; rho*ones(N,1)];       	% LP vector
	% Ax <= b

	
	for i = 1:N,
		A = allA((i-1)*n + 1:i*n, :);
		B = allB((i-1)*n + 1:i*n, :);

		Apower = eye(n);
		for t = 1:T,
			Z = Apower*B;
			Y = [Z; -Z];

			LPA( (i-1)*2*n + 1:i*2*n , (t-1)*p + 1:t*p ) = Y;			

			Apower = A*Apower;
		end;

		Z = Apower*x0;
		Y = [-Z; Z];
		LPb( (i-1)*2*n + 1 : i*2*n ) = Y;			
	end;
	    
    % Now SOLVE 
	% min LPc*uh
	% s.t. LPA*uh <= LPc
           
    sol = cplexlp(LPc,LPA,LPb,[],[],[-uL*ones(p*T,1); 0; zeros(N,1)],[uL*ones(p*T,1); Inf; Inf*ones(N,1)]);

    
    u = sol(1:p*T);
    h = sol(p*T+1);
    xi = sol(p*T+2:end);
      


end

