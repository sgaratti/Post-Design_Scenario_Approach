function [xT, x_inf_norm] =  test_mnid(allA,allB,x0,T,u);
	
	n = size(allA, 2);
	p = size(allB, 2);
	N = floor(size(allA,1)/n);

	assert(size(allA,1) == N*n);
	assert(size(allB,1) == N*n);
	assert(length(x0) == n);
	assert(length(u) == p*T);

	R = zeros(n*N,p*T);
	S = zeros(n*N,n);   
	
	for i = 1:N
		A = allA((i-1)*n + 1:i*n, :);
		B = allB((i-1)*n + 1:i*n, :);

		Apower = eye(n);
		for t = 1:T
			R( (i-1)*n+1:i*n , (t-1)*p+1:t*p ) = Apower*B;
			Apower = A*Apower;
		end
		S( (i-1)*n+1:i*n , : ) = Apower;
	end
	x_vec = R*u + S*x0;
	
	xT = reshape(x_vec,n,N); 
	x_inf_norm = max(abs(xT),[],1);
end
