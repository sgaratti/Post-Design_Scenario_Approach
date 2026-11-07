rng(5032025) %reproducibile run


N=2000;

n=2;
g = 9.8;
Ml=9; Mu=10;
Ll=0.9; Lu=1;
dM = Ml+(Mu-Ml)*rand(N,1);
dL = Ll+(Lu-Ll)*rand(N,1);
da = dM.*(1+0.1*rand(N,1));

h=1;


Den=[ones(N,1) , da./(dM.*dL) , - g./dL];
Num=[zeros(N,1), zeros(N,1), 1./(dM.*dL) ];

Pref = poly([ [-1 -1 -3 -3]+[1i -1i 3i -3i ]]);

w = ones(1,2*n);
f = [zeros(1,2*n) w];
A = zeros(N*4*n,2*n+2*n);
B = zeros(N*4*n,1);
X = zeros(N*2*n,2*n);
Y = zeros(N*2*n,1);

for ix=1:N
    auxAD = toeplitz([Den(ix,:) zeros(1,n-1)],[Den(ix,1) zeros(1,n-1)]);
    auxAN = toeplitz([Num(ix,:) zeros(1,n-1)],[Num(ix,1) zeros(1,n-1)]);
    auxBD = [Den(ix,2:n+1) zeros(1,n)]';
    auxBR = Pref(1,2:2*n+1)';
    
    A((ix-1)*4*n+1:ix*4*n,:) = [[auxAD auxAN; -auxAD -auxAN] [-eye(2*n); -eye(2*n)]];
    B((ix-1)*4*n+1:ix*4*n) = [auxBR-auxBD; auxBD-auxBR];
    X((ix-1)*2*n+1:ix*2*n,:) = [auxAD auxAN];
    Y((ix-1)*2*n+1:ix*2*n) = auxBD;
end

[sol,fval,exitflag,output,lambda] = linprog(f,A,B);
PC = [ones(N,1) reshape(X*sol(1:2*n)+Y,[2*n N])'];

% % % %complexity analysis 
is_scenario=sum(reshape(lambda.ineqlin,4*n,N),1)>0;
support_indices=find(is_scenario==1);
complexity=length(support_indices)


%Safety check: the solution with support_indices only must be the same as before
    fred = [zeros(1,2*n) w];
    Ared = zeros(complexity*4*n,2*n+2*n);
    Bred = zeros(complexity*4*n,1);
    Xred = zeros(complexity*2*n,2*n);
    Yred = zeros(complexity*2*n,1);

    for ix=1:complexity
        auxADred = toeplitz([Den(support_indices(ix),:) zeros(1,n-1)],[Den(support_indices(ix),1) zeros(1,n-1)]);
        auxANred = toeplitz([Num(support_indices(ix),:) zeros(1,n-1)],[Num(support_indices(ix),1) zeros(1,n-1)]);
        auxBDred = [Den(support_indices(ix),2:n+1) zeros(1,n)]';
        auxBRred = Pref(1,2:2*n+1)';

        Ared((ix-1)*4*n+1:ix*4*n,:) = [[auxADred auxANred; -auxADred -auxANred] [-eye(2*n); -eye(2*n)]];
        Bred((ix-1)*4*n+1:ix*4*n) = [auxBRred-auxBDred; auxBDred-auxBRred];

        Xred((ix-1)*2*n+1:ix*2*n,:) = [auxADred auxANred];
        Yred((ix-1)*2*n+1:ix*2*n) = auxBDred;
    end

    solred = linprog(fred,Ared,Bred);

    if norm(sol-solred)>10^-10
        warning('missing support constraint?')
    end
%%%


figure; hold on;
axis([-1.5 0 -4 4])
axis([-4 0 -4 4])
for ix=1:N
    plot(real(roots(PC(ix,:))),imag(roots(PC(ix,:))),'xb')
end
hold on

MinDumping=0.5;
MaxRealPart=-.7; %main result in the paper
%MaxRealPart=-.6; %uncomment for additional result in the paper
%MaxRealPart=-.5; %uncomment for additional result in the paper
%MaxRealPart=-.4; %uncomment for additional result in the paper

for ix=1:N
    rePC=real(roots(PC(ix,:)));
    imPC=imag(roots(PC(ix,:)));    
    selPC=((rePC>MaxRealPart) | (imPC.^2>rePC.^2*(1/MinDumping^2-1)) );
    plot(rePC(selPC),imPC(selPC),'xr')
end
box on
hold on
%%%

plot(real(roots(Pref)),imag(roots(Pref)),'.k','MarkerSize',20); 
x=[-20:0.001:MaxRealPart];
plot(x,sqrt(x.^2*(1/MinDumping^2-1)),'k:')
plot(x,-sqrt(x.^2*(1/MinDumping^2-1)),'k:')
plot([MaxRealPart,MaxRealPart],[-sqrt(MaxRealPart.^2*(1/MinDumping^2-1)),sqrt(MaxRealPart.^2*(1/MinDumping^2-1))],'k:')
grid on
inappropriateness2=zeros(N,1);
RealViol=zeros(N,1);
DumpViol=zeros(N,1);
for ix=1:N
    real_parts=real(roots(PC(ix,:)));
    imaginary_parts=imag(roots(PC(ix,:)));
    RealViol(ix)=(max(real_parts)>MaxRealPart);
    DumpViol(ix)=sum(imaginary_parts.^2>real_parts.^2*(1/MinDumping^2-1))>0;
    inappropriateness2(ix)=(RealViol(ix)||DumpViol(ix));
end
complexity2=length(find(inappropriateness2>0));

instrumental_complexity= length(union(find(inappropriateness2>0), support_indices))


