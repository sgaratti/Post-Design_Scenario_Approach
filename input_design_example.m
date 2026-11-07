%REQUIREMENT: CPLEXLP
%Tested with CPLEX Optimization Studio V12.8, MATLAB R2020a


N = 1000;
T = 5;
u_lim = 10;

Anom = [0.8 -1; 0 -0.9];
Bnom = [0; 0.25];

assert(size(Anom,1)==size(Anom,2));
assert(size(Anom,1)==size(Bnom,1));
n = size(Anom,1);
p = size(Bnom,2);

sig = 0.05;
tag = 0.97;
c1 = 2;
c2 = 2;


allA = repmat(Anom,N,1)+(sig+2*sig*(rand(n*N,n)>tag)).*randn(n*N,n);
allB = repmat(Bnom,N,1);

x0 = [1; 1];


rhoL = [0.05 1];
L = length(rhoL);

U = zeros(p*T,L);
H = zeros(1,L);
XI = zeros(N,L);
IX = cell(1,L);

for i = 1:L
    rho = rhoL(i)
    [u, h, xi] = solve_mnid(allA,allB,x0,T,rho,u_lim);
    U(:,i) = u; %designed input
    H(i) = h; %h*_N
    XI(:,i) = xi;
end


tot = 100*N;
testA = repmat(Anom,tot,1)+(sig+2*sig*(rand(n*tot,n)>tag)).*randn(n*tot,n);
testB = repmat(Bnom,tot,1);



[epsL, epsU] = epsLU_fast(1e-7,N);

gridN = 200;
lvl_rel = [1:gridN]'/gridN*0.3 -0.3;
IVa = -0.1;
IVb = 0.2;
IV = [IVa:0.001:IVb];
for i = 1:L
    u = U(:,i);
    h = H(:,i);
    xi = XI(:,i);
    nv = sum(xi>0);
    
    lvl = lvl_rel+h;
        
    [Z, Z_norm] = test_mnid(allA,allB,x0,T,u);
    [Zn Zn_ix] = sort(Z_norm);
    nsclvl = zeros(gridN,1);
    for k = 1:gridN
        if lvl(k)==h
            kk =  sum(Zn >= lvl(k));
        else
            kk = sum(Zn > lvl(k));
        end
        nsclvl(k) = kk;
    end
    
    [W, W_norm] = test_mnid(testA,testB,x0,T,u);
    prIV = zeros(length(IV),1);
    for k = 1:length(IV)
        prIV(k) = sum(W_norm<=IV(k))/tot;
    end
    
    figure; hold on; grid on;
    axis([IVa IVb 0 1]);
    yticks([0:0.5:1]);
    xlabel('||\eta(T)||_\infty')
    ylabel('CDF')
    stairs([IVa; lvl; IVb],[1-epsL(nsclvl+1).*(lvl<=h); 1; 1],'k')
    

    stairs([IVa; lvl; IVb],[0; 1-epsU(nsclvl+1); 1-epsU(nsclvl(end)+1)],'k') 

    stairs(IV,prIV,'--k')    
    title(['rho=',num2str(rhoL(i))]);
end



    


    


