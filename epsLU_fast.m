function [outL, outU] = epsLU_fast(bet,N,varargin)

if nargin >= 4
    r = varargin{2};    
else
    r = 2;
end

if nargin >= 3
    range = varargin{1};    
else
    range = [0:1:N];
end

outL = zeros(length(range),1);
outU = zeros(length(range),1);

for k = range
    
    t1 = k/N;
    t2 = 1;
    while t2-t1 > 1e-10         
        
        t = (t1+t2)/2;
        val1 = (bet/2-bet/2/(r-1))*betainc(t,k+1,N-k) + bet/2/(r-1)*betainc(t,k+1,r*N+1-k);
        val2 = (1+bet/2/(r-1))*t*N*(betainc(1-t,N-k,k+1)-betainc(1-t,N-k+1,k));
        val = val1-val2;
        
        if val > 0
            t2=t;
        else
            t1=t;
        end
        
    end
    outU(k-range(1)+1) = t2;
    
    t1 = 0;
    t2 = k/N;
    while t2-t1 > 1e-10         
        
        t = (t1+t2)/2;
        val1 = (bet/2-bet/2/(r-1))*betainc(t,k+1,N-k) + bet/2/(r-1)*betainc(t,k+1,r*N+1-k);
        val2 = (1+bet/2/(r-1))*t*N*(betainc(1-t,N-k,k+1)-betainc(1-t,N-k+1,k));
        val = val1-val2;
        
        if val > 0
            t1=t;
        else
            t2=t;
        end
        
    end
    outL(k-range(1)+1) = t1;
end