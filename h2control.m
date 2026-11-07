%REQUIREMENT: CVX - Software for Disciplined Convex Programming 
%https://cvxr.com/cvx/doc/solver.html
%Tested with Version 2.2, Build 1148, MATLAB 9.8 (R2020a)

rng(04032025) %reproducible run

nx = 4;
nu = 2;
Lp1 = -2.93; Lp2= -1;
Lb1 = -4.75; Lb2 = -73.14;
Lr1 = 0.78; Lr2 = 3.18;
Nbp1 = 0; Nbp2 = 0.1;
Np1 = -0.042; Np2 = 0.086;
Nb1 = 2.59; Nb2 = 8.94;
Nr1 = -0.39; Nr2 = -0.29; 

A = [0 1 0 0; 0 Lp1 Lb1 Lr1; 0.086 0 -0.11 -1; 0.086*Nbp1 Np1 Nb1-0.11*Nbp1 Nr1-Nbp1];
B = [0 0; 0 -3.91; 0.035 0; -2.53 0.31];


x0 = ones(nx,1);
R = eye(nu);
S = 0.01*eye(nx);
Z = zeros(nx,nu);
aux = sqrtm([S Z; Z' R]);
C = aux(:,1:nx);
D = aux(:,nx+1:nx+nu);

%%% number of optimization variables
d = nx*(nx+1)/2 + nu*nx;

cvx_clear;
cvx_solver sdpt3
oflag = cvx_quiet(true);
tol = 1e-6;

N=2000;

%%% optimization outputs
P = zeros(nx,nx); 
K = zeros(nu,nx);
J = zeros(1);
 

%%% randomization over uncertainty for design
dLp = (Lp2-Lp1)*rand(1,N)+Lp1;
dLb = (Lb2-Lb1)*rand(1,N)+Lb1;
dLr = (Lr2-Lr1)*rand(1,N)+Lr1;
dgV = 0.086*ones(1,N);
dYb = -0.11*ones(1,N);
dNbp = (Nbp2-Nbp1)*rand(1,N)+Nbp1;
dNp = (Np2-Np1)*rand(1,N)+Np1;
dNb = (Nb2-Nb1)*rand(1,N)+Nb1;
dNr = (Nr2-Nr1)*rand(1,N)+Nr1;

par = [dLp; dLb; dLr; dNbp; dNp; dNb; dNr]';

dA = zeros(nx,nx,N);
dA(1,2,:) = ones(1,N);
dA(3,4,:) = -ones(1,N);
dA(2,2,:) = dLp;
dA(2,3,:) = dLb;
dA(2,4,:) = dLr;
dA(3,1,:) = dgV;
dA(3,3,:) = dYb;
dA(4,1,:) = dNbp.*dgV;
dA(4,2,:) = dNp;
dA(4,3,:) = dNb+dNbp.*dYb;
dA(4,4,:) = dNr-dNbp;


cvx_begin sdp

variable Y(nx,nx) symmetric;
variable X(nu,nx);

minimize matrix_frac(x0,Y);
 subject to

    Y >= 0;

    for h=1:N 
       [ -dA(:,:,h)*Y-Y*dA(:,:,h)'-B*X-X'*B' (C*Y+D*X)'; C*Y+D*X eye(nu+nx) ] >= 0;
    end

cvx_end

P(:,:) = inv(Y);
K(:,:) = X*inv(Y);
J(1) = x0'*inv(Y)*x0;
J(1)

cvx_quiet(oflag);

smalleig=zeros(1,N);
for h=1:N
    smalleig(h)=min(eig([ -dA(:,:,h)*Y-Y*dA(:,:,h)'-B*X-X'*B' (C*Y+D*X)'; C*Y+D*X eye(nu+nx) ]));
end
indices_of_support_elements=find(smalleig<tol);
complexity=length(indices_of_support_elements)
dAsupport=dA(:,:,indices_of_support_elements);

gamma=0.8;
rob_smalleig=zeros(1,N);
for h=1:N
    rob_smalleig(h)=min(eig([ gamma*(-dA(:,:,h)*Y-Y*dA(:,:,h)'-B*X-X'*B') (C*Y+D*X)'; C*Y+D*X eye(nu+nx) ]));
end
bound_instrumental_complexity=length(find(rob_smalleig<tol))
