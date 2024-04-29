%This code has been used to generate the numerical example in the paper.
%The code requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. 
%Any SDP solver can be used. Here we used SDPT3 freely avaialbe at https://github.com/SQLP/SDPT3
clear all
%%system's data
A=[0.8 0.5; -0.4 1.2];
B=[0;1];
n=2;
m=1;
ubar=5;
%%
%Data collection
T=50;
delta=sqrt(.01);
U=(-1+ (1+1).*rand(1,T));
D=delta*(-1+ (1+1).*rand(2,T));
X(:,1)=rand(2,1);
for i=1:T
X(:,i+1)=A*X(:,i)+B*min(abs(U(:,i)),ubar)*sign(U(:,i))+D(:,i);
end
%%
%Definition of parameters for the LMIs. In this example, l=n,i.e., Htilde=H; 
Q=Tfinal*delta^2*eye(n);
Rcal=[Xd', Ud'];
Term=Xd_plus*(eye(Tfinal)-Rcal*inv(Rcal'*Rcal)*Rcal')*Xd_plus';
Htilde=-1/2*(Term+Term')+Q;
H=Htilde;
T_H=eye(2);
[~,r]=size(null(Htilde));
l=n-r;
E=T_H';
Sigma0=pinv(Rcal)*Xd_plus';
M=Rcal'*Rcal;
L=[zeros(n,m);eye(m)];
%%
%%Solutions to the LMIs
tol=1e-6;
W=sdpvar(n,n,'symmetric');
M_W=sdpvar(n,n,'symmetric');
sigma=sdpvar(1,1,'symmetric');
Z=sdpvar(m,n,'full');
Y=sdpvar(m,n,'full');
Omega= diag(sdpvar(m,1));
Psi_hat=[-W, -Z', zeros(n,l),[W Y']*Sigma0, [W Y'], zeros(n,n+m);

          -Z, -2*Omega, zeros(m,l), Omega*L'*Sigma0, Omega*L', zeros(m,n+m);

          zeros(n,l)', zeros(m,l)', -sigma*inv(H), sigma*E, zeros(l,n+m), zeros(l,n+m);

          ([W Y']*Sigma0)', (Omega*L'*Sigma0)', sigma*E', -W, zeros(n,n+m), zeros(n,n+m);

          [W Y']', L*Omega,zeros(l,n+m)' ,zeros(n,n+m)', -sigma*M,zeros(n+m,n+m);

          zeros(n,n+m)', zeros(m,n+m)', zeros(l,n+m)',zeros(n,n+m)', zeros(n+m,n+m), -sigma*M
];
    
Poly=[W, (Y-Z)'; (Y-Z), ubar^2];

Opti=[M_W, eye(n); eye(n), W];

Constraints=[Psi_hat<=-tol*eye(size(Psi_hat)), Poly>=0,Opti>=0];

options=sdpsettings('solver','sdpt3','verbose',2);

solution=optimize(Constraints,trace(M_W),options);

W=double(W);
Y=double(Y);
K=Y*inv(W);
plotellisa(inv(W),[0;0],'-','b');%plot of the region of attraction
hold on
grid on;
