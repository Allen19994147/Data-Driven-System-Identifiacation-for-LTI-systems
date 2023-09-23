%%% System Identification using Hankel Matrix%%%
%%% It works for single input system with constraints:
%%% Impulse response, zero IC, D=0

clc
clear

%%% First, generate a discrete random model
n = 5; % system order
p = 2; % #output
m = 1; % #input     % must be one here!!
sys = drss(n,p,m);  % random discrete SS generation
A = sys.A;
B = sys.B;
C = sys.C;
sys = ss(A,B,C,zeros(p,m),-1); % Enforce D=0

%%% Generate data
x0 = zeros(n,1); % enforced zero IC
run_sample = 100; % Maunally adjust!
x = x0; % system states
Y0 = zeros(p,run_sample*m);
u = zeros(m,1); % Impulse reponse
for i = 1:run_sample
    
    if(i==1)
        u = ones(m,1);
    else
        u = zeros(m,1);
    end
    x = A*x+B*u;
    Y0(:,(i-1)*m+1:i*m) = C*x; % No D term
end
%%% Data analysis, State Space Reconstruction

H_row = n+p; % Manually Adjusted
H_col = n+p;

H0 = Hankel_Matrix(Y0,H_row,H_col,p,m,1); % starting from 
[U,Sin_value,V] = svd(H0); %H0 = UDV';
Vt = V';
Ns = rank(Sin_value);% dimension of states are known here bc no noise
O_n = U(:,1:Ns);
C_n = Sin_value(1:Ns,1:Ns)*Vt(1:Ns,:);
O_n_linv = (O_n'*O_n)\O_n'; % left inverse
C_n_rinv = C_n'/(C_n*C_n'); % right inverse
H1 = Hankel_Matrix(Y0,H_row,H_col,p,m,2); % one shift

% Constructed A,B,C. Denoted Ac,Bc,Cc
Ac = O_n_linv*H1*C_n_rinv;
Bc = C_n(:,1:m);
Cc = O_n(1:p,:);
Dc = zeros(p,m); % Enforced to be zero here

sys_c = ss(Ac,Bc,Cc,Dc,-1);
bode(sys,sys_c) % Compare by bode plots

%%
%%% Impulse Data verification %%%
x0 = zeros(n,1); % zero IC
x = x0; % system states
Yc = zeros(p,run_sample);
u = zeros(m,1); % Impulse reponse
for i = 1:run_sample   
    if(i==1)
        u = ones(m,1);
    else
        u = zeros(m,1);
    end
    x = Ac*x+Bc*u;
    Yc(:,i) = Cc*x;
end
%%% plot, assuming SISO
t = 1; % examine t-th channel
figure
hold on
plot(Yc(t,:))
plot(Y0(t,:))
legend('Yc',"Y0")
