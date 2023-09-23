%%% system identification for Full State Measurement %%%
%%% Allen Lee

%%% Assuming full access to state and output
%%% random initial condition, w or w/o noise
%%% Full identification to A,B,C,D
clc
clear

%%% Generate SS sys
n = 4; %system order
p = 1; %#output
m = 1; %#input, doesn't matter here
noisy = false;
%%%%%%%%%%%%%%%
sys = drss(n,p,m);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

sys_real = ss(sys.A,sys.B,sys.C,sys.D,-1);
sys_real = tf(sys_real);
run_sample = (2*n+p+m); % estimated maximum lag, long enough
Y = zeros(p,run_sample);
U = zeros(m,run_sample);
X = zeros(n,run_sample); % collect states
x = -5+10*rand(n,1); % col:# ICs % Manual
u = -1+2*rand(m,1);
for i = 1:run_sample
    X(:,i) = x;
    u = -1+2*rand(m,1);
    Y(:,i) = C*x+D*u;
    U(:,i) = u;
    x = A*x+B*u; 
end


if(noisy)
    % If add measurement noise
    noise = -0.05+0.1*rand(size(Y)); % noise:+-0.05 % Manual
    Y = Y + noise; % Corrupted Measurement
    noise = -0.05+0.1*rand(size(X)); % noise:+-0.05 % Manual
    X = X + noise; % Corrupted States
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction Ac,Bc, Cc,Dc
%%% s.t. [X(shifted);Y] = [A B;C D]*[X;U]
%%% denoted as NN = SS*MM
% Find right inverse of MM:
NN = [X(:,2:end);Y(:,1:end-1)];
MM = [X(:,1:end-1);U(:,1:end-1)];
SS = NN*MM'/(MM*MM'); % [A B;C D] matrix, right_inv of MM

Ac = SS(1:n,1:n);
Bc = SS(1:n,n+1:end);
Cc = SS(n+1:end,1:n);
Dc = SS(n+1:end,n+1,end);

[sort(eig(A)) sort(eig(Ac))]; % examine, should be the same

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bode Plots for comparison
% should be identical when no noise
sys_rcv = ss(Ac,Bc,Cc,Dc,-1);
sys_rcv = tf(sys_rcv);
figure 
bode(sys_real,sys_rcv)
legend('real','recovered')








