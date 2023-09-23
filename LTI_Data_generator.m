%%% Data generator of LTI systems for system identification
%%% Allen Lee
clc
clear

%%% Using 1 response, with different channel of inputs & outputs
%%% Input&Output are separated by columns, only stable systems in this case.
%%% Unstable systems are unadvisable to be identified using time-domain
%%% data, try frequency domain identification when systems are unstable.

%%% System definition
n = 8;  % system states
p = 1;  % output dimension
m = 2;  % input dimension
run_sample = 2048; % Data Length

%%% System generation
sys = drss(n,p,m); % only stable discrete LTI
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

%%% Process data generation
Y = zeros(run_sample,p); % output collection
x0 = -1+2*rand(n,1); % random IC
u = -1+2*rand(m,1); % control preallocation
U = zeros(run_sample,m); % control input collection
x1 = x0(:,1); % single states

for i = 1:run_sample
    u = -1+2*rand(m,1); % random input
    Y(i,:) = (C*x1+D*u)';
    x1 = A*x1+B*u;
    U(i,:) = u';
end

noise = -0.025+0.05*rand(size(Y));% Manual
Yn = Y + noise; % Corrupted Measurement

W = [U Yn Y];
size(W);

% save('LTI_data_MIMO.mat','W');
save('LTI_data_MIMO_noisy.mat','W');

sort(eig(A));

%% Generate system for comparison in SysID
A
B
C
D
x0
sys = ss(A,B,C,D,-1);
figure
plot(Y,'.')
figure
plot(U,'.')
