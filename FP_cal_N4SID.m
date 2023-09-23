%%% Caculate free responses of a system, using N4SID
%%% Allen Lee
clear
clc
%%% Given random input and unknown IC to generate arbitrary trajectory,
%%% Calculate Free response to Find A,C matrices


%%% Generate SS sys
n = 4; % system order
p = 3; % #output
m = 1; % #input % works only when m=1
noisy = false; % if data are corrupted
%%%%%%%%%%%%%%%
%%% Generate response (in&output data)
lag = 2*n+m; % long enough, at lease 2*n
run_sample = 9*n+4*m+5; % more general
beta = 5;
sys = drss(n,p,m);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;


x0_num = 1; %fixed 1!
Y = zeros(p,run_sample*beta);
U = zeros(m,run_sample*beta);
x0 = -5+10*rand(n,1); % col:# ICs % Manual
x1 = zeros(n,1); % preallocated memory
u = -1+2*rand(m,1);
for j = 1:x0_num
    x1 = x0(:,j); 
    for i = 1:run_sample*beta
        u = -1+2*rand(m,1);
        Y(:,i) = C*x1+D*u;
        x1 = A*x1+B*u; % random response
        U(:,i) = u;
    end
end
if(noisy)
    % If add measurement noise
    noise = -0.05+0.1*rand(size(Y)); % noise:+-0.05 % Manual
    Y = Y + noise; % Corrupted Measurement
end

%%% Calculation of free responses inpu/output data
H_r = 2*(lag+1);
H_c = H_r*max(m,p)+n; % just a little more column
U_hk = Hankel_Matrix(U,H_r,H_c,m,1,1);% m*H_row,H_col
Y_hk = Hankel_Matrix(Y,H_r,H_c,p,1,1);% p*H_row, H_col

Up = U_hk(1:m*(lag+1),:);
Uf = U_hk(m*(lag+1)+1:end,:); % Make sure Uf is long enough
Yp = Y_hk(1:p*(lag+1),:);
Yf = Y_hk(p*(lag+1)+1:end,:);
%%% free response: Y0 = Yf*G, by N4SID algorithm
W = [Up;Yp];
W0 = [W;zeros(size(Uf))];
NS = [W*W' W*Uf';Uf*W' Uf*Uf'];
Gc = [W' Uf']*pinv(NS)*W0;
Y0 = Yf*Gc;
Y0 = Y0(:,1:n); % n free responses are everything we need
% sum(sum(abs([W;Uf]*Gc-[W;zeros(size(Uf))]))) % testing, should be zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction Ac, Cc from free response
%%% s.t. Y1*A = Y2
Y1 = Y0(1:end-p,:); % rm last block row
Y2 = Y0(p+1:end,:); % rm first row

Ac = (Y1'*Y1)\Y1'*Y2;
Cc = Y1(1:p,:);

%%% Comparison 
% the eigenvalues of A and Ac should be the same
[sort(eig(A)) sort(eig(Ac))]

