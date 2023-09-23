%%% System Identification from Extended Observability Matrix %%%
%%% Allen Lee

%%% Free responses, random initial condition, w or w/o noise
clc
clear

%%% Generate SS sys
n = 4; %system order
p = 2; %#output
m = 1; %#input, doesn't matter here
noisy = true;
%%%%%%%%%%%%%%%
sys = drss(n,p,m);
A = sys.A;
B = sys.B;
C = sys.C;
D = zeros(p,m);

sys = ss(sys.A,sys.B,sys.C,zeros(p,m),-1);
sys = tf(sys);
run_sample = (2*n+p+m); % estimated maximum lag, long enough
x0_num = n;
Y = zeros(p*run_sample,x0_num);
x0 = -5+10*rand(n,n); % col:# ICs % Manual
x = zeros(n,1); % preallocated memory
for j = 1:x0_num
    x = x0(:,j); % +-5 , IC
    for i = 1:run_sample
        Y((i-1)*p+1:p*i,j) = C*x;
        x = A*x; % free responses
    end
end

if(noisy)
    % If add measurement noise
    noise = -0.05+0.1*rand(size(Y)); % noise:+-0.05 % Manual
    Y = Y + noise; % Corrupted Measurement
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction Ac, Cc
%%% s.t. Y1*A = Y2
Y1 = Y(1:end-p,:); % rm last block row
Y2 = Y(p+1:end,:); % rm first row

Ac = (Y1'*Y1)\Y1'*Y2;
Cc = Y(1:p,:);

% Construct new ICs, x0c
r = rank(Ac); % sys state dim
Obs = zeros(p*r,n); % Free responses matrix, Observability Matrix
for i = 0:r-1
    Obs(i*p+1:(i+1)*p,:) = Cc*Ac^(i);
end
% Obs*x0c = Y ,And X0c must be Identity matrix
x0c = (Obs'*Obs)\Obs'*Y(1:p*n,:);% constructed ICs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% responses reconstrution
%%% Verification

Yc = zeros(p*run_sample,x0_num);
for j = 1:x0_num
    x = x0c(:,j); % +-5 , IC
    for i = 1:run_sample
        Yc((i-1)*p+1:p*i,j) = Cc*x;
        x = Ac*x; % free response
    end
end



for i = 1:size(x0,2) % # ICs
    figure
    hold on
    % two plots should be identical when no noise
    plot(Y(:,i),'r.')
    plot(Yc(:,i),'b.')
    hold off
    legend('real','recovered')
end






