%%% System Identification Using Subspace Method
%%% Allen Lee

%%% The task is to find A,B,C,D state space matrices that are most
%%% compatible with Input&Output data without the knowledge of the system.
%%% The algorithm supports MIMO systems, unknown initial condition,
%%% and inputs diverse enough for identification.

% clc
% clear

% System Parameters 
n = 8; % system order (if unknown, guess or iterate)
p = 1; % #output
m = 2; % #input

nn = 2*m; % break 1 into 2m responses
data = load('LTI_data_MIMO_noisy.mat');
U = data.W(:,1:m); % input
Yn = data.W(:,m+1:m+p); % noisy output
Y = data.W(:,m+p+1:end); % clean output

[len,n_] = size(U); % U,Y same length

%%% First, create free response to calculate C,A
lag = 2*n+m; % longer, the better
run_sample = 9*n+4*m+5; % more general
Ut = U(1:run_sample*m,:); % for training
Yt = Y(1:run_sample*m,:); % for training
% Free response calculation
%%% Calculation of free responses inpu/output data
H_r = 2*(lag+1); % row of Hankel matrix
H_c = H_r*max(m,p)+n; % just a little more column
U_hk = Hankel_Matrix(Ut',H_r,H_c,m,1,1);% m*H_row,H_col
Y_hk = Hankel_Matrix(Yt',H_r,H_c,p,1,1);% p*H_row, H_col

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
% testing correctness of G,should be zero
% error = sum(sum(abs([W;Uf]*Gc-W0)))/sum(sum(abs([W;Uf])));
% fprintf("Free response calculation errors: %f\n" ...
%     ,error/prod(size(W0)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction Ac, Cc from free response
%%% s.t. Y1*A = Y2
Y1 = Y0(1:end-p,:); % rm last block row
Y2 = Y0(p+1:end,:); % rm first row
Ac = (Y1'*Y1)\Y1'*Y2;
Cc = Y1(1:p,:);
sort(eig(Ac));  % check eig of Ac
%%% Construct Observability Matrix using Ac,Cc and A,C
Oc = Cc;
for i = 1:run_sample-1
    Oc = [Oc;Oc(end-p+1:end,:)*Ac]; % Observability Matrix, recovered
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify B,D,x0
%%% Solve B,D, and X0 using Kronecker Product

Unn = zeros(run_sample*m,nn);% for training
Ynn = zeros(run_sample*p,nn);% for training
for j=1:nn
    shift = (j-1)*run_sample;
    for i=1:run_sample
        Unn((i-1)*m+1:i*m,j) = U(shift+i,:)';
        Ynn((i-1)*p+1:i*p,j) = Y(shift+i,:)';
    end
end

Uk = kron(Unn',eye(run_sample*p));
Ok = kron(eye(nn),Oc);
Yk = zeros(p*run_sample*nn,1);
for i=1:size(Ynn,2)
    Yk((i-1)*p*run_sample+1:i*p*run_sample,1) = Ynn(:,i);
end
Uk1 = Uk;
%%% Upper zeros
for k = 0:m-1
    for i=1:run_sample-1
        shift = i*p*run_sample*m+k*p*run_sample;
        Uk1(:,shift+1:shift+i*p) = zeros(p*run_sample*nn,i*p);
    end
end

%%% Uk reduction
Ukr = zeros(p*run_sample*nn,p*run_sample*m);
dummy1 = zeros(p*run_sample*nn,p);
for i=1:m*run_sample
    j = (i-1)*p+1;
    dummy1 = zeros(p*run_sample*nn,p);
    while(j<=p*run_sample^2*m)
        dummy1 = dummy1 + Uk1(:,j:j+p-1);
        j = j+p*m*run_sample+p;
    end
    Ukr(:,(i-1)*p+1:i*p) = dummy1;
end

M_ = [Ok Ukr];
if(rank(M_)<size(M_,2))
    warning("M_ is not full column rank!")
    fprintf("The rank of M_ is %d\n",rank(M_))
end
Sol = (M_'*M_)\M_'*Yk; % left inverse
x0c = zeros(n,nn);
for i=1:nn
    x0c(:,i) = Sol((i-1)*n+1:i*n);
end

Gk = Sol(nn*n+1:end);
% sum(abs(Yk - Ok*Sol(1:n*nn) - Ukr*Gk)) % testing, should be zero

Gc = zeros(p*run_sample,m);
for i = 1:m
    Gc(:,i) = Gk((i-1)*p*run_sample+1:i*p*run_sample,1);
end


Dc = Gc(1:p,1:m); % D term
CAB = Gc(p+1:end,1:m);% first column, but no D
%%% formula: CAB = Or*B = Oc*Bc;
Bc = (Oc(1:end-p,:)'*Oc(1:end-p,:))\Oc(1:end-p,:)'*CAB;

%%% responses reconstrution
%%% Verification
Yc = zeros(len,p); % Output collection, nn set
x1 = x0c(:,1); % first IC is the true IC

for i = 1:len
    u = U(i,:)';
    Yc(i,:) = Cc*x1+Dc*u;
    x1 = Ac*x1+Bc*u;
end



%% Data points comparison
for i=1:p
    figure
    hold on
    plot(Y(:,i),'r.')
    plot(Yc(:,i),'b.') % exploded when sys unstable

    hold off
    legend('real','recovered')
end
%% Frequency Response Comparison
sys_rec = ss(Ac,Bc,Cc,Dc,-1);
[num_rec,den_rec] = tfdata(sys_rec,'v');
sys_rec = tf(num_rec,den_rec,-1);
for i =1:m
    figure
    hold on
    for j = 1:p
        bode(sys_rec(j,i),sys(j,i))% sys is from LIT_Data_generator
        legend('true','recovered')
    end

end
hold off