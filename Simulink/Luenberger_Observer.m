T = out.data.Time;
Id_s = out.data.Data(:, 1);
Iq_s = out.data.Data(:, 2);
speed = out.data.Data(:, 3);
theta = out.data.Data(:, 4);
Valpha_s = out.data.Data(:, 5);
Vbeta_s = out.data.Data(:, 6);
I_alpha = out.data.Data(:, 7);
I_beta = out.data.Data(:, 8);

th = 0;
w = 0;
X = [0; 0; 0; 0]; 
Ts = 5e-5;

Rs = 3.4;
L = 0.01215;
lambda = 0.25;
J = 2.9e-4;
p = 4;

k = length(T)/10;

a1 = -Rs/L;
a2 = w;
b1 = -w;
b2 = -Rs/L;
b3 = -lambda/L;

Fx = [a1 0 1/L 0;
     0 b2 0 1/L;
     0 0 0 b1;
     0 0 a2 0];

A = eye(4, 4) + Fx*Ts;

B = [1/L 0;
     0 1/L;
     0 0;
     0 0];

B = B*Ts;

C = [1 0 0 0;
     0 1 0 0];

theta_est = [];
w_est = [];

L_gain = [1 0;
          1E-3 0;
          1 1E-5;
          1 1E-3];

Q = 1E15 * eye(4);
R = 1E15;

tic
for i = 1:1:k
    ia           = I_alpha(i);
    ib           = I_beta(i);
    va           = Valpha_s(i);
    vb           = Vbeta_s(i);

    th_new       = atan(X(4)/X(3));
    th_new
    if isnan(th_new)
        th_new = 0;
    end
 
    w            = (th_new-th)/(Ts);
    %w
    th           = th_new;
    theta_est = [theta_est; th];
    w_est = [w_est; w*30/pi];

    a1 = -Rs/L;
    a2 = w;
    b1 = -w;
    b2 = -Rs/L;
    b3 = -lambda/L;
    
    Fx = [a1 0 1/L 0;
         0 b2 0 1/L;
         0 0 0 b1;
         0 0 a2 0];
    
    A = eye(4, 4) + Fx*Ts;

    U = [va; vb];
    X_hat   = (A * X) + (B * U);

    y = C*X;
    error   = [ia; ib] - y;

    co = ctrb(A', C');
    ob = obsv(A, C);
    
    %p = [-1.0, -1.0, -0.995, -0.99];
    %temp = place(A', C', p);
    K = dlqr(A', C', Q, R);
    K

    L_gain = K';
    %e = eig(A-L_gain*C);
    %e

    X = X_hat + L_gain*error;
end
toc

plot(theta,'k');
hold on
plot(theta_est,'r');
hold off
RMSE = rmse(mod(theta_est, 2*pi), theta(1:k));


