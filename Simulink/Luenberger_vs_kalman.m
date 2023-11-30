load data1.mat

T = out.data.Time;
Id_s = out.data.Data(:, 1);
Iq_s = out.data.Data(:, 2);
speed = out.data.Data(:, 3);
theta = out.data.Data(:, 4);
Vd_s = out.data.Data(:, 5);
Vq_s = out.data.Data(:, 6);

Tload = 0;

th = 0;
w = 0;
w_luen = 0;
w_kalman = 0;
w_AEKF = 0;

X = [0; 0; 0; 0];
X_kalman = X;
X_luen = X;
X_AEKF = X;

Ts = 5e-5;

Rs = 3.4;
L = 0.01215;
lambda = 0.25;
J = 2.9e-4;
p = 4;

k = length(T/30);

R = [1 0; 0 1];
P_kalman = eye(4, 4);

P_AEKF = eye(4, 4);

Q_kalman = [0.4 0 0 0 ;
     0 0.004 0 0 ;
     0 0 200 0 ;
     0 0 0 2 ];

Q_AEKF = [0.4 0 0 0 ;
          0 0.004 0 0 ;
          0 0 200 0 ;
          0 0 0 2 ];

a1 = -Rs/L;
a2 = w;
b1 = -w;
b2 = -Rs/L;
b3 = -lambda/L;
c1 = (3/2)*lambda/J;

Fx = [a1 a2 0 0;
     b1 b2 b3 0;
     0 c1 0 0;
     0 0 1 0 ];

A = eye(4, 4) + Fx*Ts;
%A = Fx;

B = [1/L 0;
     0 1/L;
     0 0;
     0 0];

B = B*Ts;

C = [1 0 0 0;
     0 1 0 0];

theta_est = [];
w_est = [];
w_est_kalman = [];
w_est_luen = [];
w_est_AEKF = [];

L_gain = [0.01 0;
          0 0.01;
          0.05 0.0;
          0.0 0.05];


tic
for i = 1:1:k
    id           = Id_s(i);
    iq           = Iq_s(i);
    vd           = Vd_s(i);
    vq           = Vq_s(i);

    w            = X(3);
    w_kalman     = X_kalman(3);
    w_luen       = X_luen(3);
    w_AEKF       = X_AEKF(3);

    th           = X(4);
    theta_est = [theta_est; th];

    w_est = [w_est; w*30/pi];
    w_est_kalman = [w_est_kalman; w_kalman*30/pi];
    w_est_luen = [w_est_luen; w_luen*30/pi];
    w_est_AEKF = [w_est_AEKF; w_AEKF*30/pi];

    a1 = -Rs/L;
    a2 = w;
    b1 = -w;
    b2 = -Rs/L;
    b3 = -lambda/L;
    c1 = (3/2)*lambda/J;
    
    Fx = [a1 a2 0 0;
         b1 b2 b3 0;
         0 c1 0 0;
         0 0 1 0];

    Fx1 = [a1 w_kalman 0 0;
         -w_kalman b2 b3 0;
         0 c1 0 0;
         0 0 1 0];

    Fx2 = [a1 w_luen 0 0;
         -w_luen b2 b3 0;
         0 c1 0 0;
         0 0 1 0];

    Fx3 = [a1 w_AEKF 0 0;
         -w_AEKF b2 b3 0;
         0 c1 0 0;
         0 0 1 0];
    
    A = eye(4, 4) + Fx*Ts;
    A_kalman = eye(4, 4) + Fx1*Ts;
    A_luen = eye(4, 4) + Fx2*Ts;
    A_AEKF = eye(4, 4) + Fx3*Ts;
    % A = Fx;

    U = [vd; vq];
    X   = (A * X) + (B * U);
    X_kalman   = (A_kalman * X) + (B * U);
    X_luen   = (A_luen * X) + (B * U);
    X_AEKF = (A_AEKF * X) + (B * U);

    y = C*X;
    y_kalman = C*X_kalman;
    y_luen = C*X_luen;
    y_AEKF = C*X_AEKF;

    error   = [id; iq] - y;
    error1   = [id; iq] - y_kalman;
    error2   = [id; iq] - y_luen;
    error3   = [id; iq] - y_AEKF;


    P_kalman   = (A_kalman * P_kalman * A_kalman') + Q_kalman;
    KalmanGain_x_EKF = (P_kalman) * (C') * (inv((C * P_kalman * C') + (R)));
    X_kalman            = X_kalman + (KalmanGain_x_EKF * error1);
    P_kalman            = (eye(4,4) - (KalmanGain_x_EKF * C)) * P_kalman;

    P_AEKF   = (A_AEKF * P_AEKF * A_AEKF') + Q_AEKF;
    KalmanGain_x_AEKF = (P_AEKF) * (C') * (inv((C * P_AEKF * C') + (R)));
    X_AEKF            = X_AEKF + (KalmanGain_x_AEKF * error3);
    P_AEKF            = (eye(4,4) - (KalmanGain_x_AEKF * C)) * P_AEKF;
    Q_AEKF          = KalmanGain_x_AEKF*(error3*error3')*transpose(KalmanGain_x_AEKF);

    X_luen = X_luen + (L_gain*error2);

end
toc

plot(speed(1:k)*30/pi,'k');
hold on
plot(w_est, 'r');
plot(w_est_luen, 'g');
plot(w_est_kalman, 'b');
plot(w_est_AEKF, 'y');
hold off
RMSE = rmse(mod(theta_est, 2*pi), theta(1:k));


