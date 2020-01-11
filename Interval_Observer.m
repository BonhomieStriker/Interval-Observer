clear all
close all
clc
%% Iteration Parameters
iter = 200;
n = 1:iter;
t = n * 0.01;

%% System Parameters
A = [0.5  1
     0   -0.5];
B = [1 0
     0 1];
E = [1  1/2]';
C = [1 -1
     0  1];
F = E;
L = [0.25 0.25
     0   -0.75];    

%% Interval Observer Parameters
[S, R] = eig(A - L *C);


%% Kalman Filter Parameters


%% Initialize State Space Model
%Input
u1 = 2 + sin(2*pi*t);
u2 = 2 + cos(2*pi*t);
u1_l = min(u1);
u1_u = max(u1);
u2_l = min(u2);
u2_u = max(u2);
u = [u1;u2];
u_l = [u1_l;u2_l];
u_u = [u1_u;u2_u];
%System
x = zeros(2,iter);
y = zeros(2,iter);
%Luenberger Observer
xo = zeros(2,iter);
yo = zeros(2,iter);
%Kalman Filter Observer


%Interval Observer
z = zeros(2,iter);
x_l = x;
x_u = x;
z_l = z;
z_u = z;
y_l = y;
y_u = y;
%Observation Error
eo = zeros(2,iter);
e_l = eo;
e_u = eo;

%White Noise Generation
w = randn(1,iter);
w = w/std(w);
w = w - mean(w);        
% v = 0.3 + 0.01 * w;
w = 0.5 + 0.1 * w;
figure(1)
plot(t,w,'linewidth',2)
grid on
xlabel('t / s')
ylabel('Value of white noise')
legend('white noise v')
%% Iterative Observation
for k = 1:iter    

    
    %System
    y(:,k) = C * x(:,k) + F * w(k);
    x(:,k+1) = A * x(:,k) + B * u(:,k);
    
    %Luenberger Observer
    yo(:,k) = C * xo(:,k);
    xo(:,k+1) = (A - L * C) * xo(:,k) + B * u(:,k) + L * y(:,k);
    
    %Kalman Filter

    
    %Interval Observer
%     z_l(:,k+1) = R * z_l(:,k) + S * B * u(:,k) + (S * E) * min(w) - (S * L * F) * min(v);
%     z_u(:,k+1) = R * z_u(:,k) + S * B * u(:,k) + (S * E) * max(w) - (S * L * F) * max(v);
%     x_l(:,k+1) = abs(inv(S)) * z_l(:,k) - (abs(inv(S)) - S) * z_u(:,k);
%     x_u(:,k+1) = abs(inv(S)) * z_u(:,k) - (abs(inv(S)) - S) * z_l(:,k);
%     y_l(:,k) = abs(C) * x_l(:,k) - (abs(C) - C) * x_u(:,k) + abs(-F) * max(v) - (abs(-F)+F) * max(v);
%     y_u(:,k) = abs(C) * x_u(:,k) - (abs(C) - C) * x_l(:,k) + abs(-F) * max(v) - (abs(-F)+F) * min(v);
    x_l(:,k+1) = A * x_l(:,k) + u(:,k) + L * (y(:,k) - C * x_l(:,k)) - (abs(L) - L)*[1;1] * max(w);
    x_u(:,k+1) = A * x_u(:,k) + u(:,k) + L * (y(:,k) - C * x_u(:,k)) + (abs(L) - L)*[1;1] * max(w);
    y_l(:,k) = abs(C) * x_l(:,k) - (abs(C) - C) * x_u(:,k) + abs(-F) * min(w) - (abs(-F)+F) * max(w);
    y_u(:,k) = abs(C) * x_u(:,k) - (abs(C) - C) * x_l(:,k) + abs(-F) * max(w) - (abs(-F)+F) * min(w);
end
%Residual Error
eo = y - yo;
e_l = y - y_u;
e_u = y - y_l;

%% Plotting
figure(2)
subplot(211)
plot(t,y(1,:),t,y_l(1,:),'-.b',t,y_u(1,:),'--r','linewidth', 2)
grid on
xlabel('t / s')
ylabel('x_1')
legend('x_1', 'lower bound', 'upperbound')
subplot(212)
plot(t,y(2,:),t,y_l(2,:),'-.b',t,y_u(2,:),'--r','linewidth', 2)
grid on
xlabel('t / s')
ylabel('x_2')
legend('x_2', 'lower bound', 'upperbound')