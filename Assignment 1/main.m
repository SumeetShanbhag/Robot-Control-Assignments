clear; 
clc; 
close all;

%% Defining / Creating Symbols
syms l1 l2 r1 r2 h1 h2 theta_1(t) theta_2(t) theta_ddot1 theta_ddot2 T1 T2 g 'real'
syms m1 m2 v1 v2 w1 w2 I1 I2 'real'

sympref('AbbreviateOutput',false);
sympref('MatrixWithSquareBrackets',true);
sympref('PolynomialDisplayStyle','ascend');

%% Declaring Position of Link 1 and Link 2
x1 = r1*sin(theta_1);
y1 = r1*cos(theta_1);
x2 = l1*sin(theta_1) + r2*sin(theta_1 + theta_2);
y2 = l1*cos(theta_1) + r2*cos(theta_1 + theta_2);

%% Differentiating Position values with respect to time 
x1_dot = diff(x1,t);
x2_dot = diff(x2,t);
y1_dot = diff(y1,t);
y2_dot = diff(y2,t);

%% Defining Generalised Co-ordinates
q = sym('q', [2,1]);
q(1) = theta_1;
q(2) = theta_2;

%% Defining Generalised Inputs
u = sym('u', [2,1]);
u(1) = T1;
u(2) = T2;

%% Differentiating Generalised Coordinates
q_dot = diff(q,t);

%% Kinetic Energy Equations
KE = (m1*(v1.^2))/2 + (I1*w1.^2)/2 +  (m2*(v2.^2))/2 + (I2*w2.^2)/2;
KE = subs(KE,[v1,v2,w1,w2],[sqrt(x1_dot^2 + y1_dot^2), sqrt(x2_dot^2 + y2_dot^2), diff(theta_1,t), diff(theta_1 + theta_2,t)]);
KE = simplify(KE);

%% Potential Energy Equations
PE = (m1*g*h1) + (m2*g*h2);
PE = subs(PE,[h1,h2],[y1,y2]);

%% Lagrangian 
L = KE - PE;
L = simplify(L); 

%% Finding the terms for forming equations of motion
DL_dq = jacobian(L,q);
DL_Ddq = jacobian(L,q_dot);
DL_dtDdq = diff(DL_Ddq, t);

%% Equations of motion for two links
eq1 = DL_dtDdq(1) - DL_dq(1) - u(1);
eq2 = DL_dtDdq(2) - DL_dq(2) - u(2);
display(eq1)
display(eq2)
eq1 = subs(eq1,[diff(theta_1,t,2), diff(theta_2,t,2)],[theta_ddot1, theta_ddot2]);
eq2 = subs(eq2,[diff(theta_1,t,2), diff(theta_2,t,2)],[theta_ddot1, theta_ddot2]);
sol = solve([eq1 == 0, eq2 == 0],[theta_ddot1, theta_ddot2],'ReturnConditions', true); 

%% Display these soluttions
display(sol.theta_ddot1)
display(sol.theta_ddot2)

X = sym('X',[4,1]);
% These equations give the result in terms of state space coordinates X1,
% X2, X3 and X4
theta_ddot1f = subs(sol.theta_ddot1,[diff(theta_1),diff(theta_2),theta_1,theta_2],[X(2),X(4),X(1),X(3)]);
theta_ddot2f = subs(sol.theta_ddot2,[diff(theta_1),diff(theta_2),theta_1,theta_2],[X(2),X(4),X(1),X(3)]);

%% Converting to state space and making ODE
T=10; % seconds
initial_cond = [deg2rad(30); 0; deg2rad(45); 0];
X = sym('X',[4,1]);
X(1) = theta_1; % theta_1
X(2) = diff(theta_1,t); % theta_dot1
X(3) = theta_2; % theta_2
X(4) = diff(theta_2,t); % theta_dot2
[t,x] = ode45(@(t,X) ode_variables(t,X),[0 T],initial_cond);

%% Generating Plots
plot(t,x(:,1))
xlabel('Time (seconds)')
ylabel('theta_1 (radians)')
plot(t,x(:,2))
xlabel('Time (seconds)')
ylabel('theta_1dot (radians/second)')
plot(t,x(:,3))
xlabel('Time (seconds)')
ylabel('theta_2 (radians)')
plot(t,x(:,4))
xlabel('Time (seconds)')
ylabel('theta_2dot (radians/second)')

%% Creating a function to input variable values
function dX = ode_variables(t,X)
m1=1;
m2=1;
T1=0;
T2=0;
I1=0.084;
I2=0.084;
r1=0.45;
r2=0.45;
l1=1;
l2=1;
g = 9.81;

X1 = X(1);
X2 = X(2);
X3 = X(3);
X4 = X(4);
dX = zeros(4,1);
dX(1) = X2;
dX(2) = (I2*T1 - I2*T2 + T1*m2*r2^2 - T2*m2*r2^2 + g*l1*m2^2*r2^2*sin(X1) + I2*g*l1*m2*sin(X1) + I2*g*m1*r1*sin(X1) - T2*l1*m2*r2*cos(X3) + X2^2*l1*m2^2*r2^3*sin(X3) + X4^2*l1*m2^2*r2^3*sin(X3) + X2^2*l1^2*m2^2*r2^2*cos(X3)*sin(X3) - g*l1*m2^2*r2^2*sin(X3 + X1)*cos(X3) + I2*X2^2*l1*m2*r2*sin(X3) + I2*X4^2*l1*m2*r2*sin(X3) + g*m1*m2*r1*r2^2*sin(X1) + 2*X2*X4*l1*m2^2*r2^3*sin(X3) + 2*I2*X2*X4*l1*m2*r2*sin(X3))/(I1*I2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + I2*m1*r1^2 + I1*m2*r2^2 + m1*m2*r1^2*r2^2 - l1^2*m2^2*r2^2*cos(X3)^2);
dX(3) = X4;
dX(4) = -(I2*T1 - I1*T2 - I2*T2 - T2*l1^2*m2 - T2*m1*r1^2 + T1*m2*r2^2 - T2*m2*r2^2 - g*l1^2*m2^2*r2*sin(X3 + X1) - I1*g*m2*r2*sin(X3 + X1) + g*l1*m2^2*r2^2*sin(X1) + I2*g*l1*m2*sin(X1) + I2*g*m1*r1*sin(X1) + T1*l1*m2*r2*cos(X3) - 2*T2*l1*m2*r2*cos(X3) + X2^2*l1*m2^2*r2^3*sin(X3) + X2^2*l1^3*m2^2*r2*sin(X3) + X4^2*l1*m2^2*r2^3*sin(X3) + 2*X2^2*l1^2*m2^2*r2^2*cos(X3)*sin(X3) + X4^2*l1^2*m2^2*r2^2*cos(X3)*sin(X3) - g*l1*m2^2*r2^2*sin(X3 + X1)*cos(X3) + g*l1^2*m2^2*r2*cos(X3)*sin(X1) - g*m1*m2*r1^2*r2*sin(X3 + X1) + I1*X2^2*l1*m2*r2*sin(X3) + I2*X2^2*l1*m2*r2*sin(X3) + I2*X4^2*l1*m2*r2*sin(X3) + g*m1*m2*r1*r2^2*sin(X1) + 2*X2*X4*l1*m2^2*r2^3*sin(X3) + X2^2*l1*m1*m2*r1^2*r2*sin(X3) + 2*I2*X2*X4*l1*m2*r2*sin(X3) + 2*X2*X4*l1^2*m2^2*r2^2*cos(X3)*sin(X3) + g*l1*m1*m2*r1*r2*cos(X3)*sin(X1))/(I1*I2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + I2*m1*r1^2 + I1*m2*r2^2 + m1*m2*r1^2*r2^2 - l1^2*m2^2*r2^2*cos(X3)^2);
end
