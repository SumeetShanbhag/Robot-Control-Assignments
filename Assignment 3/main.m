%VARIABLES
clc;
clear;
close all; 
syms theta1_initial theta1_final theta1dot_initial theta1dot_final theta2_initial theta2_final theta2dot_initial theta2dot_final 'real'

syms theta1 theta2  theta1_dot theta2_dot theta1_ddot theta2_ddot T1  T2 'real'
syms r1 r2 l1 l2 I1 I2 m1 m2 g
%PHYSICAL PARAMETERS
%m1=1; m2=1; l1= 0.084; l2 = 0.084; r1 = 0.45; r2= 0.45; I1=0.084; I2=0.084; g = 9.81;
syms t0 t1 t
%Time Intervals
to = 0;
tf = 10;
%%
A = [1 to to^2 to^3;
    0 1 2*to 3*to^2; 
    1 tf tf^2 tf^3;
    0 1  2*tf 3*(tf^2)];

% JOINT 1
B1 = [ pi; 0; 0; 0];
J1 = inv(A)*B1;

% JOINT 2
B2 = [ pi/2; 0; 0; 0];
J2 = inv(A)*B2;

% Trjacetories for joint1 and joint2

Joint1 = [1 t t^2 t^3]*J1;
vel1 = jacobian(Joint1,t);
acc1 = jacobian(vel1,t);


Joint2 = [1 t t^2 t^3]*J2;
vel2 = jacobian(Joint2,t);
acc2 = jacobian(vel2,t);
%%
% MANIPULATOR FORM OF EQUATION

T1 = theta1_ddot*(I1 + I2 + (m1*(2*r1^2*cos(theta1)^2 + 2*r1^2*sin(theta1)^2))/2 + (m2*(2*(r2*cos(theta2 + theta1) + l1*cos(theta1))^2 + 2*(r2*sin(theta2 + theta1) + l1*sin(theta1))^2))/2) + theta2_ddot*(I2 + (m2*(2*r2*sin(theta2 + theta1)*(r2*sin(theta2 + theta1) + l1*sin(theta1)) + 2*r2*cos(theta2 + theta1)*(r2*cos(theta2 + theta1) + l1*cos(theta1))))/2) + (m2*theta2_dot*(2*(cos(theta1 + theta2)*r2*theta2_dot + cos(theta1 + theta2)*r2*theta1_dot)*(r2*sin(theta2 + theta1) + l1*sin(theta1)) - 2*(sin(theta1 + theta2)*r2*theta2_dot + sin(theta1 + theta2)*r2*theta1_dot)*(r2*cos(theta2 + theta1) + l1*cos(theta1)) - 2*r2*sin(theta2 + theta1)*(theta1_dot*(r2*cos(theta2 + theta1) + l1*cos(theta1)) + r2*theta2_dot*cos(theta2 + theta1)) + 2*r2*cos(theta2 + theta1)*(theta1_dot*(r2*sin(theta2 + theta1) + l1*sin(theta1)) + r2*theta2_dot*sin(theta2 + theta1))))/2 - g*m2*(r2*sin(theta2 + theta1) + l1*sin(theta1)) - g*m1*r1*sin(theta1);
T1 = subs(T1,[m1,m2,r1,r2,I1,I2,g,l1,l2],[1,1,0.45,0.45,0.084,0.084,9.81,1,1]);


T2 = (m2*(2*(sin(theta1 + theta2)*r2*theta2_dot + sin(theta1 + theta2)*r2*theta1_dot)*(theta1_dot*(r2*cos(theta2 + theta1) + l1*cos(theta1)) + r2*theta2_dot*cos(theta2 + theta1)) - 2*(cos(theta1 + theta2)*r2*theta2_dot + cos(theta1 + theta2)*r2*theta1_dot)*(theta1_dot*(r2*sin(theta2 + theta1) + l1*sin(theta1)) + r2*theta2_dot*sin(theta2 + theta1))))/2 + theta2_ddot*(I2 + (m2*(2*r2^2*cos(theta2 + theta1)^2 + 2*r2^2*sin(theta2 + theta1)^2))/2) + theta1_ddot*(I2 + (m2*(2*r2*sin(theta2 + theta1)*(r2*sin(theta2 + theta1) + l1*sin(theta1)) + 2*r2*cos(theta2 + theta1)*(r2*cos(theta2 + theta1) + l1*cos(theta1))))/2) - (m2*theta2_dot*(2*r2*sin(theta2 + theta1)*(theta1_dot*(r2*cos(theta2 + theta1) + l1*cos(theta1)) + r2*theta2_dot*cos(theta2 + theta1)) - 2*r2*cos(theta2 + theta1)*(theta1_dot*(r2*sin(theta2 + theta1) + l1*sin(theta1)) + r2*theta2_dot*sin(theta2 + theta1)) - 2*r2*sin(theta2 + theta1)*(cos(theta1 + theta2)*r2*theta2_dot + cos(theta1 + theta2)*r2*theta1_dot) + 2*r2*cos(theta2 + theta1)*(sin(theta1 + theta2)*r2*theta2_dot + sin(theta1 + theta2)*r2*theta1_dot)))/2 - g*m2*r2*sin(theta2 + theta1);
T2 = subs(T2,[m1,m2,r1,r2,I1,I2,g,l1,l2],[1,1,0.45,0.45,0.084,0.084,9.81,1,1]);

G(1,1) =subs(T1,[theta1_ddot,theta2_ddot,theta1_dot,theta2_dot],[0,0,0,0]);
G(1,2) =subs(T2,[theta1_ddot,theta2_ddot,theta1_dot,theta2_dot],[0,0,0,0]);
G = simplify(G)

% isolating M MATRIX
M(1,1) = subs((T1-G(1,1)),[theta1_dot,theta1_ddot,theta2_dot,theta2_ddot],[0,1,0,0]);
M(2,1) = subs((T2-G(1,2)),[theta1_dot,theta1_ddot,theta2_dot,theta2_ddot],[0,1,0,0]);
M(1,2) = subs((T1-G(1,1)),[theta1_dot,theta1_ddot,theta2_dot,theta2_ddot],[0,0,0,1]);
M(2,2) = subs((T2-G(1,2)),[theta1_dot,theta1_ddot,theta2_dot,theta2_ddot],[0,0,0,1]);

M = simplify(M)
%M = subs(M,[theta1_ddot,theta2_ddot],[1,1])

C(1,1) = T1 - M(1,1)*theta1_ddot - M(1,2)*theta2_ddot- G(1,1);
C(2,1) = T2 - M(2,1)*theta1_ddot - M(2,2)*theta2_ddot - G(1,2);
C = simplify(C); % equation in the form of theta1_dot and theta2_dot
% No need to isolate C matrix
%% 

% defining feedback linearization using eigen value placement
v = [theta1_ddot;theta2_ddot];
X = [ theta1;theta2;theta1_dot;theta2_dot];
%state space form
dX = [theta1_dot;theta2_dot;v(1);v(2)];

A = jacobian(dX,X);
A = double(A);
B = jacobian(dX,v);
B = double(B);

%EIGEN VALUE PLACEMENT TO CHOOSE K VALUES
lamdas = [ -10,-8,-7,-6];
K = place(A,B,lamdas);

%%
% CALLING ODE FUNCTION

init = [deg2rad(200), deg2rad(125), 0, 0];
[t, y] = ode45(@ode_rrbot, [0, 10], init);
 

% Desired Trajectory 
int  = 0:0.5:10;

%Trajectory for Joint 1
Joint1 = (pi*int.^3)/500 - (3*pi*int.^2)/100 + pi;
des_vel1 = (3*pi*int.^2)/500 -(3*pi*int)/50;
des_acc1 = (3*pi*int)/250 -(3*pi)/50;

%Trajectory for Joint 2
Joint2 = (pi*int.^3)/1000 - (3*pi*int.^2)/200 + pi/2;
des_vel2 = (3*pi*int.^2)/1000 -(3*pi*int)/100;
des_acc2 = (3*pi*int)/500 -(3*pi)/100;

% Expected/Experimental Values 
tau = [];

for i=1:length(t)
    
    exp_theta1 = y(i,1);
    exp_theta2 = y(i,2);
    exp_theta1_dot = y(i,3);
    exp_theta2_dot = y(i,4);

    des_theta1 = pi + (pi*t(i).^3)/500 - (3*pi*t(i).^2)/100 ;
    des_theta1_dot = (3*pi*t(i).^2)/500 - (3*pi*t(i))/50;

    des_theta2 = pi/2 + (pi*t(i).^3)/1000 - (3*pi*t(i).^2)/200 ;
    des_theta2_dot = (3*pi*t(i).^2)/1000 - (3*pi*t(i))/100;

    des_acc1 = (3*pi*t(i))/250 - (3*pi)/50;
    des_acc2= (3*pi*t(i))/500 - (3*pi)/100;

    des_M = double(subs(M,[theta1,theta2,theta1_dot,theta2_dot],[exp_theta1,exp_theta2,exp_theta1_dot,exp_theta2_dot]));
    des_C = double(subs(C,[theta1,theta2,theta1_dot,theta2_dot],[exp_theta1,exp_theta2,exp_theta1_dot,exp_theta2_dot]));
    des_G = double(subs(G,[theta1,theta2,theta1_dot,theta2_dot],[exp_theta1,exp_theta2,exp_theta1_dot,exp_theta2_dot]));

    total_acc = [des_acc1;des_acc2];
    des_v = - (K * ([exp_theta1;exp_theta2;exp_theta1_dot;exp_theta2_dot]-[des_theta1;des_theta2;des_theta1_dot;des_theta2_dot])) + total_acc;

    des_u = des_M*des_v + des_C + des_G;
    tau = [tau; des_u(1) des_u(2)];
end    

disp("MATLAB TRAJECTORIES");
subplot(1,3,1)

plot(t,y(:,1))
title('Theta1(Expected)',fontsize=8)

subplot(1,3,2)
plot(int,Joint1)
title('Theta1(desired)',fontsize=8)

subplot(1,3,3)
plot(t,y(:,1))
hold on
plot(int ,Joint1)
hold off
title('Expected & desired theta1',fontsize=8)

figure;
subplot(1,3,1)
plot(t,y(:,2))
title('Theta2(Expected)',fontsize=8)

subplot(1,3,2)
plot(int,Joint2)
title('Theta2(desired)',fontsize=8)

subplot(1,3,3)
plot(t,y(:,2))
hold on
plot(int ,Joint2)
hold off
title('Expected & desired theta2',fontsize=8)

figure;

subplot(1,3,1)
plot(t,y(:,3))
title('Theta1dot(Expected)',fontsize=8)

subplot(1,3,2)
plot(int,des_vel1)
title('Theta1dot(desired)',fontsize=8)

subplot(1,3,3)
plot(t,y(:,3))
hold on
plot(int ,des_vel1)
hold off
title('Expected & desired Theta1dot',fontsize=8)

figure;

subplot(1,3,1)
plot(t,y(:,4))
title('Theta2dot(Expected) ',fontsize=8)

subplot(1,3,2)
plot(int,des_vel2)
title('Theta2dot(desired)',fontsize=8)

subplot(1,3,3)
plot(t,y(:,4))
hold on
plot(int ,des_vel2)
hold off
title('Expected & Desired Theta2dot',fontsize=8)

figure;
subplot(1,2,1)
plot(t,tau(:,1))
title('U1')
subplot(1,2,2)
plot(t,tau(:,2))
title('U2 ( Control Input)')
%%
function dX=ode_rrbot(t,X)
    
    m1=1; 
    m2=1; 
    l1=1;
    l2=1;
    r1=0.45; 
    r2=0.45;
    I1=0.084; 
    I2=0.084;
    g = 9.81;
   
    
    dX=zeros(4,1);
    X=num2cell(X);
    
    [theta1, theta2, theta1_dot1, theta2_dot1]=deal(X{:});

    if abs(theta1)>2*pi
        theta1= mod(theta1,2*pi);
    end

    if abs(theta2)>2*pi
        theta2 = mod(theta2,2*pi) ;   
    end  

    K =  [30 0 11 0;  0 12 0 7 ];

   Joint1 = (pi*t^3)/500 - (3*pi*t^2)/100 + pi;
   des_vel1 = (3*pi*t^2)/500 -(3*pi*t)/50;
   net_acc1 = (3*pi*t)/250 -(3*pi)/50;

%Trajectory for Joint 2
   Joint2 = (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2;
   des_vel2 = (3*pi*t^2)/1000 -(3*pi*t)/100;
   net_acc2 = (3*pi*t)/500 -(3*pi)/100;

   net_acc = [ net_acc1;net_acc2];
   des_X = [theta1 - Joint1;theta2 - Joint2;theta1_dot1-des_vel1;theta2_dot1-des_vel2];

   %manipulator form 
   M = [    (81*cos(theta1)^2)/400 + (81*sin(theta1)^2)/400 + ((9*cos(theta1 + theta2))/20 + cos(theta1))^2 + ((9*sin(theta1 + theta2))/20 + sin(theta1))^2 + 21/125, (9*cos(theta1 + theta2)*((9*cos(theta1 + theta2))/20 + cos(theta1)))/20 + (9*sin(theta1 + theta2)*((9*sin(theta1 + theta2))/20 + sin(theta1)))/20 + 21/250 ; (9*cos(theta1 + theta2)*((9*cos(theta1 + theta2))/20 + cos(theta1)))/20 + (9*sin(theta1 + theta2)*((9*sin(theta1 + theta2))/20 + sin(theta1)))/20 + 21/250,                                                                                       (81*cos(theta1 + theta2)^2)/400 + (81*sin(theta1 + theta2)^2)/400 + 21/250];
   C = [-(theta2_dot1*((9*sin(theta1 + theta2)*(theta1_dot1*((9*cos(theta1 + theta2))/20 + cos(theta1)) + (9*theta2_dot1*cos(theta1 + theta2))/20))/10 + ((9*cos(theta1 + theta2))/20 + cos(theta1))*((9*theta1_dot1*sin(theta1 + theta2))/10 + (9*theta2_dot1*sin(theta1 + theta2))/10) - ((9*theta1_dot1*cos(theta1 + theta2))/10 + (9*theta2_dot1*cos(theta1 + theta2))/10)*((9*sin(theta1 + theta2))/20 + sin(theta1)) - (9*cos(theta1 + theta2)*(theta1_dot1*((9*sin(theta1 + theta2))/20 + sin(theta1)) + (9*theta2_dot1*sin(theta1 + theta2))/20))/10))/2;
    (((9*theta1_dot1*sin(theta1 + theta2))/10 + (9*theta2_dot1*sin(theta1 + theta2))/10)*(theta1_dot1*((9*cos(theta1 + theta2))/20 + cos(theta1)) + (9*theta2_dot1*cos(theta1 + theta2))/20))/2 - (theta2_dot1*((9*sin(theta1 + theta2)*(theta1_dot1*((9*cos(theta1 + theta2))/20 + cos(theta1)) + (9*theta2_dot1*cos(theta1 + theta2))/20))/10 - (9*cos(theta1 + theta2)*(theta1_dot1*((9*sin(theta1 + theta2))/20 + sin(theta1)) + (9*theta2_dot1*sin(theta1 + theta2))/20))/10 - (9*sin(theta1 + theta2)*((9*theta1_dot1*cos(theta1 + theta2))/20 + (9*theta2_dot1*cos(theta1 + theta2))/20))/10 + (9*cos(theta1 + theta2)*((9*theta1_dot1*sin(theta1 + theta2))/20 + (9*theta2_dot1*sin(theta1 + theta2))/20))/10))/2 - ((theta1_dot1*((9*sin(theta1 + theta2))/20 + sin(theta1)) + (9*theta2_dot1*sin(theta1 + theta2))/20)*((9*theta1_dot1*cos(theta1 + theta2))/10 + (9*theta2_dot1*cos(theta1 + theta2))/10))/2];
   G =  [- (8829*sin(theta1 + theta2))/2000 - (28449*sin(theta1))/2000;
    -(8829*sin(theta1 + theta2))/2000];
    
   v = [-K(1,:)*des_X + net_acc1;-K(2,:)*des_X + net_acc2 ];

   %v = - (K * ([theta1;theta2;theta1_dot1;theta2_dot1]-[Joint1;Joint2;des_vel1;des_vel2])) + net_acc;
   v;
   U = M*v + C + G;
   t1 = U(1);
   t2 = U(2);

   theta1_dot2=(I2*t1 - I2*t2 + m2*r2^2*t1*cos(theta1 + theta2)^2 - m2*r2^2*t2*cos(theta1 + theta2)^2 + m2*r2^2*t1*sin(theta1 + theta2)^2 - m2*r2^2*t2*sin(theta1 + theta2)^2 + g*I2*l1*m2*sin(theta1) + g*I2*m1*r1*sin(theta1) - l1*m2^2*r2^3*theta1_dot1^2*cos(theta1 + theta2)^3*sin(theta1) + l1*m2^2*r2^3*theta1_dot1^2*sin(theta1 + theta2)^3*cos(theta1) - l1*m2^2*r2^3*theta2_dot1^2*cos(theta1 + theta2)^3*sin(theta1) + l1*m2^2*r2^3*theta2_dot1^2*sin(theta1 + theta2)^3*cos(theta1) - l1*m2*r2*t2*cos(theta1 + theta2)*cos(theta1) - l1*m2*r2*t2*sin(theta1 + theta2)*sin(theta1) + g*l1*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1) - g*l1*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1) + l1^2*m2^2*r2^2*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 - 2*l1*m2^2*r2^3*theta1_dot1*theta2_dot1*cos(theta1 + theta2)^3*sin(theta1) + 2*l1*m2^2*r2^3*theta1_dot1*theta2_dot1*sin(theta1 + theta2)^3*cos(theta1) - l1^2*m2^2*r2^2*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 - l1^2*m2^2*r2^2*theta1_dot1^2*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) + l1^2*m2^2*r2^2*theta1_dot1^2*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) - I2*l1*m2*r2*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1) + I2*l1*m2*r2*theta1_dot1^2*sin(theta1 + theta2)*cos(theta1) - I2*l1*m2*r2*theta2_dot1^2*cos(theta1 + theta2)*sin(theta1) + I2*l1*m2*r2*theta2_dot1^2*sin(theta1 + theta2)*cos(theta1) + g*m1*m2*r1*r2^2*cos(theta1 + theta2)^2*sin(theta1) + l1*m2^2*r2^3*theta1_dot1^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) + l1*m2^2*r2^3*theta2_dot1^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) + g*m1*m2*r1*r2^2*sin(theta1 + theta2)^2*sin(theta1) - l1*m2^2*r2^3*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - l1*m2^2*r2^3*theta2_dot1^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - 2*l1*m2^2*r2^3*theta1_dot1*theta2_dot1*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - 2*I2*l1*m2*r2*theta1_dot1*theta2_dot1*cos(theta1 + theta2)*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot1*theta2_dot1*sin(theta1 + theta2)*cos(theta1) + 2*l1*m2^2*r2^3*theta1_dot1*theta2_dot1*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1))/(I1*I2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + l1^2*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + l1^2*m2^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*sin(theta1)^2 - 2*l1^2*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)*sin(theta1));
   theta2_dot2 = (I1*t2 - I2*t1 + I2*t2 - m2*r2^2*t1*cos(theta1 + theta2)^2 + m2*r2^2*t2*cos(theta1 + theta2)^2 - m2*r2^2*t1*sin(theta1 + theta2)^2 + m2*r2^2*t2*sin(theta1 + theta2)^2 + l1^2*m2*t2*cos(theta1)^2 + m1*r1^2*t2*cos(theta1)^2 + l1^2*m2*t2*sin(theta1)^2 + m1*r1^2*t2*sin(theta1)^2 + g*I1*m2*r2*sin(theta1 + theta2) - g*I2*l1*m2*sin(theta1) - g*I2*m1*r1*sin(theta1) + l1*m2^2*r2^3*theta1_dot1^2*cos(theta1 + theta2)^3*sin(theta1) - l1*m2^2*r2^3*theta1_dot1^2*sin(theta1 + theta2)^3*cos(theta1) + l1^3*m2^2*r2*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1)^3 - l1^3*m2^2*r2*theta1_dot1^2*sin(theta1 + theta2)*cos(theta1)^3 + l1*m2^2*r2^3*theta2_dot1^2*cos(theta1 + theta2)^3*sin(theta1) - l1*m2^2*r2^3*theta2_dot1^2*sin(theta1 + theta2)^3*cos(theta1) - l1*m2*r2*t1*cos(theta1 + theta2)*cos(theta1) + 2*l1*m2*r2*t2*cos(theta1 + theta2)*cos(theta1) - l1*m2*r2*t1*sin(theta1 + theta2)*sin(theta1) + 2*l1*m2*r2*t2*sin(theta1 + theta2)*sin(theta1) - g*l1*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1) + g*l1^2*m2^2*r2*sin(theta1 + theta2)*cos(theta1)^2 + l1^3*m2^2*r2*theta1_dot1^2*cos(theta1 + theta2)*cos(theta1)^2*sin(theta1) + g*l1*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1) - l1^3*m2^2*r2*theta1_dot1^2*sin(theta1 + theta2)*cos(theta1)*sin(theta1)^2 - 2*l1^2*m2^2*r2^2*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 - l1^2*m2^2*r2^2*theta2_dot1^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 - g*l1^2*m2^2*r2*cos(theta1 + theta2)*cos(theta1)*sin(theta1) + 2*l1*m2^2*r2^3*theta1_dot1*theta2_dot1*cos(theta1 + theta2)^3*sin(theta1) - 2*l1*m2^2*r2^3*theta1_dot1*theta2_dot1*sin(theta1 + theta2)^3*cos(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 + l1^2*m2^2*r2^2*theta2_dot1^2*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 + 2*l1^2*m2^2*r2^2*theta1_dot1^2*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) + l1^2*m2^2*r2^2*theta2_dot1^2*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) - 2*l1^2*m2^2*r2^2*theta1_dot1^2*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) - l1^2*m2^2*r2^2*theta2_dot1^2*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) + I1*l1*m2*r2*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1) - I1*l1*m2*r2*theta1_dot1^2*sin(theta1 + theta2)*cos(theta1) + I2*l1*m2*r2*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1) - I2*l1*m2*r2*theta1_dot1^2*sin(theta1 + theta2)*cos(theta1) + I2*l1*m2*r2*theta2_dot1^2*cos(theta1 + theta2)*sin(theta1) - I2*l1*m2*r2*theta2_dot1^2*sin(theta1 + theta2)*cos(theta1) - g*m1*m2*r1*r2^2*cos(theta1 + theta2)^2*sin(theta1) + g*m1*m2*r1^2*r2*sin(theta1 + theta2)*cos(theta1)^2 - l1*m2^2*r2^3*theta1_dot1^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) - l1*m2^2*r2^3*theta2_dot1^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) - g*m1*m2*r1*r2^2*sin(theta1 + theta2)^2*sin(theta1) + g*m1*m2*r1^2*r2*sin(theta1 + theta2)*sin(theta1)^2 + l1*m2^2*r2^3*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) + l1*m2^2*r2^3*theta2_dot1^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) + 2*l1*m2^2*r2^3*theta1_dot1*theta2_dot1*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - 2*l1^2*m2^2*r2^2*theta1_dot1*theta2_dot1*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 + 2*l1^2*m2^2*r2^2*theta1_dot1*theta2_dot1*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 + 2*l1^2*m2^2*r2^2*theta1_dot1*theta2_dot1*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) - 2*l1^2*m2^2*r2^2*theta1_dot1*theta2_dot1*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot1*theta2_dot1*cos(theta1 + theta2)*sin(theta1) - 2*I2*l1*m2*r2*theta1_dot1*theta2_dot1*sin(theta1 + theta2)*cos(theta1) + l1*m1*m2*r1^2*r2*theta1_dot1^2*cos(theta1 + theta2)*sin(theta1)^3 - l1*m1*m2*r1^2*r2*theta1_dot1^2*sin(theta1 + theta2)*cos(theta1)^3 - 2*l1*m2^2*r2^3*theta1_dot1*theta2_dot1*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) - g*l1*m1*m2*r1*r2*sin(theta1 + theta2)*sin(theta1)^2 - g*l1*m1*m2*r1*r2*cos(theta1 + theta2)*cos(theta1)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot1^2*cos(theta1 + theta2)*cos(theta1)^2*sin(theta1) - l1*m1*m2*r1^2*r2*theta1_dot1^2*sin(theta1 + theta2)*cos(theta1)*sin(theta1)^2)/(I1*I2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + l1^2*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + l1^2*m2^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*sin(theta1)^2 - 2*l1^2*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)*sin(theta1));
   
   dX(1)= theta1_dot1;
   dX(2)= theta2_dot1;
   dX(3)= theta1_dot2;
   dX(4)= theta2_dot2;
   
   
end