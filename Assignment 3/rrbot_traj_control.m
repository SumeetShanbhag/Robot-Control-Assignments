rosshutdown;
clear; close; clc;
% ROS Setup
rosinit;
syms t0 t1 
syms theta1 theta1_dot1 theta1_dot2 
syms theta2 theta2_dot1 theta2_dot2


m1 = 1; m2=1; l1=1; l2=1 ; r1=0.45; r2=0.45; 
I1=0.084 ; I2 = 0.084 ; g = 9.81;

j1_effort = rospublisher("/rrbot/joint1_effort_controller/command");
j2_effort = rospublisher("/rrbot/joint2_effort_controller/command");

JointStates = rossubscriber("/rrbot/joint_states");

tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);

tau1.Data = 0;
tau2.Data = 0;


send(j1_effort,tau1);
send(j2_effort,tau2);
    
client = rossvcclient("/gazebo/set_model_configuration");
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);

tic;
t = 0;
counter=1;

while(t < 10)
    t = toc;
    % Trajectory for Joint one
    Joint1 = (pi*t^3)/500 - (3*pi*t^2)/100 + pi;
    des_vel1 = (3*pi*t^2)/500 -(3*pi*t)/50;
    net_acc1 = (3*pi*t)/250 -(3*pi)/50;

%Trajectory for Joint 2
    Joint2 = (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2;
    des_vel2 = (3*pi*t^2)/1000 -(3*pi*t)/100;
    net_acc2 = (3*pi*t)/500 -(3*pi)/100; 

    % read the joint states
    jointData = receive(JointStates);
    X = [(jointData.Position(1));(jointData.Position(2));jointData.Velocity(1);jointData.Velocity(2)];
   g1(counter)=jointData.Position(1); % X(1)
   g2(counter)=jointData.Position(2); % X(2)
   g3(counter)=jointData.Velocity(1); % X(3)
   g4(counter)=jointData.Velocity(2); % % X(4)
    
   a1(counter) = Joint1;
   a2(counter) = des_vel1;
   a3(counter) = Joint2;
   a4(counter) = des_vel2;
   

    %K =  [30 0 11 0;  0 12 0 7 ];
    K = [48, 0 , 14, 0; 0, 90, 0, 19]; 
    theta1 = X(1);
    theta2 = X(2);
    theta1_dot1 = X(3);
    theta2_dot1 = X(4);

    des_X = [jointData.Position(1) - Joint1; jointData.Position(2)- Joint2;jointData.Velocity(1)-des_vel1; jointData.Velocity(2)-des_vel2];

    M = [    (81*cos(theta1)^2)/400 + (81*sin(theta1)^2)/400 + ((9*cos(theta1 + theta2))/20 + cos(theta1))^2 + ((9*sin(theta1 + theta2))/20 + sin(theta1))^2 + 21/125, (9*cos(theta1 + theta2)*((9*cos(theta1 + theta2))/20 + cos(theta1)))/20 + (9*sin(theta1 + theta2)*((9*sin(theta1 + theta2))/20 + sin(theta1)))/20 + 21/250 ; (9*cos(theta1 + theta2)*((9*cos(theta1 + theta2))/20 + cos(theta1)))/20 + (9*sin(theta1 + theta2)*((9*sin(theta1 + theta2))/20 + sin(theta1)))/20 + 21/250,                                                                                       (81*cos(theta1 + theta2)^2)/400 + (81*sin(theta1 + theta2)^2)/400 + 21/250];
    C = [-(theta2_dot1*((9*sin(theta1 + theta2)*(theta1_dot1*((9*cos(theta1 + theta2))/20 + cos(theta1)) + (9*theta2_dot1*cos(theta1 + theta2))/20))/10 + ((9*cos(theta1 + theta2))/20 + cos(theta1))*((9*theta1_dot1*sin(theta1 + theta2))/10 + (9*theta2_dot1*sin(theta1 + theta2))/10) - ((9*theta1_dot1*cos(theta1 + theta2))/10 + (9*theta2_dot1*cos(theta1 + theta2))/10)*((9*sin(theta1 + theta2))/20 + sin(theta1)) - (9*cos(theta1 + theta2)*(theta1_dot1*((9*sin(theta1 + theta2))/20 + sin(theta1)) + (9*theta2_dot1*sin(theta1 + theta2))/20))/10))/2;
     (((9*theta1_dot1*sin(theta1 + theta2))/10 + (9*theta2_dot1*sin(theta1 + theta2))/10)*(theta1_dot1*((9*cos(theta1 + theta2))/20 + cos(theta1)) + (9*theta2_dot1*cos(theta1 + theta2))/20))/2 - (theta2_dot1*((9*sin(theta1 + theta2)*(theta1_dot1*((9*cos(theta1 + theta2))/20 + cos(theta1)) + (9*theta2_dot1*cos(theta1 + theta2))/20))/10 - (9*cos(theta1 + theta2)*(theta1_dot1*((9*sin(theta1 + theta2))/20 + sin(theta1)) + (9*theta2_dot1*sin(theta1 + theta2))/20))/10 - (9*sin(theta1 + theta2)*((9*theta1_dot1*cos(theta1 + theta2))/20 + (9*theta2_dot1*cos(theta1 + theta2))/20))/10 + (9*cos(theta1 + theta2)*((9*theta1_dot1*sin(theta1 + theta2))/20 + (9*theta2_dot1*sin(theta1 + theta2))/20))/10))/2 - ((theta1_dot1*((9*sin(theta1 + theta2))/20 + sin(theta1)) + (9*theta2_dot1*sin(theta1 + theta2))/20)*((9*theta1_dot1*cos(theta1 + theta2))/10 + (9*theta2_dot1*cos(theta1 + theta2))/10))/2];
    G =  [- (8829*sin(theta1 + theta2))/2000 - (28449*sin(theta1))/2000;
    -(8829*sin(theta1 + theta2))/2000];
    
    
    
    v = [-K(1,:)*des_X + net_acc1;-K(2,:)*des_X + net_acc2 ];
    tau = M*v + C + G;
  
    tau1.Data = tau(1);
    tau2.Data = tau(2);
     
   send(j1_effort,tau1);
   send(j2_effort,tau2);
   
   t1(counter)=tau1.Data;
   t2(counter)=tau2.Data;  
   time(counter)=t;
   counter=counter+1;
end

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

disp("GAZEBO.....TRAJECTORIES")
subplot(1,3,1)

plot(time,g1)
title('Theta1(actual) ',fontsize=8)

subplot(1,3,2)
plot(time,a1)
title('Theta1(desired)',fontsize=8)

subplot(1,3,3)
plot(time,g1)
hold on
plot(time ,a1)
hold off
%legend('actual','desired')
title('actual & desired theta1',fontsize=8)
 figure;
subplot(1,3,1)
plot(time, g2)
title('Theta2(actual) ',fontsize=8)

subplot(1,3,2)
plot(time,a3)
title('Theta2(desired)',fontsize=8)

subplot(1,3,3)
plot(time,g2)
hold on
plot(time ,a3)
hold off
%legend('actual traj2','desired traj2')
title('actual & desired theta2',fontsize=8)

figure;

subplot(1,3,1)
plot(time, g3)
title('Theta1dot(actual) ',fontsize=8)

subplot(1,3,2)
plot(time,a2)
title('Theta1dot(desired)',fontsize=8)

subplot(1,3,3)
plot(time,g3)
hold on
plot(time ,a2)
hold off
%legend('actual traj3','desired traj3')
title('actual & desired theta1dot',fontsize=8)

figure;

subplot(1,3,1)
plot(time,g4)
title('Theta2dot(actual)',fontsize=8)

subplot(1,3,2)
plot(time,a4)
title('Theta2dot(desired)',fontsize=8)

subplot(1,3,3)
plot(time,g4)
hold on
plot(time ,a4)
hold off
%legend('actual traj4','desired traj4')
title('actual & desired theta2dot',fontsize=8)

figure;
subplot(1,2,1)
plot(time,t1(1,:))
title('U1 (Control Input)')

subplot(1,2,2)
plot(time,t2(1,:))
title('U2 ( Control Input)')
rosshutdown;