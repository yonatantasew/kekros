%link length
L1=0;L2=150;L3=200;L4=100;
%link parameters
a1=0;alpha1=0;d1=L1;theta1=0;
a2=0;alpha2=pi/2;d2=L2;theta2=0;
a3=0;alpha3=pi/2;d3=0;theta3=0;
a4=L3;alpha4=0;d4=0;theta4=0;
a5=L4;alpha5=0;d5=0;theta5=0;
%Eulers integration method
%theta_new = theta_old + (theta_old*delta_T)
L1_new = 0; theta2_new = 0; theta3_new = 0; theta4_new = 0;
x_desir = [1;2;3];
time_taken = 0.3;

while 1

    %assume change in time = 2 sec
    L1 = L1 + (L1_new*2);
    theta2 = theta2 + (theta2_new*time_taken);
    theta3 = theta3 + (theta3_new*time_taken);
    theta4 = theta4 + (theta4_new*time_taken);
    %Forward Kinematics
    T01 = DH_of(a1,alpha1,d1,theta1);
    T12 = DH_of(a2,alpha2,d2,theta2);
    T23 = DH_of(a3,alpha3,d3,theta3);
    T34 = DH_of(a4,alpha4,d4,theta4);
    T45 = DH_of(a5,alpha5,d5,theta5);
    T02 = T01*T12;
    T03 = T02*T23;
    T04 = T03*T34;
    T05 = T04*T45;
%   cla;
%    line([0 0], [0 0], [0 L1],'color', [1 0 0], 'LineWidth', 3);hold on;
%    line([0 T01(1,4)], [0 T01(2,4)],[L1 T01(3, 4)],'color',[1 0 0], 'LineWidth',3);
%    line([T01(1,4) T02(1,4)], [T01(2,4) T02(2,4)],[T01(3,4) T02(3,4)],'color',[0 1 1], 'LineWidth',3);
%    line([T02(1,4) T03(1,4)], [T02(2,4) T03(2,4)],[T02(3,4) T03(3,4)],'color',[1 1 0], 'LineWidth',3);
%    line([T03(1,4) T04(1,4)], [T03(2,4) T04(2,4)],[T03(3,4) T04(3,4)],'color',[1 0 0], 'LineWidth',3);
%    line([T04(1,4) T05(1,4)], [T04(2,4) T05(2,4)],[T04(3,4) T05(3,4)],'color',[0 0 1], 'LineWidth',3);

%    plot3(T05(:,1),T05(:,2),T05(:,3),'.m','LineWidth',2)
%    view([-54 16])
%    axis([-450 450 -300 300 0 200])
%    grid on
    %pause(0.1)

    %Velocity Kinematics
    Z0=[0;0;1];p0=[0;0;0];pe=T05(1:3,4);o=[0;0;0];
    Z1=T01(1:3,3);p1=T01(1:3,4);
    Z2=T02(1:3,3);p2=T02(1:3,4);
    Z3=T03(1:3,3);p3=T03(1:3,4);
    Z4=T04(1:3,3);p4=T04(1:3,4);
    J1 = Z0;
    J2 = cross(Z1,(pe-p1));
    J3 =cross(Z2,(pe-p2));
    J4 = cross(Z3,(pe-p3));
    J = [J1 J2 J3 J4;o Z1 Z2 Z3];

    %Psuedo_inverse of The Jacobian
    % q_dot = pinv(J)*x_dot
    % position_error = x_desired - x_current
    x_current = T05(1:3,4);
    x_desired = x_desir;
    x_desired = x_desired - x_current;
    q_dot = pinv(J)*x_desired;
    phi_current = atan2d(x_current(1,1),x_current(2,1));
    phi_desired = atan2d(x_desired(1,1),x_desired(2,1));
    distance_error=sqrt((x_desired(1,1))^2+(x_desired(2,1))^2)- sqrt((x_current(1,1))^2+(x_current(2,1))^2);
    if abs(distance_error)<=0.2 && abs(x_desired)<=2
       break;
    end


    L1_new = q_dot(1,1)
    tht2 = q_dot(2,1);
    tht3 = q_dot(3,1);
    tht4 = q_dot(4,1);
    theta2_new = radtodeg(tht2)
    theta3_new = radtodeg(tht3)
    theta4_new = radtodeg(tht4)

end
    
    

