close all;clear;
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
desired_location = [11;22;33];
time_taken = 0.3;
a = 0;b=0;
while a==0
    while b==0;
        %assume change in time 2 sec
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
        P0=[0 0];
        P1=transpose(T01(1:2,4));
        P2=transpose(T02(1:2,4));
        P3=transpose(T03(1:2,4));
        P4=transpose(T04(1:2,4));
        P5=transpose(T05(1:2,4));
        
        Q1=[P0(1,1) P1(1,1) P2(1,1) P3(1,1) P4(1,1) P5(1,1)];
        Q2=[P0(1,2) P1(1,2) P2(1,2) P3(1,2) P4(1,2) P5(1,2)];
        %Q3=[P0(1,1) P1(1,1) P2(1,1) P3(1,1)];
        plot(Q1,Q2,'-o','LineWidth',4);
        axis([-31,31,-31,31]);
        grid on;
        
        %Velocity Kinematics
        Z0=[0;0;1];p0=[0;0;0];pe=T05(1:3,4);
        Z1=T01(1:3,3);p1=T01(1:3,4);
        Z2=T02(1:3,3);p2=T02(1:3,4);
        Z3=T03(1:3,3);p3=T03(1:3,4);
        Z4=T04(1:3,3);p4=T04(1:3,4);
        J1 = cross(Z0,(pe-p0));
        J2 = cross(Z1,(pe-p1));
        J3 =cross(Z2,(pe-p2));
        J4 = cross(Z3,(pe-p3));
        J = [J1 J2 J3 J4];
        
        x_current = P5(1,3);y_current=P5(2,3);
        x_desired = desired_location(1,1);y_desired = desired_location(2,1);
        x_velocity = x_desired - x_current;
        y_velocity = y_desired - y_current;
        q_dot = pinv(J)*[x_velocity;y_velocity;0];
        phi_desired = atan2d(y_desired,x_desired);
        phi_current = atan2d(y_current,x_current);
        distance_error=sqrt(x_desired^2+y_desired^2)- sqrt(x_current^2+y_current^2);
        phi_error=phi_desired-phi_current;
        if abs(distance_error)<=0.2 && abs(phi_error)<=2
            b=1;
        end
        L1_new = q_dot(1,1);
        tht2 = q_dot(2,1);
        tht3 = q_dot(3,1);
        tht4 = q_dot(4,1);
        theta2_new = radtodeg(tht2);
        theta3_new = radtodeg(tht3);
        theta4_new = radtodeg(tht4);
         text(P4(1,1),P4(1,2),['  (', num2str(P4(1,1),3), ', ', num2str(P4(1,2),3), ')']);
        text(-25,-17,'Orinerror:','Color','red','FontSize',12)
        text(-25,-20,num2str(abs(phi_error),3),'Color','red','FontSize',12)
        text(-25,-23,'diserror:','Color','red','FontSize',12)
        text(-25,-26,num2str(abs(distance_error),3),'Color','red','FontSize',12)
        pause(0.01);
    end
    if b==1
       [desired_location(1,1),desired_location(2,1),buttons] = ginput;
       b=0;
    end
end
