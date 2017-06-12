%--------------------------TRAJECTORY PLANNING--------------------------
close all;clear;
%link length
L1=0;L2=0.15;L3=0.2;L4=0.1;
%link parameters
a1=0;alpha1=0;d1=L1;theta1=0;
a2=0;alpha2=pi/2;d2=L2;theta2=0;
a3=0;alpha3=pi/2;d3=0;theta3=0;
a4=L3;alpha4=0;d4=0;theta4=0;
a5=L4;alpha5=0;d5=0;theta5=0;
init_pos = [0,0,0,0];
home_position = [0.1,pi/2,5*pi/3,pi/2];
q = init_pos;
qf = home_position;
%ti=0;
%tht2 = theta2;
%tht3 = theta3;
%tht4 = theta4;
ti = 0;
ts = 0.01;% 10ms
tf = 3; 
pos =[];
for t=ti:ts:tf
    %Forward Kinematics
    T = forwardKinematics(L1,L2,L3,L4,theta2,theta3,theta4);
    pex = T(13);
    pey = T(14);
    pez = T(15);
    %Inverse Kinematics
    q(2) = atan2(pex,pey);
    q(2) = real(q(2))
    %using cosine law of triangle
    num = pex^2+pey^2-L2^2-L3^2-L4^2;
    den = 2*L3*L4;
    A = num/den;
    A = real(A)
    D = sqrt(1-A^2);
    D = real(D)
    q(4) = atan2(A,D);
    q(4) = real(q(4))
    B = atan(pez/sqrt(pex^2+pey^2));
    C = atan(L4*sin(theta4)/(L3+L4*cos(theta4)));
    q(3) = B-C;
    q(3) = real(q(3))
    pos = [pos;q];
    hold on;
    grid on;
    plot3(pos(:,2),pos(:,3),pos(:,4),'--m')

end
    


%%
syms L1 L2 L3 L4 tht2 tht3 tht4 
J =[0 L2+L4*(cos(tht3)*sin(tht4)+cos(tht4)*sin(tht3))+L3*sin(tht3) L4*(sin(tht2)*sin(tht3)*sin(tht4)-cos(tht3)*cos(tht4)*sin(tht2))-L3*cos(tht3)*sin(tht2) -cos(tht2)*(L4*(cos(tht3)*sin(tht4)+cos(tht4)*sin(tht3))+L3*sin(tht3));
    0 L3*cos(tht2)*cos(tht3)-L4*(cos(tht2)*sin(tht3)*sin(tht4)-cos(tht2)*cos(tht3)*cos(tht4)) 0 sin(tht2)*(L4*(sin(tht2)*sin(tht3)*sin(tht4)-cos(tht3)*cos(tht4)*sin(tht2))-L3*cos(tht3)*sin(tht2))+cos(tht2)*(L4*(cos(tht2)*sin(tht3)*sin(tht4)-cos(tht2)*cos(tht3)*cos(tht4))-L3*cos(tht2)*cos(tht3));
    1 0 L3*cos(tht2)*cos(tht3)-L4*(cos(tht2)*sin(tht3)*sin(tht4)-cos(tht2)*cos(tht3)*cos(tht4)) -sin(tht2)*(L4*(cos(tht3)*sin(tht4)+cos(tht4)*sin(tht3))+L3*sin(tht3));
    0 0 0 sin(tht2);
    0 0 -1 0;
    0 1 0  -cos(tht2)];
 
