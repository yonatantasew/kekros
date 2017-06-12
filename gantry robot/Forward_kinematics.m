%--------------------------------------------------------------------------
%----------------------------Forward Kinematics----------------------------
%--------------------------------------------------------------------------
%*************************DH - Parameters**********************************
%***************   +------+------+-----------+----+----------+*************
%***************   |links | ai-1 | alpha(i-1)| di | theta(i) |*************
%***************   +------+------+-----------+----+----------+*************
%***************   |  1   |  0   |     0     | L1 |    0     |*************
%***************   +------+------+-----------+----+----------+*************
%***************   |  2   |  0   |    pi/2   | L2 |  theta2  |*************
%***************   +------+------+-----------+----+----------+*************
%***************   |  3   |  0   |    pi/2   | 0  |  theta3  |*************
%***************   +------+------+-----------+----+----------+*************
%***************   |  4   |  L3  |     0     | 0  |  theta4  |*************
%***************   +------+------+-----------+----+----------+*************
%***************   |  5   |  L4  |     0     | 0  |    0     |*************
%***************   +------+------+-----------+----+----------+*************
%%
%--------------------------------------------------------------------------
%-----------------------Transformation Matrix------------------------------
syms L1 L2 L3 L4 tht2 tht3 tht4 a b c d e f g h k px py pz theta3 theta4
T01 = [1 0 0 0;0 1 0 0;0 0 1 L1;0 0 0 1];
T12 = [cos(tht2) -sin(tht2) 0 0;0 0 -1 -L2;sin(tht2) cos(tht2) 0 0;0 0 0 1];
T23 = [cos(tht3) -sin(tht3) 0 0;0 0 -1 0;sin(tht3) cos(tht3) 0 0;0 0 0 1];
T34 = [cos(tht4) -sin(tht4) 0 L3;sin(tht4) cos(tht4) 0 0;0 0 1 0;0 0 0 1];
T45 = [1 0 0 L4;0 1 0 0;0 0 1 0;0 0 0 1];
T02 = T01*T12;
T03 = T02*T23;
T04 = T03*T34;
T05 = T04*T45



P0e = [T05(1,4);T05(2,4);T05(3,4)];
Px = T05(1,4)
Py = T05(2,4)
Pz = T05(3,4)


%%
%-----------------------Velocity Kinematics--------------------------------
%J =[0 -L3*cos(tht3)*sin(tht2)-L4*(sin(tht2)*cos(tht3)*cos(tht4)-sin(tht2)*sin(tht3)*sin(tht4)) -L3*cos(tht2)*sin(tht3)-L4*(cos(tht2)*cos(tht4)*sin(tht3)+cos(tht2)*cos(tht3)*sin(tht4)) -L4*(cos(tht2)*cos(tht3)*sin(tht4)+cos(tht2)*cos(tht4)*sin(tht3)); 0 0 -L3*cos(tht3)+L4*(sin(tht3)*sin(tht4)-cos(tht3)*cos(tht4)) L4*(sin(tht3)*sin(tht4)-cos(tht3)*cos(tht4));1 L3*cos(tht3)*cos(tht2)+L4*(cos(tht2)*cos(tht3)*cos(tht4)-cos(tht2)*sin(tht3)*sin(tht4)) -L3*sin(tht2)*sin(tht3)-L4*(sin(tht2)*sin(tht3)*cos(tht4)+sin(tht2)*cos(tht3)*sin(tht4)) -L4*(sin(tht2)*cos(tht3)*sin(tht4)+sin(tht2)*sin(tht3)*cos(tht4))]
jac=jacobian_of(Px,Py,Pz,L1,tht2,tht3,tht4)


%%
%-----------------------Differential Kinematics----------------------------
% Ve = J(theta)*theta_dot
%J(theta) = [Jp; Jo] = [Z0 Z1x(Pe-P1) Z2x(Pe-P2) Z3x(Pe-P3);0 Z1 Z2 Z3]
zero = [0;0;0];p0 = [0;0;0];pe=T05(1:3,4);Z0 = [0;0;1];
Z1=T01(1:3,3)
p1=T01(1:3,4)
Z2=T02(1:3,3)
p2=T02(1:3,4)
Z3=T03(1:3,3)
p3=T03(1:3,4)
Z4=T04(1:3,3)
p4=T04(1:3,4)
J1 = Z0
J2 = cross(Z1,(pe-p1))
J3 =cross(Z2,(pe-p2))
J4 = cross(Z3,(pe-p3))
J = [J1 J2 J3 J4;zero Z1 Z2 Z3]
%%
%------------------------Inverse Kinematics--------------------------------
%The inverse kinematics will be done using CLOSED LOOP INVERSE KINEMATICS
%ALGORITHM (CLIK): which uses Psuedo-Inverse Jacobian
%Position and orientation of the end effector will be in the form of:
%psuedo_J = simplify(J.'*(inv(J*J.')))
%Pe = [px;py;pz] 
%theta(end_effector) = theta2 + theta3 +theta4 @ a distance of L1 from the
J_psuedo = simplify(pinv(J))


