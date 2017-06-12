%-----------------------------DAYNAMICS---------------------------------*
%-----------------------Euler-Lagrangian Method----------------------------
%--------------------------------------------------------------------------
%----------------Position of centriod of each links------------------------
m1=1;m2=0.6;m3=0.9;m4=0.7;
L2=0.1;L3=0.28;L4=0.30;
Izz1 = 2.7750e-06;
Izz2 = 6.5935e-05;
Izz3 = 1.6678e-04;
Izz4 = 2.0147e-04;
g=9.8;
%L1 = 1;
%tht2 = 0;
%tht3 = 0;
%tht4 = 0;
syms L1 tht2 tht3 tht4 L1_dot tht2_dot tht3_dot tht4_dot t x1 x2% m1 m2 m3 m4
pc1 = [0;0;L1];
pc2 = [0;-L2/2;L1];
pc3 = [0;-L2-L3/2*sin(tht3);L1];
pc4 = [L3*cos(tht2)*cos(tht3)-L4/2*(cos(tht2)*sin(tht3)*sin(tht4)-cos(tht2)*cos(tht3)*cos(tht4)); -L2-L4/2*(cos(tht3)*sin(tht4)+cos(tht4)*sin(tht3))-L3*sin(tht3);L1-L4/2*(sin(tht2)*sin(tht3)*sin(tht4)-cos(tht3)*cos(tht4)*sin(tht2))+L3*cos(tht3)*sin(tht2)];
%%
%---------Jacobian for the centriod of each links linear velocity----------
%------------------------Geometrical Approach------------------------------

Jv1 = [diff(pc1(1),L1) diff(pc1(1),tht2) diff(pc1(1),tht3) diff(pc1(1),tht4);diff(pc1(2),L1) diff(pc1(2),tht2) diff(pc1(2),tht3) diff(pc1(2),tht4);diff(pc1(3),L1) diff(pc1(3),tht2) diff(pc1(3),tht3) diff(pc1(3),tht4)];

Jv2 = [diff(pc2(1),L1) diff(pc2(1),tht2) diff(pc2(1),tht3) diff(pc2(1),tht4);diff(pc2(2),L1) diff(pc2(2),tht2) diff(pc2(2),tht3) diff(pc2(2),tht4);diff(pc2(3),L1) diff(pc2(3),tht2) diff(pc2(3),tht3) diff(pc2(3),tht4)]

Jv3 = [diff(pc3(1),L1) diff(pc3(1),tht2) diff(pc3(1),tht3) diff(pc3(1),tht4);diff(pc3(2),L1) diff(pc3(2),tht2) diff(pc3(2),tht3) diff(pc3(2),tht4);diff(pc3(3),L1) diff(pc3(3),tht2) diff(pc3(3),tht3) diff(pc3(3),tht4)];

Jv4 = [diff(pc4(1),L1) diff(pc4(1),tht2) diff(pc4(1),tht3) diff(pc4(1),tht4);diff(pc4(2),L1) diff(pc4(2),tht2) diff(pc4(2),tht3) diff(pc4(2),tht4);diff(pc4(3),L1) diff(pc4(3),tht2) diff(pc4(3),tht3) diff(pc4(3),tht4)];

Jv = [Jv1 Jv2 Jv3 Jv4]
%%
%------Jacobian for the centroid of each links angular velocity------------

Jw1 = [0 0 0 0;0 0 0 0;0 0 0 0];

Jw2 = [0 0 0 0;0 0 0 0;0 1 0 0]

Jw3 = [0 0 0 0;0 0 0 0;0 1 1 0];

Jw4 = [0 0 0 0;0 0 0 0;0 1 1 1];

Jv = [Jv1 Jv2 Jv3 Jv4;Jw1 Jw2 Jw3 Jw4]
%Jv = [Jv1 Jw1;Jv2 Jw2;Jv3 Jw3;Jv4 Jw4]
%%
%------------------------singularity---------------------------------------
jacobian = [0  L4*(sin(tht2)*sin(tht3)*sin(tht4)-cos(tht3)*cos(tht4)*sin(tht2))-L3*cos(tht3)*sin(tht2)  -L4*(cos(tht2)*cos(tht3)*sin(tht4)+cos(tht2)*cos(tht4)*sin(tht3))-L3*cos(tht2)*sin(tht3)  -L4*(cos(tht2)*cos(tht3)*sin(tht4)+cos(tht2)*cos(tht4)*sin(tht3));
            0  0                                                                                        -L4*(cos(tht3)*cos(tht4)-sin(tht3)*sin(tht4))-L3*cos(tht3)                                -L4*(cos(tht3)*cos(tht4)-sin(tht3)*sin(tht4));
            1  L3*cos(tht2)*cos(tht3)-L4*(cos(tht2)*sin(tht3)*sin(tht4)-cos(tht2)*cos(tht3)*cos(tht4))  -L4*(cos(tht3)*sin(tht2)*sin(tht4)+cos(tht4)*sin(tht2)*sin(tht3))-L3*sin(tht2)*sin(tht3)  -L4*(cos(tht3)*sin(tht2)*sin(tht4) + cos(tht4)*sin(tht2)*sin(tht3));
            0  0                                                                                         0                                                                                         1                                                                 ]                                    
determinant = det(jacobian)
%%
%--------------------Matrix M----------------------------------------------
M1_l = m1*(Jv1.')*Jv1;
M1_a = Jw1.'*Izz1*Jw1;
M2_l = m2*(Jv2.')*Jv2;
M2_a = Jw2.'*Izz2*Jw2;
M3_1 = m3*(Jv3.')*Jv3;
M3_a = Jw3.'*Izz3*Jw3;
M4_1 = m4*(Jv4.')*Jv4
M4_a = Jw4.'*Izz4*Jw4
M = simplify(M1_l+M1_a+M2_l+M2_a+M3_1+M3_a+M4_1+M4_a)

%% 
%--------------------CENTRIFUGAL AND CORIOLIS EFFECT-----------------------
%vector V = C(theta)[(theta_dot)^2] + B(theta)[(theta_dot)n-1(theta_dot)n)]
%C = [b111 b122 b133 b144;b211 b222 b233 b244;b311 b322 b333 b344;b411....
%     b422 b433 b444]*[(L1_dot)^2;(tht2_dot)^2;(tht3_dot)^2;(tht4_dot)^2]
%B = [2b12 2b113 2b114;2b212 2b213 2b214; 2b312 2b313 2b314;2b412 2b413...
%     2b414]*[L1_dot*tht2_dot;L1_dot*tht3_dot;L1_dot*tht4_dot]
%where: C and B = christoffel symbols
%       bijk = 1/2(mijk + mikj - mjki), m = intertia tensor
%like: b312 = simplify(1/2*(diff(M(3,1),tht2) + diff(M(3,2),L1) - diff(M(1,2),tht3)))
%       mijk = diff(m(i,j),theta_k)
%
C = [0 -(m4*sin(tht2)*(L4*cos(tht3 + tht4) + 2*L3*cos(tht3)))/2 -(m4*sin(tht2)*(L4*cos(tht3 + tht4) + 2*L3*cos(tht3)))/2 -(L4*m4*cos(tht3 + tht4)*sin(tht2))/2;0 0 0 0;0 (m4*(4*sin(2*tht3)*L3^2 + 4*sin(2*tht3 + tht4)*L3*L4 + sin(2*tht3 + 2*tht4)*L4^2))/8 -(L3^2*m3*sin(2*tht3))/8 -(L3*L4*m4*sin(tht4))/2; 0 (L4*m4*(L4*sin(2*tht3 + 2*tht4) + 2*L3*sin(tht4) + 2*L3*sin(2*tht3 + tht4)))/8 (L3*L4*m4*sin(tht4))/2 0]

B = [0 0 0;0 0 0;0 0 0;0 0 0];

%Vector V will be
V = C*[(L1_dot)^2;(tht2_dot)^2;(tht3_dot)^2;(tht4_dot)^2] + B*[L1_dot*tht2_dot;L1_dot*tht3_dot;L1_dot*tht4_dot]
%%
%---------------------------The Gravity Vector-----------------------------
% G = -[Jv1.'*m1*g + Jv2.'*m2*g + Jv3.'*m3*g + Jv4.'*m4*g]
% where: g = [g_x g_y g_z]^T
g1 = [0 -g 0].';
g2 = [0 0 g].';
g3 = [0 g*sin(tht3) 0].';
g4 = [0 g*sin(tht4) 0].';
G = -[Jv1.'*m1*g1 + Jv2.'*m2*g2 + Jv3.'*m3*g3 + Jv4.'*m4*g4]

%%
%-----------------Dynamic Model for the PRRR Robotic arm-------------------
% Torque/Force = M(theta)*theta_double_dot + V(theta, theta_dot) + G(theta)
T = M + V + G

%%
%--------------------------potential energy--------------------------------
U = [0; L1*g*m2; -m3*g*sin(tht3)*(L2+L3/2*sin(tht3)); g*cos(tht4)*(L3*cos(tht2)*cos(tht3)-L4/2*(cos(tht2)*sin(tht3)*sin(tht4))) + g*sin(tht4)*(-L2-L4/2*(cos(tht3)*sin(tht4))-L3*sin(tht3))];
%-------------------------------lagrangian---------------------------------
%individual lagrange components
q_dot= [L1_dot;tht2_dot;tht3_dot;tht4_dot];
q = [L1;tht2;tht3;tht4];
T_1 = M(1)+M(5)+M(9)+M(13)+U(1)
T_2 = M(2)+M(6)+M(10)+M(14)+U(2);
T_3 = M(3)+M(7)+M(11)+M(15)+U(3);
T_4 = M(4)+M(8)+M(12)+M(16)+U(4);
Lagrangian = [T_1;T_2;T_3;T_4];
partial_derivative_of_lagrangian_L1_dot = diff(Lagrangian,q_dot(1));
partial_derivative_of_lagrangian_tht2_dot = diff(Lagrangian,q_dot(2));
partial_derivative_of_lagrangian_tht3_dot = diff(Lagrangian,q_dot(3));
partial_derivative_of_lagrangian_tht4_dot = diff(Lagrangian,q_dot(4));
second_derivative_L1 = diff(partial_derivative_of_lagrangian_L1_dot,t);
second_derivative_tht2 = diff(partial_derivative_of_lagrangian_tht2_dot,t);
second_derivative_tht3 = diff(partial_derivative_of_lagrangian_tht3_dot,t);
second_derivative_tht4 = diff(partial_derivative_of_lagrangian_tht4_dot,t);
partial_derivative_of_lagrangian_L1 = diff(Lagrangian,q(1));
partial_derivative_of_lagrangian_tht2 = diff(Lagrangian,q(2));
partial_derivative_of_lagrangian_tht3 = diff(Lagrangian,q(3));
partial_derivative_of_lagrangian_tht4 = diff(Lagrangian,q(4));
torque_1 = second_derivative_L1 - partial_derivative_of_lagrangian_L1
torque_2 = second_derivative_tht2 - partial_derivative_of_lagrangian_tht2
torque_3 = second_derivative_tht3 - partial_derivative_of_lagrangian_tht3
torque_4 = second_derivative_tht4 - partial_derivative_of_lagrangian_tht4
Torque = [torque_1;torque_2;torque_3;torque_4]

%%
%______________________controller__________________________________________
%-------------------------state_space equation-----------------------------
%x1 = tht;
%x2 = tht_dot;
syms L1 tht tht_dot 

 M = [16/5                                 105*tht*tht+28/1000                                                                                                        -(7*tht)*(15*(tht-tht))+(28*tht)/1000            -(21*(tht-tht)*(tht))/200                                  
     (105*tht*tht + 28)/1000               (7*(15-tht*tht + 28)^2)/100000 + (7*tht^2*(15-tht*tht + 28)^2)/100000 + 32037198302574327/73786976294838206464              3.6825e-04                                       2.0147e-04
     -(7*tht*15*(tht+tht)+28*tht)/1000     3.6825e-04                                                                                                                  (147)/2500 + (441)/25000 +  0.0710               (147)/5000 +  0.0160
     -21*(tht+tht)*tht/200                 2.0147e-04                                                                                                                  (147)/5000 + 0.0160                              0.0160];
C = [0 -(m4*sin(tht)*(L4*cos(tht + tht) + 2*L3*cos(tht)))/2 -(m4*sin(tht)*(L4*cos(tht + tht) + 2*L3*cos(tht)))/2 -(L4*m4*cos(tht + tht)*sin(tht))/2;0 0 0 0;0 (m4*(4*sin(2*tht)*L3^2 + 4*sin(2*tht + tht)*L3*L4 + sin(2*tht + 2*tht)*L4^2))/8 -(L3^2*m3*sin(2*tht))/8 -(L3*L4*m4*sin(tht))/2; 0 (L4*m4*(L4*sin(2*tht + 2*tht) + 2*L3*sin(tht) + 2*L3*sin(2*tht + tht)))/8 (L3*L4*m4*sin(tht4))/2 0]

B = [0 0 0;0 0 0;0 0 0;0 0 0];

%Vector V will be
V = C*[(tht_dot)^2;(tht_dot)^2;(tht_dot)^2;(tht_dot)^2] + B*[tht_dot*tht_dot;tht_dot*tht_dot;tht_dot*tht_dot]

G = [-147/25;
    0;
    (49*sin(tht)*((49*cos(tht))/250 + (21*cos(tht)*cos(tht))/200 - (21*sin(tht)*sin(tht))/200))/5 + (3087*cos(tht)*sin(tht))/2500;
    (49*sin(tht)*((21*cos(tht)*cos(tht))/200 - (21*sin(tht)*sin(tht))/200))/5];
 

D =-inv(M)*(V+G)
f = simplify([x2;D])
J = [diff(f(1),tht) diff(f(1),tht_dot);
    diff(f(2),tht) diff(f(2),tht_dot)]