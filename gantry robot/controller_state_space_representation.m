%--------------------CONTROLLER_STATE_SPACE_REPRESENTATIOM-----------------
%--------------------------state space components--------------------------
m1=1;m2=0.6;m3=0.9;m4=0.7;
L2=0.1;L3=0.28;L4=0.30;
L1 = 1;
tht2 = 0;
tht3 = 0;
tht4 = 0;
M = [ 16/5                                                           (7*cos(tht2)*(15*cos(tht3 + tht4) + 28*cos(tht3)))/1000                                                                                                     -(7*sin(tht2)*(15*sin(tht3+tht4)+28*sin(tht3)))/1000                                                 -(21*sin(tht3+tht4)*sin(tht2))/200;
    (7*cos(tht2)*(15*cos(tht3+tht4)+28*cos(tht3)))/1000            (7*cos(tht2)^2*(15*cos(tht3+tht4)+28*cos(tht3))^2)/100000+(7*sin(tht2)^2*(15*cos(tht3+tht4)+28*cos(tht3))^2)/100000+32037198302574327/73786976294838206464  13586027010287085/36893488147419103232                                                               7432931057060527/36893488147419103232;
    -(7*sin(tht2)*(15*sin(tht3 + tht4) + 28*sin(tht3)))/1000        13586027010287085/36893488147419103232                                                                                                                      (147*cos(tht4))/2500 + (441*cos(tht3)^2)/25000 + 8185540921445307332113/115292150460684697600000     (147*cos(tht4))/5000 + 73563171172363925363/4611686018427387904000;
    -(21*sin(tht3 + tht4)*sin(tht2))/200                            7432931057060527/36893488147419103232                                                                                                                       (147*cos(tht4))/5000 + 73563171172363925363/4611686018427387904000                                   73563171172363925363/4611686018427387904000];
%finally
M = [3.2000    0.3010         0         0;
     0.3010    0.1299    0.0004    0.0002;
          0    0.0004    0.1474    0.0454;
          0    0.0002    0.0454    0.0160]
 %%
 syms L1 tht2 tht3 tht4 L1_dot tht2_dot tht3_dot tht4_dot t% m1 m2 m3 m4

 M = [16/5                                     105*tht3*tht4+28/1000                                                                                                        -(7*tht2)*(15*(tht3-tht4))+(28*tht3)/1000            -(21*(tht3-tht4)*(tht2))/200                                  
     (105*tht3*tht4 + 28)/1000               (7*(15-tht3*tht4 + 28)^2)/100000 + (7*tht2^2*(15-tht3*tht4 + 28)^2)/100000 + 32037198302574327/73786976294838206464                3.6825e-04                                       2.0147e-04
      -(7*tht2*15*(tht3+tht4)+28*tht3)/1000     3.6825e-04                                                                                                                  (147)/2500 + (441)/25000 +  0.0710                  (147)/5000 +  0.0160
      -21*(tht3+tht4)*tht2/200                  2.0147e-04                                                                                                                  (147)/5000 + 0.0160                                  0.0160] 
 %           
 %%
 %--------------------------coriolis and centrifugal-----------------------
 %coriolis term
%C = [0 -(m4*sin(tht2)*(L4*cos(tht3 + tht4) + 2*L3*cos(tht3)))/2 -(m4*sin(tht2)*(L4*cos(tht3 + tht4) + 2*L3*cos(tht3)))/2 -(L4*m4*cos(tht3 + tht4)*sin(tht2))/2;0 0 0 0;0 (m4*(4*sin(2*tht3)*L3^2 + 4*sin(2*tht3 + tht4)*L3*L4 + sin(2*tht3 + 2*tht4)*L4^2))/8 -(L3^2*m3*sin(2*tht3))/8 -(L3*L4*m4*sin(tht4))/2; 0 (L4*m4*(L4*sin(2*tht3 + 2*tht4) + 2*L3*sin(tht4) + 2*L3*sin(2*tht3 + tht4)))/8 (L3*L4*m4*sin(tht4))/2 0];
%finally
C = zeros(4);
%centrifugal term
B = [0 0 0;
    0 0 0;
    0 0 0;];
%finally, the V term becomes
%V = C*[(L1_dot)^2;(tht2_dot)^2;(tht3_dot)^2;(tht4_dot)^2] + B*[L1_dot*tht2_dot;L1_dot*tht3_dot;L1_dot*tht4_dot];
%which results
V = [0 0 0 0]'
%%
%-------------------------the Gravity Term--------------------------
G =  [                                                                                                                          -147/25;
                                                                                                                                     0;
 (49*sin(tht4)*((49*cos(tht3))/250 + (21*cos(tht3)*cos(tht4))/200 - (21*sin(tht3)*sin(tht4))/200))/5 + (3087*cos(tht3)*sin(tht3))/2500;
                                                        (49*sin(tht4)*((21*cos(tht3)*cos(tht4))/200 - (21*sin(tht3)*sin(tht4))/200))/5]
 %which results
 G = [-5.8800 0 0 0]'
%%

%a = [-m_inv(1)*(V(1)+G(1)) -m_inv(1)*(V(1)+G(1)) -m_inv(1)*(V(1)+G(1)) -m_inv(1)*(V(1)+G(1))] 

m_inv =inv(M);
f = -m_inv*(V+G)

J1 = [diff(f,L1) diff(f,tht2) diff(f,tht3) diff(f,tht4)]
J2 = [diff(f,L1_dot) diff(f,tht2_dot) diff(f,tht3_dot) diff(f,tht4_dot)]
A = [zeros(4) eye(4);J1 J2]

A_new = A
%%
B = [zeros(4) eye(4)]'
%%
%for the linear equation
%V1 = ([(diff(V,L1))' (diff(V,tht2))' (diff(V,tht3))' (diff(V,tht4))'])
V1 = [0 0 0 0]';
%V2 = ([(diff(V,L1))' (diff(V,tht2))' (diff(V,tht3))' (diff(V,tht4))'])
V2 = [0 0 0 0]';
%m_1 = ([(diff(M(1,1:4),L1))' (diff(M(2,1:4),tht2))' (diff(M(3,1:4),tht3))' (diff(M(4,1:4),tht4))'])
m_1 = zeros(4);
W = G %+ V1 + m_1;
a11 = zeros(4);
a12 = eye(4);
a21 = -(inv(M)*W)'
a22 = -(inv(M)*V2)';
A = [a11 a12;a21 a22]
B = [zeros(4) inv(M)]'
