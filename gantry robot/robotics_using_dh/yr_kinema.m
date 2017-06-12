Syms theta1 l1 theta2 l2
T1=[cos(theta1) -sin(theta1) 0 0;
    sin(theta1) cos(theta1) 0 0;
    0 0 1 l1;
    0 0 0 1];

%%
T2=[1 0 0 0;
    0 0 -1 0;
    0 1 0 0;
    0 0 0 1];


T3= [cos(theta2) -sin(theta2) 0 l2*cos(theta2);
     sin(theta2) cos(theta2) 0 l2*sin(theta2);
    0 0 1 0;
    0 0 0 1];  

T =T1*T2*T3;