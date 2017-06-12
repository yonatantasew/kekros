function T_oe = forwardKinematics(L1,L2,L3,L4,q2,q3,q4)
T01 = [1 0 0 0;
       0 1 0 0;
       0 0 1 L1;
       0 0 0 1];
T12 = [cos(q2) -sin(q2)  0  0;
             0        0 -1 -L2;
       sin(q2)  cos(q2)  0  0;
             0        0  0  1];
T23 = [cos(q3) -sin(q3)  0  0;
             0        0 -1  0;
       sin(q3)  cos(q3)  0  0;
             0        0  0  1];
T34 = [cos(q4) -sin(q4)  0  L3;
       sin(q4)  cos(q4)  0  0;
               0          0  1  0;
               0          0  0  1];
T45 = [1 0 0 L4;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];
T_oe = T01*T12*T23*T34*T45;
end