%%--------------------------JOINT SPACE TRAJECTORY--------------------------
%function q = joint_space_trajectory(qi,qvi,qf,qvf,t)
%Third-order polinomials with a single via-point
% q(t) = a3*t^3 + a2*t^2 + a1*t + a0
qi = 5;
qvi = 0;
qf = 15;
qvf = 0
ts = 0.1;
ti = 0;
tf = 6;
pos = [];
q_both = [];
for t = ti:ts:tf
    a0 = qi;
    a1 = qvi; 
    a2 = 1/t^3*(3*(qf-qi)-t*(2*qvi+qvf));
    a3 = 1/t^3*(2*(qi-qf)+t*(qvi+qvf));
    q =  a3*t^3 + a2*t^2 + a1*t + a0;
    q_dot = 3*a3*t^2+2*a2*t+a1;
    pos = [pos;q];
    q_both = [q_both;q_dot];
    
    
    hold on;
    plot(pos(1),q_both(1),'-b')
end