function dh = DH_of(a,alpha,d,theta)
%returns the transformation matrix using the approach from craig.J Book
dh = [cos(theta) -sin(theta) 0 a;
    sin(theta)*cos(alpha) cos(theta)*cos(alpha) -sin(alpha) -sin(alpha)*d;
    sin(theta)*sin(alpha) cos(theta)*sin(alpha) cos(alpha) cos(alpha)*d;
    0 0 0 1]
end

