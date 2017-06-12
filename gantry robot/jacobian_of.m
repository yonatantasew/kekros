function J = jacobian_of(px,py,pz,L1,tht2,tht3,tht4)
J = [diff(px,L1) diff(px,tht2) diff(px,tht3) diff(px,tht4);
    diff(py,L1) diff(py,tht2) diff(py,tht3) diff(py,tht4);
    diff(pz,L1) diff(pz,tht2) diff(pz,tht3) diff(pz,tht4)]
