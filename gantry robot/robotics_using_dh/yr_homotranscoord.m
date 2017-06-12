function homo_trans=htm(rot_mtrx,px,py,pz)
homo_trans=[rot_mtrx(1,1) rot_mtrx(1,2) rot_mtrx(1,3) px;
            rot_mtrx(2,1) rot_mtrx(2,2) rot_mtrx(2,3) py;
            rot_mtrx(3,1) rot_mtrx(3,2) rot_mtrx(3,3) pz;
            0 0 0 1 ];
end