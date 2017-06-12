L1=0.3;
L2=0.5;
L3 =0.4;
P0=[];
P1=[];
P2=[];
P3=[];
for theta_1=pi:pi/4:2*pi;
    theta_2=pi:pi/2:2*pi;
    theta_3=pi:pi/2:2*pi;

    p0=dh_mat(theta_1,pi/2,0,0);
    p0xyz=[p0(1,4) p0(2,4) p0(3,4)];
    P0=[P0;p0xyz];

    p1=dh_mat(theta_2,0,0,L1);
    p1xyz=[p1(1,4) p1(2,4) p1(3,4)]
    P1=[P1;p1xyz];
    p2=dh_mat(theta_3,0,0,L2);
    p2xyz=[p2(1,4) p2(2,4) p2(3,4)]
    P2=[P2;p2xyz];
    p_e = dh_mat(theta_3,0,L3*cos(theta_3),L3*sin(theta_3));
    p3=p0*p1*p2*p_e;
    p3xyz=[p3(1,4) p3(2,4) p3(3,4)];
    P3=[P3;p3xyz];

    cla;
   line([p0(1,4) p1(1,4)], [p0(2,4) p1(2,4)], [p0(3,4) p1(3,4)] ,'color',[1 0 0],'linewidth',3);
    hold on
   line([p1(1,4) p2(1,4)], [p1(2,4) p2(2,4)], [p1(3,4) p2(3,4)] ,'color',[0 1 0],'linewidth',3);
    hold on
    line([p2(1,4) p3(1,4)], [p2(2,4) p3(2,4)], [p2(3,4) p3(3,4)] ,'color',[0 0 1],'linewidth',3);
    plot3(0,0,0,'om',p2(1,4),p2(2,4),p2(3,4),'ob',p3(1,4),p3(2,4),p3(3,4),'or')
    plot3(P3(:,1),P3(:,2),P3(:,3),'--m','linewidth',2)
    [x_grid,y_grid]=meshgrid(-1:0.1:1,-1:0.1:1);
    surf(x_grid,y_grid,ones(size(x_grid)).*0,'Edgecolor','cyan','FaceColor','m','Ambientstrength',0.25,'facealpha',.2)
    axis([-1 1 -1 1 0 1])
    grid on 
    view(3)
    pause(0.1)

   
end