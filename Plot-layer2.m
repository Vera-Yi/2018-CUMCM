x=0.006/90:0.006/90:0.006;
y=1:60:5400;
[X,Y]=meshgrid(x,y);
Z=-38381*X.^2+0.0020241*Y+36.991;
mesh(X,Y,Z);
xlabel('x(m)');
ylabel('t(s)');
zlabel('T(��)');
title('��II���¶ȷֲ�ͼ');