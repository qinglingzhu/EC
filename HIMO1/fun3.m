% x=-5:0.1:5;
% y=-5:0.1:5;
% z=-5:0.1:5;
% f1=-10*exp(-0.2*sqrt(x.^2+y.^2))-10*exp(-0.2*sqrt(y.^2+z.^2));
% 
% f2=((abs(x)).^0.8+5*sin((x).^3))+((abs(y)).^0.8+5*sin((y).^3))+((abs(z)).^0.8+5*sin((z).^3));
% pa=[f1,f2];
x=linspace(-2, 2, 25); % 在x轴上取25点

y=linspace(-2, 2, 25); % 在y轴上取25点

[xx,yy]=meshgrid(x, y); % xx和yy都是21x21的矩阵

zz=xx.*exp(-xx.^2-yy.^2); % 计算函数值，zz也是21x21的矩阵

mesh(xx, yy, zz); % 画出立体网状图