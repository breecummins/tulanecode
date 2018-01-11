clear
close all

s=0:0.01:1;

%closed curve

l1x = (1-s)*cos(pi/6)+s*cos(2*pi/3);
l1y = (1-s)*sin(pi/6)+s*sin(2*pi/3);
l2x = (1-s)*cos(2*pi/3)+s*cos(7*pi/6);
l2y = (1-s)*sin(2*pi/3)+s*sin(7*pi/6); 
l3x = (1-s)*cos(7*pi/6)+s*cos(5*pi/3);
l3y = (1-s)*sin(7*pi/6)+s*sin(5*pi/3);
l4x = (1-s)*cos(5*pi/3)+s*cos(pi/6);  
l4y = (1-s)*sin(5*pi/3)+s*sin(pi/6);  
set(0,'DefaultAxesFontSize',24)

figure
hold on

plot(l1x,l1y,'k-','LineWidth',2)
plot(l2x,l2y,'k-','LineWidth',2)
plot(l3x,l3y,'k-','LineWidth',2)
plot(l4x,l4y,'k-','LineWidth',2)
axis([-1.25,1.25,-1.25,1.25])
axis off

text(cos(pi/6)+0.025,sin(pi/6)+0.025,'x_0','FontSize',24)
text(cos(2*pi/3)-0.025,sin(2*pi/3)+0.09,'x_1','FontSize',24)
text(cos(7*pi/6)-0.11,sin(7*pi/6)-0.1,'x_2','FontSize',24)
text(cos(5*pi/3)+0.01,sin(5*pi/3)-0.2,'x_3','FontSize',24)

text(cos(5*pi/12)-0.2,sin(5*pi/12)-0.15,'S_0','FontSize',24)
text(cos(11*pi/12),sin(11*pi/12),'S_1','FontSize',24)
text(cos(17*pi/12),sin(17*pi/12),'S_2','FontSize',24)
text(cos(23*pi/12)-0.2,sin(23*pi/12)-0.1,'S_3','FontSize',24)

%open curve

l1x = (1-s)*0+s*0.5;
l1y = (1-s)*0+s*0.5;
l2x = (1-s)*0.5+s*1;
l2y = (1-s)*0.5+s*0; 
l3x = (1-s)*1+s*1.5;
l3y = (1-s)*0+s*0.5;
set(0,'DefaultAxesFontSize',24)

figure
hold on

plot(l1x,l1y,'k-','LineWidth',2)
plot(l2x,l2y,'k-','LineWidth',2)
plot(l3x,l3y,'k-','LineWidth',2)
axis([-0.25,1.75,-0.5,0.75])
axis off

text(0,0-0.1,'x_0','FontSize',24)
text(0.5,0.5+0.025,'x_1','FontSize',24)
text(1,0-0.1,'x_2','FontSize',24)
text(1.5,0.5+0.025,'x_3','FontSize',24)

text(0.25-0.2,0.25,'S_0','FontSize',24)
text(0.75+0.025,0.25,'S_1','FontSize',24)
text(1.25+0.1,0.25,'S_2','FontSize',24)








