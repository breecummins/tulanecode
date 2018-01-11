clear
close all

s=0:0.01:1;

%open curve

l1x = (1-s)*0+s*0.25;
l1y = (1-s)*0+s*0.25;
l2x = (1-s)*0.25+s*0.25*(1+sqrt(2));
l2y = 0.25*ones(size(s)); 
set(0,'DefaultAxesFontSize',24)

figure
hold on

plot(l1x,l1y,'k-','LineWidth',2)
plot(l2x,l2y,'k-','LineWidth',2)
plot([0,0.25,0.25*(1+sqrt(2)),0.4],[0,0.25,0.25,0.4],'k.','MarkerSize',24)
axis([-0.05,0.75,-0.05,0.5])
axis off

text(0,0-0.05,'x_{k-1}','FontSize',24)
text(0.25,0.25+0.03,'x_k','FontSize',24)
text(0.25*(1+sqrt(2))+0.02,0.25,'x_{k+1}','FontSize',24)
text(0.4,0.4+0.025,'x_j','FontSize',24)

text(0.25-0.25,0.18,'S_{k-1}','FontSize',24)
text(0.4,0.18,'S_k','FontSize',24)

set(gcf,'PaperPosition',[0,0,8.0,6.0])
set(gcf,'PaperSize',[8.0,6.0])








