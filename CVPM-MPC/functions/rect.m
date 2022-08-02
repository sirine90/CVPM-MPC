function [x,y]=rectangle(center,width,height)

x=[center(1)-width/2,center(1)+width/2,center(1)+width/2,center(1)-width/2];
y=[center(2)+height/2,center(2)+height/2,center(2)-height/2,center(2)-height/2];
plot([x x(1)],[y y(1)],'r-');
