function [cspace,K]=build(lengths)
% Build configuration space of the robot
% lengths is the vector of lengths of the joints of the robot,
% for example: lengths=[4,4] corresponds to a robot with 2 Dof of length 4

close all
figure(1)
hold on
axis([-20 20 -20 20]);

%Obstacle bounds
[x,y]=rect([3 8],8,8);

%
q=[160 250]; % [q1 q2]
Rob(1,:)=[0,0];
qq=0;
Q1=[];
Q2=[];
for i=1:2
  qq=qq+q(i);
  Rob(i+1,:)=[Rob(i,1)+lengths(i)*cos(qq*pi/180) ...
                Rob(i,2)+lengths(i)*sin(qq*pi/180)];
end
plot(Rob(:,1),Rob(:,2),'b-');
hold off


figure(2)
cspace=zeros(360,360);
for q1=0:.5:359
  for q2=0:.5:359
    %current joint positions
    Rob(2,:)=[Rob(1,1)+lengths(1)*cos(q1*pi/180) ...
                Rob(1,2)+lengths(1)*sin(q1*pi/180)];
    Rob(3,:)=[Rob(2,1)+lengths(2)*cos((q1+q2)*pi/180) ...
                Rob(2,2)+lengths(2)*sin((q1+q2)*pi/180)];
    
    % Check if end-effector collides
	if ((discretize(Rob(3,1), [x(1),x(2)])==1)==1) && ((discretize(Rob(3,2), [y(3),y(2)])==1) ==1)
        cspace(fix(q1+1),fix(q2+1))=1;
        Q1=[Q1; [q1 q2]];
    end
    
  end
end


%Build C-Obstacle
imshow((1-(cspace'))*255,'InitialMagnification','fit')
set(gca,'YDir','normal') 
h = gca;
h.Visible = 'On';
hold on
P1=Polyhedron(Q1);
P1.plot('color','black','alpha',0.3);
P1.minVRep;

%% c-space constraints
% Build the collision-free probabilistic set Pxp
Ap=[1,0;
    -1,0;
    0,1;
    0,-1];
bp=[360;0;360;0];
P = Polyhedron('A',Ap,'b',bp);
S=P\P1;
max=volume(S(1));
k=1;
for i=1:1:size(S)
    if volume(S(i))>max
        max=volume(S(i));
        k=i;
    end
end
K=S(k);
hold off




      
