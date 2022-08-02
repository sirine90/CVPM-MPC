function Rob=createRobot()

L(1) = Link([ 0   0   4  0], 'standard');
L(1).m=15;
L(1).r=[4 0 0];
L(1).I=[0 0 0;
       0 0 0; 
       0 0 0.5 ] ;
L(1).Jm=0;
L(1).G=1; 

L(2) = Link([ 0   0   4  0], 'standard');
L(2).m=15;
L(2).r=[-4 0 0];
L(2).I=[0 0 0;
       0 0 0; 
       0 0 0.5 ] ;
L(2).Jm=0;
L(2).G=1; 

Rob= SerialLink(L,'name', 'robot');

Rob.gravity=[0;9.81;0];