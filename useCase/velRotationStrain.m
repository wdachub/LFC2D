function vel = velRotationStrain(time, xy, T)
omega= 2*pi;
x = xy(:,1);
y = xy(:,2);
%move the center to (1,1) in order to perfrom the condition number test. 
%
x=x-1;
y=y-1;
vel = [x+omega*y; -omega*x+y];