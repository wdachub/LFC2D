function dy = velocityFieldComp(t,Pos)

omega= 2*pi;
half = numel(Pos)/2;
x = Pos(1:half); y = Pos(half+1:end);
dy = [x+omega*y; -omega*x+y];

end