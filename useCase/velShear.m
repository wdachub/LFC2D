function vel = velShear(t, xy, T)


x = xy(:,1);
y = xy(:,2);

vel = [(y-1); zeros(size(y))];
end