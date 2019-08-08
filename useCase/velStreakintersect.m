function vel = velStreakintersect(time, xy,T)

x = xy(:,1);
y = xy(:,2);

vel = [-y x-time+5];

end
