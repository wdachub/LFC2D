function edl = edgeLength(pts)
%input the coordinate (x,y) of each point. N*2 matrix
%output: the distance between each pair of adjacent point. (n-1)*1 matrix
% determine the length of edges of adjacent points
tmp= diff(pts,1,1);
edl =sqrt(sum(tmp.^2,2));
end
