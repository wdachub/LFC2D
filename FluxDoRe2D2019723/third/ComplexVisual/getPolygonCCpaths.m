function paths = getPolygonCCpaths (p,cc)
% extrahiert aus ZK-Darstellung die ZK-Randpfade
% function paths = getPolygonCCpaths (p,cc)
%
% Input:
%       p                Polygon
%       cc               Zusammenhangskomponenten, wie etwa von
%                        calcConnectedComps geliefert
%
% Output:
%       paths            cell der Größe (1,numel(cc)); jedes Cell-Element
%                        ist der Polygonzug, der den Rand einer ZK im
%                        mathematisch positiven Sinne umläuft.
%
% See also 
%   CALCCONNECTEDCOMPS
%   COMPLEXVISUALLICENSE
%

if (nargin<2)
    error('Not enough input arguments.');
end
if (nargin>2)
    error('Too many input arguments.');
end
paths=cell(1,numel(cc));

for k=1:numel(cc)
    M=cc{k};
    D=diff(M);
    vorz=sign(D);
    anz=abs(D)+1;
    boundary=zeros(1,sum(anz));              % Speicher allokieren
    pos=1;
    for j=1:numel(anz)
        boundary(pos:pos+anz(j)-1)=p( M(1,j):vorz(j):M(2,j) );
        pos=pos+anz(j);
    end
    paths{k}=boundary;
end

end