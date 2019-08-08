function w=calcWindingNumber (polygon,pts)
% Berechnet Umlaufzahl von polygon für Punkt(e) pts 
% function w=calcWindingNumber (polygon,pts)
%
% Input:
%       polygon          geschlossene Polyline (Streckenzug)
%
% Output:
%       pts              Punkte, bei denen die Umlaufzahl von polygon
%                        ermittelt werden soll
%
% See also 
%   COMPLEXVISUALLICENSE
%

w=zeros(size(pts));

for k=1:numel(pts)
    p=polygon-pts(k);
    
    index=find( imag(p(1:end-1))<=0 & imag(p(2:end))>0 );
    w(k)=w(k)+sum(imag(p(index).*conj(p(index+1)))<0);
    
    index=find( imag(p(1:end-1))>0 & imag(p(2:end))<=0 );
    w(k)=w(k)-sum(imag(p(index).*conj(p(index+1)))>0);
end

end