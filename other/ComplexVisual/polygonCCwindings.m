function [wNo,p,cc]=polygonCCwindings(p,cc,opt)
% Umlaufzahlen für Zusammenhangskomponenten des Polygons p ermitteln
% function [wNo,p,cc]=polygonCCwindings(p,cc,opt)
%
% Mit Hilfe des Polygons p und den Zusammenhangskomponenten aus cc, siehe
% auch calcConnectedComps, werden für die ZK die Umlaufzahlen ermittelt.
% Der schwierige Teil dabei ist, pro ZK einen Punkt im INNEREN der ZK zu 
% finden, der dann in calcWindingNumber zur Berechnung der Umlaufzahl 
% verwendet wird. Robuste Wahlen für Punkte im Inneren einer ZK sind Punkte 
% auf dem Skelett der ZK. 
% Dies wurde hier NICHT umgesetzt. Stattdessen wird versucht, durch die
% size(cc{k},2) vielen Rand-kurven der ZK heuristisch Punkte im Inneren zu
% finden (siehe guessInsidePoint im Quelltext von polygonCCwindings) und 
% von diesen wird die Umlaufzahl bestimmt. Falls alle diese so bestimmten 
% Umlaufzahlen für eine ZK übereinstimmen, wird dieses Ergebnis 
% zurückgeliefert.
%
% Input:
%       p                Polygon
%       cc               (optionale) ZK-Ränder, wie in calcConnectedComps
%                        beschrieben. Falls cc=[] oder cc weggelassen
%                        wurde, wird cc mittels 
%                          [p,~,intersections]=augmentPolylinePoI(...
%                              p,[],intersectPolylines(p,[]));
%                          cc=calcConnectedComps(intersections);
%                        berechnet.
%       opt              (optionale) Optionen
%                          insideDelta: Default: 10*eps
%                                       Parameter für Heuristik, um
%                                       innenliegende Punkte zu finden
%
% Ouput:
%       wNo              (1,n) Vektor, wobei wNo(k) die Umlaufzahl für
%                        die k-te ZK ist oder NaN, falls keine bestimmt
%                        werden konnte (vgl. Heuristik oben).
%       p                Polygon; falls das Eingabeargument cc nicht
%                        übergeben wurde und demnach erst berechnet werden
%                        musste, kann sich p verändert haben, da dort dann
%                        die Kreuzungspunkte eingefügt wurden, siehe
%                        augmentPolylinePoI
%       cc               übergebenen bzw. berechneten
%                        Zusammenhangskomponenten
%
% See also 
%   CALCCONNECTEDCOMPS
%   AUGMENTPOLYLINEPOI
%   CALCWINDINGNUMBER
%   COMPLEXVISUALLICENSE
%

%% Eingabe-Argumente
if (nargin<1)
    error('Not enough input arguments.');
end
if (nargin>3)
    error('Too many input arguments.');
end
if (nargin<3)
    opt=struct;
end

if (nargin<2)
    cc=[];
end

if (isempty(cc)) % cc nicht angegeben => berechnen
    [p,temp,intersections]=...
        augmentPolylinePoI(p,[],intersectPolylines(p,[])); %#ok<ASGLU>
    cc=calcConnectedComps(intersections);
end

delta=processOptions(opt,{'insideDelta'},{10*eps});

%% Speicher vorbelegen
wNo=NaN*ones(1,numel(cc));
ccsizes=cellfun('size',cc,2);
pts=(NaN+1i*NaN)*ones(1,sum(ccsizes));    % speichert innere Punkt

%% Heuristik, um Punkt(e) im Inneren einer jeden ZK zu finden
pos=1;
for k=1:numel(cc)
    vorz=sign(diff(cc{k}));      % findet auf jeder Teil-Umrandung der ZK
    vorz(vorz==0)=1;
    segment=NaN*ones(2,numel(vorz));
    for j=1:numel(vorz)
        range=cc{k}(1,j):vorz(j):cc{k}(2,j);
        [val,index]=max(abs(diff(p(range))));
        if isempty(val) || val==0
            warning('CVTB:polygonCCwindings:CurveIsPoint',...
                ['Bei ',num2str(k),'-ter ZK ist ',num2str(j),...
                 '-te Teilkurve nur ein Punkt!']);
            segment(:,j)=cc{k}(1,j)*[1;1];
        else
            segment(:,j)=range(index:index+1)';            
        end
        
    end
    %segment=repmat(floor(mean(cc{k})),2,1); % Punkte in Mitte der Kurve,
    %segment(2,:)=segment(2,:)+1;            % um dort heuristisch mit
    segment=reshape(p(segment),2,ccsizes(k)); % guessInsidePoint einen Punkt
    pts(pos:pos+ccsizes(k)-1)=...           % im Inneren zu finden
        guessInteriorPoint(segment(1,:),segment(2,:),delta);
    pos=pos+ccsizes(k);
end

%% determine winding numbers of test points and evaluate
w=calcWindingNumber(p,pts);
pos=1;
for k=1:numel(cc)
    if w(pos:pos+ccsizes(k)-1)==w(pos)
        % nur wenn ALLE heuristisch ermittelten Punkte in der dieser ZK
        % dieselbe Umlaufzahl ergeben, dann diese zur�ckliefern
        wNo(k)=w(pos);
    else
        simpleLoops = getPolygonCCpaths(p,cc);
        simpleLoop = simpleLoops{k};
        [pts,flag] = constructInteriorPoint(real(simpleLoop),imag(simpleLoop));
        if flag
            wtemp = calcWindingNumber(p,pts);
            if numel(unique(wtemp))==1
                wNo(k) = wtemp(1);
            end
        end
    end
    pos=pos+ccsizes(k);
end

end

function pin=guessInteriorPoint (A,B,delta)
% Bezeichnet A und B zwei Punkte in der komplexen Ebene, etwa
%                 B                  --<--B           B
%                 |                  |    |           |
%               * |                * |    |           *
%                 |                  |    |           |
%                 A                  -->--A           A
%                Idee              Problem (a)      Problem (b)
%
% dann wird versucht den Punkt "*" links neben der Strecke [A,B] wie folgt
% zu finden: 
% Es wird der EINHEITS-Vektor der senkrecht auf AB steht und nach links
% zeigt ermittelt und dann zum Mittelwert (A+B)/2 mit dem Faktor
% delta*(1+abs(A)) addiert.
%
% Probleme: bei allg. Kurven ist nicht garantiert, dass "*" noch innerhalb
% der ZK liegt. 
% (a) Ist delta*(1+abs(A)) zu groß, kann "*" schon längst
%     wieder außerhalb der ZK liegen. 
% (b) Ist delta*(1+abs(A)) zu klein, so kann (durch die beschränkte 
%     Mantissenlänge) * AUF der Strecke [A,B] liegen.
%
% Diese Methode ist so geschrieben, dass sie für mehrere A und B auch
% funktioniert. In diesem Fall müssen A und B liegende Vektoren sein.
D=diff([A;B]);      % B-A
D=1i*D./abs(D);     % Einheitsvektor senkrecht auf B-A, nach links 
pin=mean([A;B])+delta.*D.*(1+abs(A)); % zum Streckenmittelpunkt addieren
end

function [pts, flag] = constructInteriorPoint(loopx,loopy)
% This function finds all branching points of the skeleton of the simple
% loop, whose vertices are given by loopx (x-coordinates) and loopy
% (y-coordinates).
% Input:
% -------------------------------------------------------------------------
%       loopx   |   (real) x-coordinates of the vertices of the simple loop
% -------------------------------------------------------------------------
%       loopy   |   (real) y-coordinates of the vertices of the simple loop
% -------------------------------------------------------------------------
% Output:
% -------------------------------------------------------------------------
%       pts     |   (complexified) branching points of the skeleton of the
%               |   simple loop, these are guaranteed to be interior points
%               |   of the simple loop
% -------------------------------------------------------------------------
%       flag    |   flag variable indicating whether the code went through 
%               |   (true), if not, then probably because the polygon 
%               |   contains very thin filaments
% -------------------------------------------------------------------------
% 
flag = true;
% size of the bitmap figure in which the simple loop is drawn
resolution = 400;
% % figure generation
h = figure('Position',[200 200 resolution resolution]);
hold on
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'position',[0 0 1 1],'units','normalized')
h1 = fill(loopx,loopy,'k');
set(h1,'EdgeColor','None');
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
set(gca,'Ticklength',[0 0]);
axis tight% off
ax = gca;
xaxis = get(gca,'xlim');
yaxis = get(gca,'ylim');
F = getframe(ax);
close(gcf)
% % turn frame into image
[A,~] = frame2im(F);
% % convert image to grayscale
A = rgb2gray(A);
% % this may be obsolete
[xresol,yresol] = size(A);
% % get the physical coordinates corresponding to the figure
x = linspace(min(xaxis),max(xaxis),xresol);
y = linspace(min(yaxis),max(yaxis),yresol);
[xi, yi] = meshgrid(x,y);
% % flip y-coordinates, since image and physical coordinates are upside-down
yi = flipud(yi);
% % find the black pixels (i.e. the polygon)
A = A==0;
%     figure
%     imshow(A)
try
    % % construct the skeleton
    B = bwmorph(A,'skel',Inf);
    %     figure
    %     imshow(B)
    % % extract the branchpoints of the skeleton
    C = bwmorph(B,'branchpoints');
    %     figure
    %     imshow(C)
    % % find the linear indices of the branch points
    index = find(C==1);
    pts = xi(index) + 1i*yi(index);
catch
    warning('Polygon probably too thin to build its skeleton.');
    flag = false;
    pts = 0;
end
end
