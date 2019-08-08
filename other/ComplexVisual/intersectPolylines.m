function intersections=intersectPolylines (p1,p2,opt)
% Streckenzüge p1 und p2 schneiden oder Selbstüberschneidungen von p1 finden
% function intersections=intersectPolylines (p1,p2,opt)
%
% Falls p2=[], so wird bei p1 nach Selbstüberschneidungen von p1 gesucht.
% Sollte in diesem Fall p1 geschlossen sein, so wird der Anfangs- und 
% Endpunkt NICHT als Selbstüberschneidung gezählt/gefunden.
%
% Beispiel 1: Schnitt zweier Streckenzüge
%        D                 Ist p1=[A,B,C] und p2=[D,E,F] mit den Punkten 
%          \               A-F in der komplexen Ebene, wie links in der 
%           E              Skizze zu sehen. 
%           |              Dann liefert is=intersectPolylines(p1,p2):
%     A-----+-----B        is.ind1=[1] und is.ind2=[2] und is.pts=[X], 
%           |      \       wobei X die Schnittstelle "+" ist.
%           F       \      Das bedeutet, dass sich die Strecke
%                    C     p1(is.ind1(1):is.ind1(1)+1) (also [A,B]) 
%                          am Punkt is.pts(1) (also X) mit der Strecke
%                          p2(is.ind2(1):is.ind2(1)+1) (also [E,F])
%   schneidet. Bei mehreren Strecken, die sich schneiden, kann man anhand 
%   von is.pts ermitteln, welche der Strecken sich im gleichen Schnittpunkt 
%   treffen.
%
% Beispiel 2: Selbstüberschneidungen
%        D                 Ist p1=[A,B,C,F,E,D] und p2=[] mit den Punkten 
%          \               A-F in der komplexen Ebene, wie links in der 
%           E              Skizze zu sehen. 
%           |              Dann liefert is=intersectPolylines(p1,p2):
%     A-----+-----B        is.ind1=[1,4] und is.ind2=[] und is.pts=[X,X], 
%           |      \       wobei X die Schnittstelle "+" ist.
%           F--     \      Das bedeutet, dass sich die Strecke
%              \ --- C     p1(is.ind1(1):is.ind1(1)+1) (also [A,B]) 
%                          am Punkt is.pts(1) (also X) mit der Strecke
%                          p2(is.ind1(2):is.ind1(1)+1) (also [F,E])
%   schneidet. Bei mehreren Strecken, die sich schneiden, kann man anhand 
%   von is.pts ermitteln, welche der Strecken sich im gleichen Schnittpunkt 
%   treffen.
%
% Vorsicht: 
% Diese Methode verwendet einen "BRUTE-FORCE" subdivison Algorithmus.
%
% Für schnellere Algorithmen, siehe etwa
% B. Chazelle,  H. Edelsbrunner: 
%    An Optimal Algorithm for Intersecting Line Segments in the Plane.
%    J. ACM 39, 1-54 (1992)
%                   
% Input:
%       p1          erster Streckenzug (Polyline)
%       p2          zweiter Streckenzug (oder [])
%       opt         (optionale) Optionen
%                     sec1: Indizes in p1 für die erste Unterteilung
%                           Default: in 10-er Schritten von 1 bis
%                           numel(p1)
%                     sec2: Indizes in p2 für die erste Unterteilung
%                           Default: in 10-er Schritten von 1 bis
%                           numel(p2)
%                           Falls p2=[] und sec2 verschieden von sec1
%                           ist, ist das Verhalten undefiniert.
%                     rect_thres:
%                           Threshold: Wenn sich zwei Bounding-Boxes 
%                           horizontal bzw. vertikal um mehr als rect_thres
%                           überlappen, dann wird es als Schnitt gezählt
%                           Default: 0
%
% Output:
%       intersections       Struct mit folgenden Feldern
%                    ind1:  Teilstrecken mit Schnittpunkten (darauf)
%                           für k=1:numel(ind1) gilt
%                             p1(ind1(k):ind1(k)+1) schneidet
%                             im Punkt pts(k)
%                           eine andere Teilstrecke, deren Index natürlich
%                           auch in ind2 (bzw. ind1, wenn p2=[]) vorkommt.
%                           Im Fall p2=[] ist ind1 aufsteigend sortiert.
%                           Für den Fall, dass auf einer Strecke mehrere
%                           Schnittpunkte liegen, können in ind1 Indizes
%                           mehrfach auftreten.
%                    ind2:  leer, falls p2=[]
%                           Teilstrecken mit Schnittpunkten (darauf)
%                           für k=1:numel(ind2) gilt
%                             p1(ind2(k):ind2(k)+1) schneidet
%                             im Punkt pts(k)
%                           eine andere Teilstrecke.
%                           Für den Fall, dass auf einer Strecke mehrere
%                           Schnittpunkte liegen, können in ind2 Indizes
%                           mehrfach auftreten.
%                     pts:  Schnittpunkte
%
% See also 
%   AUGMENTPOLYLINEPOI
%   COMPLEXVISUALLICENSE
%

%% Eingabe-Argumente
if (nargin<2)
    error('Not enough input arguments.');
end
if (nargin>3)
    error('Too many input arguments.');
end
if (nargin<3), opt=struct();end

ind1=[];ind2=[];            % gefundene Teilstrecken, die sich schneiden
pts=[];                     % mit Schnittpunkten

self=0;                     % Flag, ob p2=[], dann ...
if (isempty(p2))
    p2=p1;self=1;           % p1-Selbstüberschneidungen gesucht
end

% SHARED-Variablen (aus Performance, da Matlab ja nur call-by-value kann)
bb1=[];bb2=[];              % bounding-boxes für p1- und p2-Sektionen
x1=real(p1);y1=imag(p1);    % x/y-Koordinaten für p1
x2=real(p2);y2=imag(p2);    % x/y-Koordinaten für p2

[sec1,sec2,rect_thres]=processOptions(opt,...
    {'sections1','sections2','rect_thres'},{[],[],0});
% sec_list: speichert Liste der noch zu vergleichenden Sektionen
sec_list=struct('sec1',sec1,'sec2',sec2,... 
    'start1',1,'end1',numel(p1),'start2',1,'end2',numel(p2));

%% Subdivison-Schleife
warn_state=warning('query','MATLAB:singularMatrix');
warning('off','MATLAB:singularMatrix');
firstLoop=1;
while (numel(sec_list)>0)
    entry=sec_list(end);      % aktuelle Subdivisonen/Sektionen
    sec1=entry.sec1;sec2=entry.sec2;
    if (isempty(sec1))        % nur Start/Ende bekannt => Sektionen bilden
        sec1=build_index_sections(entry.start1,entry.end1);
    end
    if (isempty(sec2))        % nur Start/Ende bekannt => Sektionen bilden
        sec2=build_index_sections(entry.start2,entry.end2);
    end
    sec_list(end)=[];         % aktuellen Eintrag entfernen
    
    % Sektionen (Bounding-Boxes) finden, die sich überschneiden
    % (Vorsicht bei ersten Durchlauf und gesuchter Selbstüberschneidung)
    [range1,range2]=find_intersecting_boxes(firstLoop && self,rect_thres);

    for k=1:numel(range1)     % alle Überschneidungen näher anschauen
        b1=bb1(range1(k));    % Bounding-Box b1 von p1 überlappt mit
        b2=bb2(range2(k));    % Bounding-Box b2 von p2
        if (diff(b1.range)==1) && (diff(b2.range)==1)
            % b1 UND b2 bestehen JEWEILS nur aus einer Strecke 
            % Vorsicht bei Selbstüberschneidung von direkt angrenzenden 
            % Streckenzüge:
            if (~self) || (abs(diff([b1.range(1),b2.range(1)]))>1)
                % Test, ob Schnitt vorliegt:
                fac=[ -diff(x1(b1.range)),diff(x2(b2.range));
                      -diff(y1(b1.range)),diff(y2(b2.range))] \ ...
                      [ x2(b2.range(2))-x1(b1.range(2));
                        y2(b2.range(2))-y1(b1.range(2))];
                if (sum( 0.0<=fac & fac<=1.0 )==2)
                    % Schnitt liegt vor
                    pts(end+1)=...
                        fac(1)*(x1(b1.range(1))+1i*y1(b1.range(1)))+...
                        (1-fac(1))*(x1(b1.range(2))+1i*y1(b1.range(2))); %#ok<AGROW>
                    ind1(end+1)=b1.range(1); %#ok<AGROW> % => speichern 
                    ind2(end+1)=b2.range(1); %#ok<AGROW> % => Sektion fertig
                end
            end
        else
            % b1 und/oder b2 beinhaltet mehr Teilstrecken 
            % => b1 und/oder b2 müssen weiter unterteilt werden
            % => neuer Eintrag in sec_list
            sec_list(end+1)=struct('sec1',[],'sec2',[],...
                'start1',b1.range(1),'end1',b1.range(2),...
                'start2',b2.range(1),'end2',b2.range(2)); %#ok<AGROW>
        end
    end
    firstLoop=0;
end
warning(warn_state);

if (self)
    [ind1,P]=sort([ind1,ind2]);
    pts=[pts,pts];pts=pts(P);
    ind2=[];
end
intersections=struct('ind1',ind1,'ind2',ind2,'pts',pts);

%% Hilfs-Funktionen mit SHARED Variablen
    function [f_r1,f_r2]=find_intersecting_boxes(sflag,thres)
        % FIND_INTERSECTING_BOXES: Bounding-Boxes generieren und schneiden
        % Input:
        %       sflag       Flag, ob Selbst-Überschneidungs-Test
        %       thres       Überlapp-Threshhold für Rechtecke
        % Output:
        %       f_r1,f_r2   es gilt numel(f_r1)==numel(f_r2)
        %                   für k=1:numel(f_r1) gilt:
        %                     bb1(f_r1(k)) überlappt mit bb2(f_r2(k))       
        f_r1=[];f_r2=[];          % speichert gefundene Überlappungen
        
        build_bounding_boxes1;    % Bounding-Boxes zunächst alle ...
        build_bounding_boxes2;    % erzeugen

        for f_k=1:numel(bb1)
            f_jstart=1;
            if (sflag), f_jstart=f_k+1;end % bei sflag, nur Test mit Rest
            for f_j=f_jstart:numel(bb2)                    
                if (do_rect_intersect(bb1(f_k),bb2(f_j),thres))              
                    f_r1(end+1)=f_k; %#ok<AGROW>                        
                    f_r2(end+1)=f_j; %#ok<AGROW>                    
                end                
            end            
        end
    end

    function build_bounding_boxes1
        % für sec1-Sektionen Bounding-Boxes ermitteln
        ec=cell(1,numel(sec1)-1); % pre-alloc mem
        bb1=struct('xmin',ec,'xmax',ec,'ymin',ec,'ymax',ec,'range',ec);
        for b_i=1:numel(sec1)-1
            bb1(b_i).xmin=min(x1(sec1(b_i):sec1(b_i+1)));
            bb1(b_i).xmax=max(x1(sec1(b_i):sec1(b_i+1)));
            bb1(b_i).ymin=min(y1(sec1(b_i):sec1(b_i+1)));
            bb1(b_i).ymax=max(y1(sec1(b_i):sec1(b_i+1)));
            bb1(b_i).range=[sec1(b_i),sec1(b_i+1)];
        end
    end

    function build_bounding_boxes2
        % für sec2-Sektionen Bounding-Boxes ermitteln
        ec=cell(1,numel(sec2)-1); % pre-alloc mem
        bb2=struct('xmin',ec,'xmax',ec,'ymin',ec,'ymax',ec,'range',ec);
        for b_i=1:numel(sec2)-1
            bb2(b_i).xmin=min(x2(sec2(b_i):sec2(b_i+1)));
            bb2(b_i).xmax=max(x2(sec2(b_i):sec2(b_i+1)));
            bb2(b_i).ymin=min(y2(sec2(b_i):sec2(b_i+1)));
            bb2(b_i).ymax=max(y2(sec2(b_i):sec2(b_i+1)));
            bb2(b_i).range=[sec2(b_i),sec2(b_i+1)];
        end
    end

end

%% Hilfsfunktionen

function sections=build_index_sections (first,last)
% BUILD_INDEX_SECTIONS: von first bis last in 10er-Schritten
anz=last-first;
blocks = 20;
sections=first:ceil(anz/blocks):last;
if (sections(end)~=last)
    sections(end+1)=last;
end
end

function flag=do_rect_intersect (r1,r2,t)
% DO_RECT_INTERSECT: Test, ob sich Rechtecke (Bounding-Boxes) überschneiden
% dabei muss mindestens der horizontal und/oder vertikal ein Überlapp von t
% vorhanden sein, um als Schnitt zu zählen.
flag=~(  (r2.xmin-r1.xmax>=-t) || (r2.xmax-r1.xmin<=t) || ...
         (r2.ymin-r1.ymax>=-t) || (r2.ymax-r1.ymin<=t) );
end
