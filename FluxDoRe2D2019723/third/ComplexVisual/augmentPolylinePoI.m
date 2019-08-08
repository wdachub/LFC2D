function [p1,p2,intersections]=augmentPolylinePoI (p1,p2,intersections)
% Streckenzugschnittpunkte einfügen und mit Zusatzinfos versehen
% function [p1,p2,intersections]=augmentPolylinePoI (p1,p2,intersections)
% 
% intersections-struct muss mindestens die Felder haben, die in 
% intersectPolylines beschrieben werden.
%
% Bei geschl. Kurven (p2=[]) werden Start- und Endpunkt als weitere 
% Überschneidung hinzugefügt. 
% Außerdem werden bei den Streckenzügen p1 und p2 alle Schnittpunkte (von
% Teilstrecken) als Ecken in die Streckenzüge eingefügt.
% Damit wird p1 und p2 größer und Schnittpunkte sind stets an Ecken.
%
% Für jeden Schnittpunkt werden ferner die Winkel approximiert, mit
% dem die Kurve in diesem Punkt ankommt ("Einfallswinkel", angle_in) bzw. 
% mit dem die Kurve diesen Punkt verlässt ("Ausfallswinkel", angle_out).
%
% Input:
%       p1               erster Streckenzug
%       p2               zweiter Streckenzug oder []
%       intersections    struct mit Daten zu Selbstüberschneidungen,
%                        siehe intersect_polygons
%
% Output:
%       p1               erster Streckenzug mit Schnittpunkten als
%                        zusätzliche Ecken
%       p2               zweiter Streckenzug mit Schnittpunkten als
%                        zusätzliche Ecken oeer []
%       intersections    struct mit Daten zu Selbstüberschneidung,
%                        ergänzt um die Ein- und Ausfallwinkel angle1_in
%                        bzw. angle1_out und angle2_in und angle2_out.
%
% See also 
%   INTERSECTPOLYLINES
%   CALCCONNECTEDCOMPS
%   COMPLEXVISUALLICENSE
%

%% in p1 und p2 Schnittpunkte einfügen
% Beim Einfügen müssen wir genau aufpassen. Es gibt folgende Fallstricke
% a) die Daten in ind1 bzw. ind2 müssen nicht sortiert sein
%    (nur bei Selbstüberschneidungen liefert intersectPolylines ein
%     sortiertes ind1)
% b) in ind1 bzw. ind2 können Indizes MEHRFACH auftauchen. Dieser Fall
%    tritt dann ein, wenn auf einer Strecke mehrere verschiedene
%    Schnittpunkte liegen.
% => Mit sortIndPts erst mal alles sortieren

[ind,pts,perm]=sortIndPtsP1(intersections.ind1,intersections.pts);
for k=1:numel(ind)
    pos=ind(k);
    p1=[p1(1:pos),pts(k),p1(pos+1:end)]; % Punkt eingefügt
    ind(k:end)=ind(k:end)+1;             % Offset bei Rest anpassen    
end
intersections.ind1(perm)=ind;

[ind,pts,perm]=sortIndPtsP2(intersections.ind2,intersections.pts);
for k=1:numel(ind)
    pos=ind(k);
    p2=[p2(1:pos),pts(k),p2(pos+1:end)];
    ind(k:end)=ind(k:end)+1;
end
intersections.ind2(perm)=ind;

%% bei geschl. Kurven: Start- und Endpunkt hinzufügen
if (isempty(p2))
    intersections.ind1=[1,intersections.ind1,numel(p1)];
    av=mean(p1([1,end]));
    p1([1,end])=av*ones(1,2);      % Durchschnitt bei beiden speichern
    intersections.pts=[av,intersections.pts,av];
end

%% Winkel berechnen
intersections.angle1_in  = NaN*ones(size(intersections.ind1));
intersections.angle1_out = intersections.angle1_in;
intersections.angle2_in  = NaN*ones(size(intersections.ind2));
intersections.angle2_out = intersections.angle2_in;
for k=1:numel(intersections.ind1)
    angles=calc_angles1(intersections.ind1(k));
    intersections.angle1_in(k) =angles(1);
    intersections.angle1_out(k)=angles(2);
end
for k=1:numel(intersections.ind2)
    angles=calc_angles2(intersections.ind2(k));
    intersections.angle2_in(k) =angles(1);
    intersections.angle2_out(k)=angles(2);
end


%% Hilfsfunktionen
% Aus Performance-Gründen, werden p1 und p2 niemals als Argumente
% übergeben, da Matlab NUR Call-By-Values kann. Deshalb gibt es
% inner-functions, die p1 und p2 aus dem Namespace der Hauptfunktion
% verwenden.
    function [ind,pts,perm]=sortIndPtsP1(ind,pts)
        % Sortiert Indizes aufsteigend; wenn in ind Werte mehrfach
        % vorkommen, dann werden die pts so sortiert, dass die
        % Schnittpunkte "chronologisch", d.h. in der Reihenfolge, wie sie
        % durch die Kurvenparametrisierung durchlaufen werden, sortiert
        % werden.
        [ind,perm]=sort(ind);         % erst mal ind und pts sortieren
        pts=pts(perm);
        % jetzt kümmern wir uns um die Indizes, die mehrfach auftreten
        D=diff([ind,NaN]);
        if (any(D==0))
            % es gibt also Strecken mit multiplen Schnittpunkten
            for sort_k=unique(ind(D==0))
                % k ist jetzt ein Wert, der mehrfach in ind auftritt
                index=find(ind==sort_k); % hier tritt k auf
                [val,P]=sort(abs(pts(index)-p1(sort_k))); %#ok<ASGLU>
                pts(index)=pts(index(P));
            end
        end
    end

    function [ind,pts,perm]=sortIndPtsP2(ind,pts)
        % Sortiert Indizes aufsteigend; wenn in ind Werte mehrfach
        % vorkommen, dann werden die pts so sortiert, dass die
        % Schnittpunkte "chronologisch", d.h. in der Reihenfolge, wie sie
        % durch die Kurvenparametrisierung durchlaufen werden, sortiert
        % werden.
        [ind,perm]=sort(ind);         % erst mal ind und pts sortieren
        pts=pts(perm);
        % jetzt kümmern wir uns um die Indizes, die mehrfach auftreten
        D=diff([ind,NaN]);
        if (any(D==0))
            % es gibt also Strecken mit multiplen Schnittpunkten
            for sort_k=unique(ind(D==0))
                % k ist jetzt ein Wert, der mehrfach in ind auftritt
                index=find(ind==sort_k); % hier tritt k auf
                [val,P]=sort(abs(pts(index)-p2(sort_k))); %#ok<ASGLU>
                pts(index)=pts(index(P));
            end
        end
    end

    function angles=calc_angles1(pos)
        angles=[NaN,NaN];
        for s = [-1,1]
            q=pos+s;
            while (q>0) && (q<=numel(p1))
                diff=p1(q)-p1(pos);
                if (abs(diff)>1e-15)
                    angles((s+3)/2)=angle(diff);
                    break;
                end
                q=q+s;
            end
        end
    end

    function angles=calc_angles2(pos)
        angles=[NaN,NaN];
        for s = [-1,1]
            q=pos+s;
            while (q>0) && (q<=numel(p2))
                diff=p2(q)-p2(pos);
                if (abs(diff)>1e-15)
                    angles((s+3)/2)=angle(diff);
                    break;
                end
                q=q+s;
            end
        end
    end
end