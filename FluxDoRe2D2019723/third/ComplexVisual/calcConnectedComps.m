function cc=calcConnectedComps (is)
% Zusammenhangskomponenten eines Polygonzugs ermitteln
% function cc=calcConnectedComps (is)
%
% Notwendig ist, dass is die struct mit Daten zu SELBSTÜBERSCHNEIDUNGEN
% enthält und nicht Daten über Schnittpunkte zweier unterschiedlicher
% Streckenzüge.
%
% Input:
%       is               struct mit Daten zu Selbstüberschneidungen,
%                        siehe augmentPolylinePoI
%
% Output:
%       cc               (1,n) cell, wobei n die Anzahl der
%                        Zusammenhangskomponenten (ZK) ist. Der Rand einer
%                        jeden ZK wird im math. positiven Sinn umlaufen:
%                        Dazu ist cc{k} eine Matrix mit 2 Zeilen. In jeder
%                        Spalte von cc{k} ist der Start- und End-Index
%                        im Polygonzug gespeichert, der ein Kurvenstück
%                        für den ZK-Rand bildet.
%
% See also 
%   AUGMENTPOLYLINEPOI
%   COMPLEXVISUALLICENSE
%

%% Zuordnung: ind1-Index -> pts-Index
pts=unique(is.pts);                          % Schnittpunkte 
ind2pts=NaN*ones(size(is.ind1));
for k=1:numel(ind2pts)
    ind2pts(k)=find(pts==is.pts(k),1);
end

%% Schnittpunkt-Daten befüllen
% ptsData hat so viele Einträge, wie es Schnittpunkte gibt.
% Pro Schnittpunkt werden ALLE ausgehenden Kanten notiert in
% der Form:
%     Start, Ende, Ausfall-Winkel bei Start, Einfall-Winkel bei Ende, ZK
% In ZK wird gespeichert, in welcher Zusammenhangskomponente die Kante
% schon verwendet wurde (oder 0, falls sie noch neu ist)
% ptsData{k} ist eine (m_k,5) Matrix, wobei m_k die Anzahl der ausgehenden
% Kanten ist.
ptsData=cell(1,numel(pts));

for k=1:numel(is.ind1)    
    p=ind2pts(k);
    if (k>1) % ausgehende Kante zu k-1
        ptsData{p}(end+1,:)=[k,k-1,is.angle1_in(k),is.angle1_out(k-1),0];
    end
    if (k<numel(is.ind1)) % ausgehende Kante zu k+1
        ptsData{p}(end+1,:)=[k,k+1,is.angle1_out(k),is.angle1_in(k+1),0];
    end
end

%% ConnComp ermitteln
cc={};                                   % hier kommen alle ZK-Umläufe rein
for point=1:numel(pts)
    for entry=1:size(ptsData{point},1)
        if (ptsData{point}(entry,5)==0)
            walk_cc_boundary;            % noch nicht verw. Kante => start
        end
    end
end

    function walk_cc_boundary
        % beginnt mit der Kante ptsData{point}(entry,:) und
        % läuft bis zur nächsten Kreuzung, um dann "links" 
        % weiter zu laufen, solange bis ZK "umrundet".
        cc_num=numel(cc)+1;              % Nummer der neuen ZK
        boundary=[];                     % Indizes für Rand-Kurven-Stücke
        w_p=point;index=entry;           % aktueller Punkt und Eintrag
        while (~isempty(index))
            edge=ptsData{w_p}(index,:);  % um diese Kante geht es gerade
            boundary(:,end+1)=edge(1:2)'; %#ok<AGROW> % speichern
            ptsData{w_p}(index,5)=cc_num;% Flag: Kante für akt. ZK verw.
                        
            w_p=ind2pts(edge(2));        % da endet aktuelle Kante
            index=find_next_boundary(... % schauen, wo es weiter geht
                edge(2:-1:1),edge(4),w_p); % oder index==[] => stop
        end
        % Speichern Kurven-Indizes für ZK-Umlauf:
        cc{end+1}=reshape(is.ind1(boundary),2,size(boundary,2));
    end

    function next_b = find_next_boundary (exclude,alpha,p)
        % ermittelt in ptsData{p}, bei welchem Eintrag es weiter geht.
        % exclude ist die Kante mit der wir gerade zu p gekommen sind, also
        % muss diese ausgenommen werden, da wir nicht zurücklaufen wollen.
        % alpha ist der Winkel mit der wir gerade in p "angekommen" sind.
        % Gesucht ist also die zu exclude nächstgelegene linke Kante, die
        % von p wegführt.
        % next_b ist der Index dieser Kante oder [], falls es keine
        % weitere (unverbrauchte) Kante gibt.
        next_b=[];
        % Winkel-Abstand zu allen wegführenden Kanten
        delta=ptsData{p}(:,3)+ ...
            2*pi*ceil( (alpha+eps-ptsData{p}(:,3))/(2*pi) )-alpha;
        % exclude-Kante wird "ausgeblendet":
        delta(ptsData{p}(:,1)==exclude(1) & ...
              ptsData{p}(:,2)==exclude(2))=-inf;
        % Jetzt Maximum der Abstände
        % "linkeste" Kante wird beim max(delta) angenommen.
        [val,index]=max(delta);
        if isfinite(val) && (ptsData{p}(index,5)==0)
            next_b=index;       % Kante gefunden und Kante noch NEU
        end
    end

end

% function debug_ptsData (is,pts,ptsData)
% % Zum Anzeigen des aktuellen Kanten-Zwischenstands
% for k=1:numel(pts)
%     fprintf('Punkt %i: %10.2f +1i*%10.2f\n',k,real(pts(k)),imag(pts(k)));    
%     for j=1:size(ptsData{k},1)
%         edge=ptsData{k}(j,:);
%         c=' ';
%         if (edge(5)~=0), c='X';end
%         fprintf('    %c(%2i -> %2i) [%3i -> %3i]:  %6.1f° -> %6.1f°\n',...
%             c,edge(1:2),is.ind1(edge(1:2)),edge(3:4)*180/pi);
%     end  
% end
% fprintf('===\n');
% 
% end