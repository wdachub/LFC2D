function dy= RHS2(t,phi,xg,vel,pfun)
N=length(xg);
yg=xg;
h=xg(2)-xg(1);

% we assume the partition on x and y are same.
phi=reshape(phi,sqrt(numel(phi)),sqrt(numel(phi)));

[xxg,yyg]=meshgrid(xg,yg);
xxg=xxg';
yyg=yyg';
% we don't know the time derivative of face average, so we compute the face
% average from cell average whenever we compute the time derivative of cell
% average.

C.VA=zeros(N);
C.BL=zeros(N);
C.BR=zeros(N);
C.BU=zeros(N);
C.BD=zeros(N);

inn=3:N-2;
C.VA(inn,inn)=phi;


%update cell average on boundary
%2 order 

C.BL(inn,inn)=(C.VA(inn,inn)+C.VA(inn-1,inn))/2;
C.BR(inn(1:end-1),inn)=C.BL(inn(2:end),inn);
C.BD(inn,inn)=(C.VA(inn,inn)+C.VA(inn,inn-1))/2;
C.BU(inn,inn(1:end-1))=C.BD(inn,inn(2:end));


C.BL(inn(1),inn)=ByInitialFun(xxg(inn(1),inn)-h/2,yyg(inn(1),inn)-h/2,yyg(inn(1),inn)+h/2,t,vel,pfun);
C.BR(inn(end),inn)=ByInitialFun(xxg(inn(end),inn)+h/2,yyg(inn(end),inn)-h/2,yyg(inn(end),inn)+h/2,t,vel,pfun);
C.BD(inn,inn(1))=BxInitialFun(yyg(inn,inn(1))-h/2,xxg(inn,inn(1))-h/2,xxg(inn,inn(1))+h/2,t,vel,pfun);
C.BU(inn,inn(end))=BxInitialFun(yyg(inn,inn(end))+h/2,xxg(inn,inn(end))-h/2,xxg(inn,inn(end))+h/2,t,vel,pfun);


%velocity average
V.BU=BxVel(yyg+h/2,xxg-h/2,xxg+h/2,t,vel);%average upper edge
V.BD=BxVel(yyg-h/2,xxg-h/2,xxg+h/2,t,vel);
V.BR=ByVel(xxg+h/2,yyg-h/2,yyg+h/2,t,vel);
V.BL=ByVel(xxg-h/2,yyg-h/2,yyg+h/2,t,vel);


% C = InitialData(pfun,xg,yg,h,t);
phiu=C.BU(inn,inn).*V.BU(inn,inn)+C.BR(inn,inn).*V.BR(inn,inn)-C.BD(inn,inn).*V.BD(inn,inn)-C.BL(inn,inn).*V.BL(inn,inn);
dy=(phiu)/(-h);


dy=dy(:);




end

