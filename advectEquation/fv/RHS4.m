function dy= RHS4(t,phi,xg,vel,pfun)

N=length(xg);
yg=xg;
h=xg(2)-xg(1);
% dy=zeros(N-4,N-4);
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



%velocity average
V.BU=BxVel(yyg+h/2,xxg-h/2,xxg+h/2,t,vel);%average upper edge
V.BD=BxVel(yyg-h/2,xxg-h/2,xxg+h/2,t,vel);
V.BR=ByVel(xxg+h/2,yyg-h/2,yyg+h/2,t,vel);
V.BL=ByVel(xxg-h/2,yyg-h/2,yyg+h/2,t,vel);

%dirichelet

C.BL(inn(1),inn)=ByInitialFun(xxg(inn(1),inn)-h/2,yyg(inn(1),inn)-h/2,yyg(inn(1),inn)+h/2,t,vel,pfun);
C.BR(inn(end),inn)=ByInitialFun(xxg(inn(end),inn)+h/2,yyg(inn(end),inn)-h/2,yyg(inn(end),inn)+h/2,t,vel,pfun);
C.BD(inn,inn(1))=BxInitialFun(yyg(inn,inn(1))-h/2,xxg(inn,inn(1))-h/2,xxg(inn,inn(1))+h/2,t,vel,pfun);
C.BU(inn,inn(end))=BxInitialFun(yyg(inn,inn(end))+h/2,xxg(inn,inn(end))-h/2,xxg(inn,inn(end))+h/2,t,vel,pfun);
%extrapolation
%update cell average on boundary
C.VA(inn(end)+1,inn)=(-13*C.VA(inn(end),inn)+5*C.VA(inn(end)-1,inn)-C.VA(inn(end)-2,inn)+12*C.BR(inn(end),inn))/3;
C.VA(inn(end)+2,inn)=(-70*C.VA(inn(end),inn)+32*C.VA(inn(end)-1,inn)-7*C.VA(inn(end)-2,inn)+48*C.BR(inn(end),inn))/3;

C.VA(inn(1)-1,inn)=(-13*C.VA(inn(1),inn)+5*C.VA(inn(1)+1,inn)-C.VA(inn(1)+2,inn)+12*C.BL(inn(1),inn))/3;
C.VA(inn(1)-2,inn)=(-70*C.VA(inn(1),inn)+32*C.VA(inn(1)+1,inn)-7*C.VA(inn(1)+2,inn)+48*C.BL(inn(1),inn))/3;


C.VA(inn,inn(end)+1)=(-13*C.VA(inn,inn(end))+5*C.VA(inn,inn(end)-1)-C.VA(inn,inn(end)-2)+12*C.BU(inn,inn(end)))/3;
C.VA(inn,inn(end)+2)=(-70*C.VA(inn,inn(end))+32*C.VA(inn,inn(end)-1)-7*C.VA(inn,inn(end)-2)+48*C.BU(inn,inn(end)))/3;


C.VA(inn,inn(1)-1)=(-13*C.VA(inn,inn(1))+5*C.VA(inn,inn(1)+1)-C.VA(inn,inn(1)+2)+12*C.BD(inn,inn(1)))/3;
C.VA(inn,inn(1)-2)=(-70*C.VA(inn,inn(1))+32*C.VA(inn,inn(1)+1)-7*C.VA(inn,inn(1)+2)+48*C.BD(inn,inn(1)))/3;




%4 order
ex1=[inn(1)-1,inn,inn(end)+1];%index vector
inm2=inn(3:end-1);


C.BL(inn(1),ex1)=(25*C.VA(inn(1),ex1)-23*C.VA(inn(2),ex1)+13*C.VA(inn(3),ex1)-3*C.VA(inn(4),ex1))/12;
C.BL(inn(2),ex1)=(3*C.VA(inn(1),ex1)+13*C.VA(inn(2),ex1)-5*C.VA(inn(3),ex1)+1*C.VA(inn(4),ex1))/12;
C.BL(inm2,ex1)=(-C.VA(inm2-2,ex1)+7*C.VA(inm2-1,ex1)+7*C.VA(inm2,ex1)-1*C.VA(inm2+1,ex1))/12;
C.BL(inn(end),ex1)=(3*C.VA(inn(end),ex1)+13*C.VA(inn(end-1),ex1)-5*C.VA(inn(end-2),ex1)+1*C.VA(inn(end-3),ex1))/12;

C.BR(inn(1:end-1),ex1)=C.BL(inn(2:end),ex1);
C.BR(inn(end),ex1)=(25*C.VA(inn(end),ex1)-23*C.VA(inn(end-1),ex1)+13*C.VA(inn(end-2),ex1)-3*C.VA(inn(end-3),ex1))/12;


C.BD(ex1,inn(1))=(25*C.VA(ex1,inn(1))-23*C.VA(ex1,inn(2))+13*C.VA(ex1,inn(3))-3*C.VA(ex1,inn(4)))/12;
C.BD(ex1,inn(2))=(3*C.VA(ex1,inn(1))+13*C.VA(ex1,inn(2))-5*C.VA(ex1,inn(3))+1*C.VA(ex1,inn(4)))/12;
C.BD(ex1,inm2)=(-C.VA(ex1,inm2-2)+7*C.VA(ex1,inm2-1)+7*C.VA(ex1,inm2)-1*C.VA(ex1,inm2+1))/12;
C.BD(ex1,inn(end))=(3*C.VA(ex1,inn(end))+13*C.VA(ex1,inn(end-1))-5*C.VA(ex1,inn(end-2))+1*C.VA(ex1,inn(end-3)))/12;

C.BU(ex1,inn(1:end-1))=C.BD(ex1,inn(2:end));
C.BU(ex1,inn(end))=(25*C.VA(ex1,inn(end))-23*C.VA(ex1,inn(end-1))+13*C.VA(ex1,inn(end-2))-3*C.VA(ex1,inn(end-3)))/12;



%boundary condition
C.BL(inn(1),inn)=ByInitialFun(xxg(inn(1),inn)-h/2,yyg(inn(1),inn)-h/2,yyg(inn(1),inn)+h/2,t,vel,pfun);
C.BR(inn(end),inn)=ByInitialFun(xxg(inn(end),inn)+h/2,yyg(inn(end),inn)-h/2,yyg(inn(end),inn)+h/2,t,vel,pfun);
C.BD(inn,inn(1))=BxInitialFun(yyg(inn,inn(1))-h/2,xxg(inn,inn(1))-h/2,xxg(inn,inn(1))+h/2,t,vel,pfun);
C.BU(inn,inn(end))=BxInitialFun(yyg(inn,inn(end))+h/2,xxg(inn,inn(end))-h/2,xxg(inn,inn(end))+h/2,t,vel,pfun);


phiu=C.BU(inn,inn).*V.BU(inn,inn)...
    +C.BR(inn,inn).*V.BR(inn,inn)...
    -C.BD(inn,inn).*V.BD(inn,inn)...
    -C.BL(inn,inn).*V.BL(inn,inn);



gphi=(C.BU(inn+1,inn)-C.BU(inn-1,inn)).*(V.BU(inn+1,inn)-V.BU(inn-1,inn))...
    -(C.BD(inn+1,inn)-C.BD(inn-1,inn)).*(V.BD(inn+1,inn)-V.BD(inn-1,inn))...
    +(C.BR(inn,inn+1)-C.BR(inn,inn-1)).*(V.BR(inn,inn+1)-V.BR(inn,inn-1))...
    -(C.BL(inn,inn+1)-C.BL(inn,inn-1)).*(V.BL(inn,inn+1)-V.BL(inn,inn-1));

dy=(phiu+1/48*gphi)/(-h);


dy=dy(:);

end

