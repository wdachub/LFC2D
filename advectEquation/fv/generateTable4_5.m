% The purpose of this file is comparing the computation time of the FV method and LFC method.
%It will reproduce Table 4.5 in the paper
clear

addpath '../../LFC';
addpath '../../useCase';
addpath '../../other';
addpath '../useCase';%includes the exact formula for fv numerical test.



vel=@velIntersection;
pfun=@pfunIntersection;

Bound0=4;
Refine=2.^(2:4);
Ns0=32;
Nt0=10;
t0 = 0;
te=0.1;
%
p1=[1,1];
p2=[1,2];
curve = [p1;p2];%we will compute the flux pass through this curve
exactflux=computeFlux(p1, p2, t0, te, vel, pfun);
for re=1:length(Refine)
    tic
    Bound=Bound0*Refine(re);
    Ns=Ns0*Refine(re);
    %     Nt=Nt0*Refine(re);%time step size is fixed
    Nt=Nt0;
    %     Ns=Ns0;
    h=Bound/Ns;
    x=h/2:h:Bound-h/2;%the center of the cell
    y=x;
    % y=0:h:Bound;
    xg=[x(1)-2*h,x(1)-h,x,x(end)+h,x(end)+2*h]; %add the ghost cell
    yg=xg;
    
    
    
    
    
    t=t0;
    C = InitialData(xg,yg,h,t0,vel,pfun);
    inn=3:length(xg)-2;
    [xxg,yyg]=meshgrid(xg,yg);
    xxg=xxg';
    yyg=yyg';
    
    
    
    
    phi=C.VA(inn,inn);
    phi=phi(:);
    
    time=linspace(t0,te,4);
    dt=time(2)-time(1);
    
    
    %index of the curve
    indx=find(xg==p1(1)-h/2);
    indy1=find(yg==p1(2)+h/2);
    indy2=find(yg==p2(2)-h/2);
    indy=indy1:indy2;
    exindy=[indy(1)-1,indy,indy(end)+1];
    %
    V.BR=ByVel(xxg+h/2,yyg-h/2,yyg+h/2,te,vel);
    fvflux(indy)=C.BR(indx,indy).*V.BR(indx,indy)+(C.BR(indx,indy+1)-C.BR(indx,indy-1)).*(V.BR(indx,indy+1)-V.BR(indx,indy-1))/48;
    numflux(1)=sum(fvflux);%space integral over the curve at time t.
       
    for indtime=1:length(time)-1
   
    subtime=linspace(time(indtime),time(indtime+1),Nt);
    for indsubtime=1:length(subtime)-1
         phi = RungeKutta(phi, subtime(indsubtime), subtime(indsubtime+1), @(t,phi)RHS4(t,phi,xg,vel,pfun), 4);       
    end
    phi1=reshape(phi,sqrt(numel(phi)),sqrt(numel(phi))); 
    C.VA(inn,inn)=phi1;
    C.BR(indx,exindy)=(-C.VA(indx-1,exindy)+7*C.VA(indx,exindy)+7*C.VA(indx+1,exindy)-1*C.VA(indx+2,exindy))/12;
    V.BR=ByVel(xxg+h/2,yyg-h/2,yyg+h/2,te,vel);
    fvflux(indy)=C.BR(indx,indy).*V.BR(indx,indy)+(C.BR(indx,indy+1)-C.BR(indx,indy-1)).*(V.BR(indx,indy+1)-V.BR(indx,indy-1))/48;
    numflux(indtime+1)=sum(fvflux);%space integral over the curve at time t.
    end
    numflux=(numflux(1)+3*numflux(2)+3*numflux(3)+numflux(4))/8*3*dt*h;%time integral 
    display(['The computation time of FV method is', num2str(toc)])
    display(['The relative error of FV method is', num2str(abs(numflux-exactflux)/abs(exactflux))]) 
    
    tic
    order=4;
    nseg=4;
    option=[1,1];
    numflux=fluxDR2D(curve, vel, pfun,t0,te ,nseg,order,option);
    display(['The computation time of LFC method is ', num2str(toc)])
     display(['The relative error of LFC method is ', num2str(abs(numflux-exactflux)/abs(exactflux))]) 
    
    
    
end








