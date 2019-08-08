function C = InitialData(x,y,h,t0,vel,pfun)
[xx,yy]=meshgrid(x,y);
xx=xx';
yy=yy';


C.VA=VInitialFun(xx-h/2,xx+h/2,yy-h/2,yy+h/2,t0,vel,pfun);
C.BU=BxInitialFun(yy+h/2,xx-h/2,xx+h/2,t0,vel,pfun);
C.BD=BxInitialFun(yy-h/2,xx-h/2,xx+h/2,t0,vel,pfun);
C.BR=ByInitialFun(xx+h/2,yy-h/2,yy+h/2,t0,vel,pfun);
C.BL=ByInitialFun(xx-h/2,yy-h/2,yy+h/2,t0,vel,pfun);


end

