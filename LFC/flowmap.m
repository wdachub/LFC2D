function timeline = flowmap(p, t0, te, vel, dt)
% compute the final position of point which started from p,
% p can be a vector
%  each point go through p at some instant between [t0, te].
% vel the flow velocity function
% dt the time increment
% n is the number of segment between interval

%t0, te can be a vector or a number




% dt = (te-t0)/n;
Npts=size(p,1);
timeline=zeros(Npts,2);

%t=[t0,te];
lt0=size(t0,1);
lte=size(te,1);
if lt0==1&&lte~=1
    lcase=1;
elseif lt0~=1&&lte==1
    lcase=2;
elseif lt0==1&&lte==1
    lcase=3;
else
    lcase=4;
end





switch lcase
    case 1
        for k=1:Npts
            if t0==te(k)
                %if the time interval is zero, ie t0=te, the flow map is identity
                %map
                timeline(k,:)=p(k,:);
            else
                ttmp=t0:dt:te(k);
                if ttmp(end)~=te(k)%|isempty(ttmp)
                    ttmp=[ttmp,te(k)];
                end
                ptmp=p(k,:);
                for i=1:length(ttmp)-1
                    ptmp=RungeKuttaf(ptmp, ttmp(i), ttmp(i+1), vel);
                end
                timeline(k,:)=ptmp;
            end
        end
    case 2
        for k=1:Npts
            if t0(k)==te
                %if the time interval is zero, ie t0=te, the flow map is identity
                %map
                timeline(k,:)=p(k,:);
            else
                ttmp=t0(k):dt:te;
                if ttmp(end)~=te%|isempty(ttmp)
                    ttmp=[ttmp,te];
                end
                ptmp=p(k,:);
                for i=1:length(ttmp)-1
                    ptmp=RungeKuttaf(ptmp, ttmp(i), ttmp(i+1), vel);
                end
                timeline(k,:)=ptmp;
            end
        end
    case 3
       
            if t0==te
                %if the time interval is zero, ie t0=te, the flow map is identity
                %map
                for k=1:Npts
                timeline(k,:)=p(k,:);
                end
            else
                ttmp=t0:dt:te;
                if ttmp(end)~=te%|isempty(ttmp)
                    ttmp=[ttmp,te];
                end
                for k=1:Npts
                ptmp=p(k,:);
                for i=1:length(ttmp)-1
                    ptmp=RungeKuttaf(ptmp, ttmp(i), ttmp(i+1), vel);
                end
                timeline(k,:)=ptmp;
                end
            end
        

        
    case 4
        for k=1:Npts
            if t0(k)==te(k)
                %if the time interval is zero, ie t0=te, the flow map is identity
                %map
                timeline(k,:)=p(k,:);
            else
                ttmp=t0(k):dt:te(k);
                if ttmp(end)~=te(k)%|isempty(ttmp)
                    ttmp=[ttmp,te(k)];
                end
                ptmp=p(k,:);
                for i=1:length(ttmp)-1
                    ptmp=RungeKuttaf(ptmp, ttmp(i), ttmp(i+1), vel);
                end
                timeline(k,:)=ptmp;
            end
        end
end
end





function pt = RungeKuttaf(p0, t0,te, vel)
% Runge-Kutta integrating velocity field.
%
% Input:
% p0     : the initial position of a point.  vertical or horizontal
% velocity field is important
% [t0,te]: the time interval;
% vel    : the velocity history,
% order  : order of the integrator
%
% Output:
% pt     : the ending position of the point
global RK
dt = te-t0;
% this case is never happened in this program.
% if dt==0
%     pt = p0;
%     return;
% end


pt = p0;
nS = max(size(RK.b));   % number of stages
sz = size(p0);
ks = zeros(sz(1), sz(2), nS);
for i=1:nS
    pp = p0;
    for j=1:i-1
        pp = pp + RK.a(i,j)*dt*ks(:,:,j);
    end
    ks(:,:,i) = vel(t0+RK.c(i)*dt, pp);
    pt = pt + dt*RK.b(i)*ks(:,:,i);
end
end
