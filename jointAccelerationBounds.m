function [ddQMin,ddQMax] = jointAccelerationBounds(q,dq,QMin,QMax,VJntMin,VJntMax,AJntMin,AJntMax)
global Ts
dQMin = max([(QMin-q)/Ts,VJntMin,-sqrt(2*AJntMax.*abs(q-QMin))],[],2);
dQMax = min([(QMax-q)/Ts,VJntMax,sqrt(2*AJntMax.*abs(QMax-q))],[],2);

ddQMin = max([2*(QMin-q-dq*Ts)/Ts^2, (dQMin-dq)/Ts, AJntMin],[],2);
ddQMax = min([2*(QMax-q-dq*Ts)/Ts^2, (dQMax-dq)/Ts, AJntMax],[],2);

idx0 = find(ddQMax<ddQMin & ddQMax<AJntMin);
if idx0>0
    for x=1:size(idx0,1)
        if (ddQMax(idx0(x))<ddQMin(idx0(x)))
            ddQMax(idx0(x))=ddQMin(idx0(x));
        end
    end
end
idx1 = find((ddQMin>ddQMax) & (ddQMin>AJntMax));
if idx1>0
    for x=1:size(idx1,1)
        if (ddQMin(idx1(x))>ddQMax(idx1(x)))
            ddQMin(idx1(x))=ddQMax(idx1(x));
        end
    end
end
end