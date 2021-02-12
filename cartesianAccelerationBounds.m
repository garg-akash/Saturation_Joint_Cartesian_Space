function [ddXMin,ddXMax] = cartesianAccelerationBounds(X,dX,XMin,XMax,VXMin,VXMax,AXMin,AXMax)
global Ts
dXMin = max([(XMin-X)/Ts,VXMin,-sqrt(2*AXMax.*abs(X-XMin))],[],2);
dXMax = min([(XMax-X)/Ts,VXMax,sqrt(2*AXMax.*abs(XMax-X))],[],2);

ddXMin = max([2*(XMin-X-dX.*Ts)./Ts.^2, (dXMin-dX)./Ts, AXMin],[],2);
ddXMax = min([2*(XMax-X-dX.*Ts)./Ts.^2, (dXMax-dX)./Ts, AXMax],[],2);

idx0 = find((ddXMax<ddXMin) & (ddXMax<AXMin));
if idx0>0
    for x=1:size(idx0,1)
        if (ddXMax(idx0(x))<ddXMin(idx0(x)))
            ddXMax(idx0(x))=ddXMin(idx0(x));
        end
    end
end
idx1 = find((ddXMin>ddXMax) & (ddXMin>AXMax));
if idx1>0
    for x=1:size(idx1,1)
        if (ddXMin(idx1(x))>ddXMax(idx1(x)))
            ddXMin(idx1(x))=ddXMax(idx1(x));
        end
    end
end
end