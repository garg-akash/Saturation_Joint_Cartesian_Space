function [sJnt,sCar] = getAccelerationRatios(WJnt,WCar,ddq,ddQMin,ddQMax,ddX,ddXMin,ddXMax)
nJnt = size(ddq,1);
mCar = size(ddX,1);
sJnt = ones(nJnt,1);
sCar = ones(mCar,1);
for i=1:nJnt
    if WJnt(i,i) == 0
        sJnt(i) = inf;
    else
        sJnt(i) = computeAccFactor(ddQMin(i), ddQMax(i), ddq(i));
    end
end
for i=1:mCar
    if WCar(i,i) == 0
        sCar(i) = inf;
    else
        sCar(i) = computeAccFactor(ddXMin(i), ddXMax(i), ddX(i));
    end
end
end