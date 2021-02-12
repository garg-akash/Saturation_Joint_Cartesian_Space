function taskScale = computeAccFactor(low, upp, acc)
global tol;
if acc < low - tol && acc < upp && low < 0
    taskScale = (low / acc);
elseif acc > upp + tol && acc > low && upp > 0
    taskScale = (upp / acc);
elseif acc >= low - tol && acc <= upp + tol
    taskScale = 1.0;
else
    taskScale = 0.0;
end
end