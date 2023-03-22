function rel = eci2cwh(chief, dep)
    
    r = chief(1:3);
    v = chief(4:6);
    
    h = cross(r, v);
    
    exhat = r ./ norm(r);
    ezhat = h ./ norm(h);
    eyhat = cross(exhat, ezhat);
    
    exhatdot = (v - dot(exhat, v) * exhat) ./ norm(r);
    ezhatdot = zeros(3, 1);
    eyhatdot = cross(exhatdot, ezhat);
    
    C = [exhat, eyhat, ezhat]';
    Cdot = [exhatdot, eyhatdot, ezhatdot]';
    
    rot = [C, zeros(3); Cdot, C];
    rel = rot*(dep - chief);
end