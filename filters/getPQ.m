function [P,Q,inRange] = getPQ(P,Q,multiplier)
inRange = false;
maxPQ = 2^31; %this is a number imposed by matlab
P = round(P*multiplier);
Q = round(Q*multiplier);
gcd_ = gcd(P,Q);
P = P/gcd_;
Q = Q/gcd_;
PQ = P*Q;

%%
if PQ < maxPQ
    inRange = true;
end