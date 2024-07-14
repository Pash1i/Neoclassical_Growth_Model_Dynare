var a k c y i r;

varexo eps;

parameters ALPHA BETA DELTA RHO SIGMAE;

ALPHA = 0.33; // params values
BETA = 0.99;
DELTA = 0.025;
RHO = 0.9;
SIGMAE = 0.01;

model; // model equations
exp(c) + exp(k) - (1 - DELTA) * exp(k(-1)) = exp(a) * exp(k(-1))^ALPHA;
a = RHO * a(-1) + SIGMAE * eps;
c^(-1) = BETA * c(+1)^(-1) * (ALPHA * exp(a) * exp(k(-1))^(ALPHA - 1) + 1 - DELTA);
exp(y) = exp(a) * exp(k(-1))^(ALPHA);
exp(c) + exp(i) = exp(y);
r = ALPHA * exp(a) * exp(k(-1))^(ALPHA - 1);
end;

steady_state_model; // steady state values
a = 0;
k = log((1 / ALPHA * (1 / BETA - 1 + DELTA))^(1 / (ALPHA - 1)));
y = log(exp(k)^ALPHA);
i = log(DELTA * exp(k));
c = log(exp(y) - exp(i));
r = ALPHA * exp(k)^(ALPHA - 1);
end;

shocks;
var eps; stderr 1;
end;

steady;
check;

//stoch_simul(order=1, irf=0);
// To compute simulations comment the line above and uncomment the line below
stoch_simul(order=1, irf=40, periods=200);
