function sse = f_fit(params, V, f)

approx = params(1)*(params(2)-V+params(3)*exp((V-params(4))/params(3)));

Error_Vector = approx - f;
% When curvefitting, a typical quantity to
% minimize is the sum of squares error
sse = mean(Error_Vector.^2);