function [y, J] = gaussfun(beta,x)

if length(beta) ~= 3
    error('beta parameter in gaussfit of wrong dimensions');
end


nl = length(x);
y = zeros(1,nl);
J = zeros(nl,3);
for i=1:nl
    dx = x(i) - beta(2);
    if ( beta(3) <= 0)
        expval = 0.0;
    else
        expval = exp(-0.5*dx*dx/(beta(3)*beta(3)));
    end
    y(1,i) = beta(1)*expval;
    J(i,1) = expval;
    if ( beta(3) <= 0 )
        J(i,2) = 0.0;
        J(i,3) = 0.0;
    else
        J(i,2) = beta(1)*dx*expval/(beta(3)*beta(3));
        J(i,3) = beta(1)*dx*dx*expval/(beta(3)*beta(3)*beta(3));
    end
end

