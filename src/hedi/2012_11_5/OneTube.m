function ot=OneTube(p,tau1,alpha1)

% a(t)= p(t)+alpha1*b(t-tau1)
% b(t) = a(t-tau1)
% ot=bt;

n=length(p);
b=zeros(n,1);
for i=1:n;
    a(i)=p(i);
    if i>tau1;
        a(i)=a(i)-alpha1*b(i-tau1);
        b(i)=alpha1*a(i-tau1);
    end;
end;

ot=b(:);

    