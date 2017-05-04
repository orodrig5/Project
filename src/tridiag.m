function x = tridiag(n,a,b,c,f)
% x = tridiag(a,b,c,f)
% Solve an nxn tridiagonal system with sub?diagonal a, diagonal b,
% superdiagonal c, and rhs f
x = zeros(size(f));
c(1)=c(1)/b(1);
f(1)=f(1)/b(1);
% Forward elimination
for i=2:n
p = 1.0/(b(i)-c(i-1)*a(i));
c(i) = c(i)*p;
f(i) = (f(i) - a(i)*f(i-1))*p ;
end
% Back substitution
x(n) = f(n);
for i=n-1:-1:1
x(i) = (f(i) - c(i)*x(i+1));
end
