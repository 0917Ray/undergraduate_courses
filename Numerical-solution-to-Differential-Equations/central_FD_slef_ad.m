function [x, U] = central_FD_slef_ad(x0, xn, ua, a, b, c, p, q, f, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This matlab function solves the following two-point                   %
%   self-adjoint Boundary Value Problem:                                  %
%   (p(x)u')' - q(x)u(x) = f(x), x0 < x < xn,                             %
%             u(x0) = ua, au(xn) bu'(xn) = c.                             %
%   Input:                                                                %
%   x0, xn: Two end points.                                               %
%   ua: Dirichlet boundary condition at x0.                               %
%   a, b, c: Robin boundary condition at xn.                              %
%   p: diffusion term.                                                    %
%   q: react term.                                                        %
%   f: force function f(x).                                               %
%   n: number of grid points.                                             %
%   Output:                                                               %
%   x: x(1),x(2),...,x(n) are grid points.                                %
%   U: U(1),U(2),...,U(n) are approximate solution at grid points.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = (xn-x0)/n; % step length
h1 = h*h;

A = sparse(n,n);
F = zeros(n,1);

x(1) = x0 + h;
for i=2:n
    x(i) = x(i-1)+h;
end

for i=1:n-1
    A(i,i) = -(p(x(i)-h/2)+p(x(i)+h/2))/h1-q(x(i));
    A(i,i+1) = p(x(i)+h/2)/h1; A(i+1,i) = A(i,i+1);
end

% A(n-1,n-1) = a+b/h;
% A(n-1,n) = p(x(i-1)+h/2)/h1;
A(n,n-1) = -b/h; A(n,n) = a+b/h;


for i=2:n-1
    F(i) = f(x(i));
end

F(1) = f(x(1)) - p(x0+h/2)*ua/h1;
F(n) = c;

U = A\F;

return
end




