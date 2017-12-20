% This function solves (LU)*v = d where L and U are LU factors of a tridiagonal
% matrix

function v = solve_Aud(d,a,e,f)
% Input:     d: right hand side vector of the system of equations
%          e,f: vectors defining the L and U factors of the tridiagonal
%          matrix. e and f are obtained with the LU_triag function
%
% Output: v = solution vector. 


% Initialize v
n = length(d); v = zeros(n,1); 

% Forward substitution to solve L*s = d. This is possible because L is
% lower triangular
v(1) = d(1)/e(1);
for i=2:n
  v(i) = (d(i) - a(i)*v(i-1))/e(i);
end

% Backward substitution to solve U*v = s. This is possible becauseU is
% upper triangular. The elements of v are overwritten.
for i=n-1:-1:1
  v(i) = v(i) - f(i)*v(i+1);
end



