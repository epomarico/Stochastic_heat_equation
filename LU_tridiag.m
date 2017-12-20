% This function calculates the vectors defining the matrices L and U in which 
% a tridiagonal matrix can be factorized 

function [e,f] = LU_tridiag(a,b,c)
 % Input:  a,b,c: vectors defining the tridiagonal matrix to factorize  
 %                a is the subdiagonal
 %                b is the main diagonal
 %                c is the superdiagonal
 %
 % Output:  e,f: vectors defining the L and U factors of the tridiagonal matrix
 
 n = length(a);
 e = zeros(n,1);  f = zeros(n,1);
 e(1) = b(1);
 f(1) = c(1)/b(1);
 
 for i=2:n
     e(i) = b(i) - a(i)*f(i-1);
     f(i) = c(i)/e(i);
 end