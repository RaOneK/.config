clc; clear; close all;

m = 3; 
n = 3;

nu = (m+1)*n;
nv = m*(n+1);
nc = m*n;
nn = (m+1)*(n+1);

% DIVERGENCE
D = zeros(nc,nu+nv);
Du = zeros(nc,nu);
Dv = zeros(nc,nv);

colShift = 0;
rowShift = 0; 
for j = 1:n % row iter
    for i = 1:m % fix row, run column
        Du(i+rowShift,i+colShift) = -1;
        Du(i+rowShift,i+1+colShift) = 1;
    end
    colShift = colShift + m+1;
    rowShift = rowShift + m;
end
% right part
for i = 1:nc
    Dv(i,i)=-1;
    Dv(i,i+m) = 1;
end
% concat
D = [Du Dv];

% GRAD
G = -D';

% CURL
C = zeros(nu+nv,nn);
Cu = zeros(nu,nn);
Cv = zeros(nv,nn);
% upper part
for i = 1:nu
    Cu(i,i)=-1;
    Cu(i,i+m+1) = 1;
end
% lower part
colShift = 0;
rowShift = 0; 
for j = 1:n+1
    for i = 1:m % run over each row
        Cv(i+rowShift,i+colShift) = 1;
        Cv(i+rowShift,i+1+colShift) = -1;
    end
    colShift = colShift + m+1;
    rowShift = rowShift + m;
end
%concat
C = [Cu;Cv];
%compute Laplacian
L = G*D - C*C';
