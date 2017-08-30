fprintf('Reading\n')
A = dlmread('../../A.txt');
b = dlmread('../../b.txt');
maxItn = max(max(size(A,1), size(A,2)),size(b,1))*100;

fprintf('Performing LSQR\n')
tic
lsqr(A, b, 1E-12, maxItn );
toc
fprintf('Performing LSMR\n')
tic
lsmr(A, b, 0, 1E-12, 1E-12);
toc
fprintf('Performing LSMR\n')
tic
bicgstab(A, b, 1E-12, maxItn);
toc
