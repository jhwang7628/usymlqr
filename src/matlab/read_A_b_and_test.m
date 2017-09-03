% settings
folder = '../../data/random_mat/benchmark_5000';

% fprintf('Reading\n')
% % LQ and QR data should be the same
A = dlmread(sprintf('%s/A.txt', folder));
b = dlmread(sprintf('%s/b.txt', folder));

maxItn = max(max(size(A,1), size(A,2)),size(b,1))*100;
EPS = 1E-12;

fprintf('Performing LSQR\n')
tic
[lsqr_x, lsqr_flag, lsqr_resl, lsqr_ite, lsqr_resv] = lsqr(A, b, EPS, maxItn);
lsqr_t = toc
lsqr_r = A*lsqr_x - b; 
dlmwrite(sprintf('%s/lsqr_x.txt'   , folder), lsqr_x    ); 
dlmwrite(sprintf('%s/lsqr_flag.txt', folder), lsqr_flag ); 
dlmwrite(sprintf('%s/lsqr_resl.txt', folder), lsqr_resl ); 
dlmwrite(sprintf('%s/lsqr_ite.txt' , folder), lsqr_ite  ); 
dlmwrite(sprintf('%s/lsqr_t.txt'   , folder), lsqr_t    ); 
dlmwrite(sprintf('%s/lsqr_resv.txt', folder), lsqr_resv ); 
dlmwrite(sprintf('%s/lsqr_res.txt' , folder), lsqr_r    ); 
  
fprintf('Performing LSMR\n')
tic
[lsmr_x, lsmr_flag, lsmr_ite, lsmr_res] = lsmr(A, b, 0, EPS, EPS, [], maxItn, [], 1)
lsmr_t = toc
lsmr_r = A*lsmr_x - b; 
dlmwrite(sprintf('%s/lsmr_x.txt'   , folder), lsmr_x   ); 
dlmwrite(sprintf('%s/lsmr_flag.txt', folder), lsmr_flag); 
dlmwrite(sprintf('%s/lsmr_res.txt' , folder), lsmr_res ); 
dlmwrite(sprintf('%s/lsmr_ite.txt' , folder), lsmr_ite ); 
dlmwrite(sprintf('%s/lsmr_t.txt'   , folder), lsmr_t   ); 
dlmwrite(sprintf('%s/lsmr_res.txt' , folder), lsmr_r   ); 

fprintf('Performing BiCGStab\n')
tic
[bcgs_x, bcgs_flag, bcgs_resl, bcgs_ite, bcgs_resv] = bicgstab(A, b, EPS, maxItn);
bcgs_t = toc
bcgs_r = A*bcgs_x - b; 
dlmwrite(sprintf('%s/bcgs_x.txt'   , folder), bcgs_x   ); 
dlmwrite(sprintf('%s/bcgs_flag.txt', folder), bcgs_flag); 
dlmwrite(sprintf('%s/bcgs_resl.txt', folder), bcgs_resl); 
dlmwrite(sprintf('%s/bcgs_ite.txt' , folder), bcgs_ite ); 
dlmwrite(sprintf('%s/bcgs_t.txt'   , folder), bcgs_t   ); 
dlmwrite(sprintf('%s/bcgs_resv.txt', folder), bcgs_resv); 
dlmwrite(sprintf('%s/bcgs_res.txt' , folder), bcgs_r   ); 
% preconditioned BiCGStab
% [L,U] = ilu(sparse(A), struct('type','ilutp','droptol',1e-6));
% [bcgsp_x, bcgsp_flag, bcgsp_resl, bcgsp_ite, bcgsp_resv] = bicgstab(A, b, EPS, maxItn);
