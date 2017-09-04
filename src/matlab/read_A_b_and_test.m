% settings
EPS = 1E-12;
maxItnMult = 100;
folder = '../../data/random_mat/benchmark_compatible_5000_2500';

read_mat  = 1;
run_lsqr  = 1;
run_lsmr  = 1;
run_bcgs  = 0;
run_bcgsp = 0;

% Start
if read_mat
    fprintf('Reading\n')
    % % LQ and QR data should be the same
    A = dlmread(sprintf('%s/A.txt', folder));
    b = dlmread(sprintf('%s/b.txt', folder));

    A_cond = cond(A); 
    A_rank = rank(A);
    [U,S,V] = svd(A);
    Nsv = max(size(S));
    A_svs = zeros(Nsv,1); 
    for ii = 1:Nsv
        if ii > min(size(S))
            A_sys(ii) = 0.0;
        else
            A_svs(ii) = S(ii,ii);
        end
    end
    dlmwrite(sprintf('%s/A_cond.txt'           , folder), A_cond); 
    dlmwrite(sprintf('%s/A_rank.txt'           , folder), A_rank); 
    dlmwrite(sprintf('%s/A_singular_values.txt', folder), A_svs ); 
end

maxItn = max(max(size(A,1), size(A,2)),size(b,1))*maxItnMult;

if run_lsqr == 1
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
end
 
if run_lsmr == 1
    fprintf('Performing LSMR\n')
    tic
    [lsmr_x, lsmr_flag, lsmr_ite, lsmr_resl, lsmr_resv] = lsmr(A, b, 0, EPS, EPS, [], maxItn, [], 1);
    lsmr_t = toc
    lsmr_r = A*lsmr_x - b; 
    dlmwrite(sprintf('%s/lsmr_x.txt'   , folder), lsmr_x    ); 
    dlmwrite(sprintf('%s/lsmr_flag.txt', folder), lsmr_flag ); 
    dlmwrite(sprintf('%s/lsmr_resl.txt' , folder), lsmr_resl); 
    dlmwrite(sprintf('%s/lsmr_ite.txt' , folder), lsmr_ite  ); 
    dlmwrite(sprintf('%s/lsmr_t.txt'   , folder), lsmr_t    ); 
    dlmwrite(sprintf('%s/lsmr_resv.txt', folder), lsmr_resv ); 
    dlmwrite(sprintf('%s/lsmr_res.txt' , folder), lsmr_r    ); 
end

if run_bcgs == 1
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
end

if run_bcgsp == 1
    fprintf('Performing BiCGStab with ILU preconditioner\n')
    tic
    [L,U] = ilu(sparse(A), struct('type','ilutp','droptol',1e-6));
    [bcgsp_x, bcgsp_flag, bcgsp_resl, bcgsp_ite, bcgsp_resv] = bicgstab(A, b, EPS, maxItn, L, U);
    bcgsp_t = toc
    bcgsp_r = A*bcgsp_x - b; 
    dlmwrite(sprintf('%s/bcgsp_x.txt'   , folder), bcgsp_x   ); 
    dlmwrite(sprintf('%s/bcgsp_flag.txt', folder), bcgsp_flag); 
    dlmwrite(sprintf('%s/bcgsp_resl.txt', folder), bcgsp_resl); 
    dlmwrite(sprintf('%s/bcgsp_ite.txt' , folder), bcgsp_ite ); 
    dlmwrite(sprintf('%s/bcgsp_t.txt'   , folder), bcgsp_t   ); 
    dlmwrite(sprintf('%s/bcgsp_resv.txt', folder), bcgsp_resv); 
    dlmwrite(sprintf('%s/bcgsp_res.txt' , folder), bcgsp_r   ); 
end
