function y = Sparse_L_ij(num_columns,npp,N,num_VTA,l5_rec,l23_rec,l5_l23,l23_l5,l23_l23_cross,l5_l5_cross,l_23_l_VTA,l_VTA_l_VTA,noise)
%instantiation of weight matrix
L_ij1 = zeros(N,N);
rho = .3;%.0675;
rho_2 = .1;
for i = 1:num_columns
    j = (i-1)*npp*2;
    if i < num_columns
        L_ij1(j+npp+1:j+(2*npp),j+(2*npp)+1:j+npp*3) = l23_l5*spones(sprandn(npp,npp,rho));
    end
    for k = 0:2*npp:(num_columns - 1)*2*npp
        if k ~= j
            L_ij1(j+npp+1:j+2*npp,k+npp+1:k+2*npp) = l23_l23_cross*spones(sprandn(npp,npp,rho));
            L_ij1(j+1:j+npp,k+1:k+npp) = l5_l5_cross*spones(sprandn(npp,npp,rho));
        end
    end
    L_ij1(j+1:j+npp,j+1:j+npp) =  l5_rec*spones(sprandn(npp,npp,rho));
    L_ij1(j+npp+1:j+2*npp,j+npp+1:j+2*npp) =  l23_rec*spones(sprandn(npp,npp,rho));
    L_ij1(j+1:j+npp,j+npp+1:j+2*npp) =  l5_l23*spones(sprandn(npp,npp,rho));
    L_ij1(j+1:j+npp,j+npp+1:j+2*npp) =  l5_l23*spones(sprandn(npp,npp,rho));
    L_ij1(N-num_VTA+1:N,N-num_VTA+1:N) =  l_VTA_l_VTA*spones(sprandn(num_VTA,num_VTA,rho));

    a = ones(num_VTA/npp,1);
    for ii = 1: npp -1
        a = blkdiag(a,ones(num_VTA/npp,1));
    end
    a = a(randperm(size(a, 1)), :);
    L_ij1(j+npp+1:j+(2*npp),N-num_VTA+1:N) = l_23_l_VTA*(sprandn(num_VTA,num_VTA,rho_2));
    L_ij1(j+npp+1:j+(2*npp),N-num_VTA+1:N) = L_ij1(j+npp+1:j+(2*npp),N-num_VTA+1:N) + l_23_l_VTA*(L_ij1(j+npp+1:j+(2*npp),N-num_VTA+1:N)>0);
end
L_ij1 = L_ij1 + noise*randn(N).*(L_ij1>0);
% y = abs(transpose(L_ij1(1:N,1:N)));
y = sparse(abs(transpose(L_ij1(1:N,1:N))));
end
