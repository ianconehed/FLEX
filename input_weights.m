function y = input_weights(num_columns,npp,N,num_VTA,gain,VTA_gain)
%instantiation of weight matrix
rho_1 = .7;
rho_2 = .9;
L_ij1 = zeros(N,N);
ndd = npp;
for i = 1:num_columns
    j = (i-1)*npp*2;
    l = (i-1)*ndd;
    L_ij1(j+1:j+npp,l+1:l+ndd) =  gain*eye(npp);
    a = ones(num_VTA/npp,1);
    for ii = 1: npp -1
        a = blkdiag(a,ones(num_VTA/npp,1));
    end
    a = a(randperm(size(a, 1)), :);
    b = a(randperm(size(a, 1)), :);
    L_ij1(N-num_VTA+1:N,l+1:l+ndd) = .05*VTA_gain*spones(sprandn(num_VTA,npp,rho_1)).*a;
    L_ij1(N-num_VTA+1:N,N-num_VTA+1:N) = VTA_gain*(sprand(num_VTA,num_VTA,rho_2)).*eye(num_VTA);
    L_ij1(N-num_VTA+1:N,N-num_VTA+1:N) = L_ij1(N-num_VTA+1:N,N-num_VTA+1:N) + VTA_gain*(L_ij1(N-num_VTA+1:N,N-num_VTA+1:N)>0);
end
y = abs(L_ij1);
end
