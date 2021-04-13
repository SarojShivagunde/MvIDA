function [n_ij,n_i,n,mu_ij,Sjr,Djr] = mvInc_chEC(x_n,c,v,n_ij,n_i,n,L_ij,L_i,L,mu_ij,Sjr,Djr,d_j)
%%
% Chunk Increment - New samples belong to some of the existing classes
% Here, existing data samples means the data samples before the increment.
% This function takes as input-
%  1. x_n   - The set of new data samples
%  2. c     - Number of classes after increment
%  3. v     - Number of views
%  4. n_ij  - A (c_old,v) matrix of number of data samples per class per view
%  5. n_i   - A (c_old,1) vector of number of data samples per class
%  6. n     - Total number of data samples
%  7. L_ij  - A (c,v) matrix of number of data samples per class per view of new data samples
%  8. L_i   - A (c,1) vector of number of data samples per class of new data samples
%  9. L     - Total number of new data samples
% 10. mu_ij - A (max(d_j),c,v) matrix of per class per view means of existing samples
% 11. Sjr   - A (sum(d_j),sum(d_j)) within-class scatter matrix of existing samples
% 12. Djr   - A (sum(d_j),sum(d_j)) between-class scatter matrix of existing samples
% 13. d_j   - A (v,1) vector of the number of dimensions of each view

% Output arguments of this function are-
% 1. n_ij  - Updated n_ij
% 2. n_i   - Updated n_i
% 3. n     - Updated n
% 4. mu_ij - Updated mu_ij
% 5. Sjr   - Updated Sjr
% 6. Djr   - Updated Djr

%% Compute mu_x_ij - Mean per class per view of new samples
mu_x_ij = zeros(max(d_j),c,v);
for i = 1 : c
    for j = 1 : v
        if (L_ij(i,j) ~= 0)
            for k = 1 : L_ij(i,j)
                mu_x_ij(1:d_j(j),i,j) = mu_x_ij(1:d_j(j),i,j) + x_n(1:d_j(j),i,j,k);
            end
            mu_x_ij(1:d_j(j),i,j) = mu_x_ij(1:d_j(j),i,j) / L_ij(i,j);
        end
    end
end

%% Update Sjr - Within-class Scatter
for j = 1 : v
    if(j == 1)
        temp_dj = 1;
    else
        temp_dj = temp_dj + d_j(j - 1);
    end
    for r = 1 : v
        if(r == 1)
            temp_dr = 1;
        else
            temp_dr = temp_dr + d_j(r - 1);
        end
        
        Temp_Sjr = zeros(d_j(j),d_j(r));
        for i = 1 : c
            if (L_ij(i,j) ~= 0)
                Temp_Sjr1 = zeros(d_j(j),d_j(r));
                if (j == r)
                    for k = 1 : L_ij(i,j)
                        Temp_Sjr1 = Temp_Sjr1 + x_n(1:d_j(j),i,j,k) * x_n(1:d_j(j),i,j,k).';
                    end
                end
                Temp_Sjr2 = ((L_ij(i,j) * L_ij(i,r)) / L_i(i)) * (mu_x_ij(1:d_j(j),i,j) * mu_x_ij(1:d_j(r),i,r).');
                Temp_Sjr3 = (1/((n_i(i) + L_i(i)) * n_i(i) * L_i(i))) * ((L_ij(i,j) * n_i(i) * mu_x_ij(1:d_j(j),i,j) - (L_i(i) * n_ij(i,j) * mu_ij(1:d_j(j),i,j))) * ((L_ij(i,r) * n_i(i) * mu_x_ij(1:d_j(r),i,r)) - (L_i(i) * n_ij(i,r) * mu_ij(1:d_j(r),i,r))).');
                Temp_Sjr = Temp_Sjr + Temp_Sjr1 - Temp_Sjr2 + Temp_Sjr3;
            end
        end
        Sjr(temp_dj:(temp_dj + d_j(j) - 1),temp_dr:(temp_dr + d_j(r) - 1)) = Sjr(temp_dj:(temp_dj + d_j(j) - 1),temp_dr:(temp_dr + d_j(r) - 1)) + Temp_Sjr;
    end
end

%% Update mu_ij
for i = 1 : c
    for j = 1 : v
        if (L_ij(i,j) ~= 0)
            mu_ij(1:d_j(j),i,j) = ((n_ij(i,j) * mu_ij(1:d_j(j),i,j)) + (L_ij(i,j) * mu_x_ij(1:d_j(j),i,j))) / (n_ij(i,j) + L_ij(i,j));
        end
    end
end
clear mu_x_ij;

%% Update n_ij
for i = 1 : c
    for j = 1 : v
        if (L_ij(i,j) ~= 0)
            n_ij(i,j) = n_ij(i,j) + L_ij(i,j);
        end
    end
end

%% Update n_i
for i = 1 : c
    if (L_i(i) ~= 0)
        n_i(i) = n_i(i) + L_i(i);
    end
end

%% Update n
n = n + L;

%% Update Djr - Between-class Scatter
for j = 1 : v
    if(j == 1)
        temp_dj = 1;
    else
        temp_dj = temp_dj + d_j(j - 1);
    end
    for r = 1 : v
        Temp_Djr1 = zeros(d_j(j),d_j(r));
        Temp_Djr2 = zeros(d_j(j),1);
        Temp_Djr3 = zeros(d_j(r),1);
        for i = 1 : c
            Temp_Djr1 = Temp_Djr1 + ((n_ij(i,j) * n_ij(i,r)) / n_i(i)) * (mu_ij(1:d_j(j),i,j) * mu_ij(1:d_j(r),i,r).');
        end
        for i = 1 : c
            Temp_Djr2 = Temp_Djr2 + (n_ij(i,j) * mu_ij(1:d_j(j),i,j));
        end
        for i = 1 : c
            Temp_Djr3 = Temp_Djr3 + (n_ij(i,r) * mu_ij(1:d_j(r),i,r));
        end
        if(r == 1)
            temp_dr = 1;
        else
            temp_dr = temp_dr + d_j(r - 1);
        end
        Djr(temp_dj:(temp_dj + d_j(j) - 1),temp_dr:(temp_dr + d_j(r) - 1)) = Temp_Djr1 - (1 / n) * (Temp_Djr2 * Temp_Djr3.');
    end
end


