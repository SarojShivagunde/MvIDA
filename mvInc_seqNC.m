function [n_ij,n_i,n,mu_ij,Sjr,Djr] = mvInc_seqNC(x_n,NC,c,v,n_ij,n_i,n,mu_ij,Sjr,Djr,d_j)
%% 
% Sequential Increment - New sample belongs to a New class 'NC'
% Here, existing data samples means the data samples before the increment.
% This function takes as input-
%  1. x_n   - The set of new data samples
%  2. NC    - New class label (number)
%  3. c     - Number of existing classes
%  4. v     - Number of views
%  5. n_ij  - A (c,v) matrix of number of data samples per class per view
%  6. n_i   - A (c,1) vector of number of data samples per class
%  7. n     - Total number of data samples
%  8. mu_ij - A (max(d_j),c,v) matrix of per class per view means of existing samples
%  9. Sjr   - A (sum(d_j),sum(d_j)) within-class scatter matrix of existing samples
% 10. Djr   - A (sum(d_j),sum(d_j)) between-class scatter matrix of existing samples
% 11. d_j   - A (v,1) vector of the number of dimensions of each view

% Output arguments of this function are-
% 1. n_ij  - Updated n_ij
% 2. n_i   - Updated n_i
% 3. n     - Updated n
% 4. mu_ij - Updated mu_ij
% 5. Sjr   - Updated Sjr
% 6. Djr   - Updated Djr


%% Update 'n_ij' of 'NC' - Total no of samples of class 'NC' per view
for j = 1 : v
    n_ij(NC,j) = 1;
end

%% Update 'n_i' of 'NC' - Total no of samples of class 'NC'
n_i(NC) = v;

%% Update 'n' - Total no of samples
n = n + v;

%% Update 'mu_ij' where 'i = NC' - Mean of 'NC' per view
for j = 1 : v
    mu_ij(1:d_j(j),NC,j) = x_n(1:d_j(j),j);
end

%% Update 'Sjr' - Within-class Scatter
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
        if(j == r)
            Sjr(temp_dj:(temp_dj + d_j(j) - 1),temp_dr:(temp_dr + d_j(r)) - 1) = Sjr(temp_dj:(temp_dj + d_j(j) - 1),temp_dr:(temp_dr + d_j(r) - 1)) + (x_n(1:d_j(j),j) * x_n(1:d_j(j),j).' - (1 / v) * ((x_n(1:d_j(j),j) * x_n(1:d_j(j),j).')));
        else
            Sjr(temp_dj:(temp_dj + d_j(j) - 1),temp_dr:(temp_dr + d_j(r)) - 1) = Sjr(temp_dj:(temp_dj + d_j(j) - 1),temp_dr:(temp_dr + d_j(r) - 1)) - (1 / v) * (x_n(1:d_j(j),j) * x_n(1:d_j(r),r).');
        end
    end
end

%% Update 'Djr' - Between-class Scatter
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
            Temp_Djr1 = Temp_Djr1 + ((n_ij(i,j) * n_ij(i,r)) / n_i(i)) * (x_n(1:d_j(j),j) * x_n(1:d_j(r),r).');
        end
        for i = 1 : c
            Temp_Djr2 = Temp_Djr2 + (n_ij(i,j) * x_n(1:d_j(j),j));
        end
        for i = 1 : c
            Temp_Djr3 = Temp_Djr3 + (n_ij(i,r) * x_n(1:d_j(r),r));
        end
        if(r == 1)
            temp_dr = 1;
        else
            temp_dr = temp_dr + d_j(r - 1);
        end
        Djr(temp_dj:(temp_dj + d_j(j) - 1),temp_dr:(temp_dr + d_j(r) - 1)) = Temp_Djr1 - (1 / n) * (Temp_Djr2 * Temp_Djr3.');
    end
end
