import numpy


##
# Sequential Increment - New sample belongs to a New class 'NC'
# Here, existing data samples means the data samples before the increment.
# This function takes as input-
#  1. x_n   - The set of new data samples
#  2. NC    - New class label (number)
#  3. c     - Number of existing classes
#  4. v     - Number of views
#  5. n_ij  - A (c,v) matrix of number of data samples per class per view
#  6. n_i   - A (c,1) vector of number of data samples per class
#  7. n     - Total number of data samples
#  8. mu_ij - A (max(d_j),c,v) matrix of per class per view means of existing samples
#  9. Sjr   - A (sum(d_j),sum(d_j)) within-class scatter matrix of existing samples
# 10. Djr   - A (sum(d_j),sum(d_j)) between-class scatter matrix of existing samples
# 11. d_j   - A (v,1) vector of the number of dimensions of each view

# Output arguments of this function are-
# 1. n_ij  - Updated n_ij
# 2. n_i   - Updated n_i
# 3. n     - Updated n
# 4. mu_ij - Updated mu_ij
# 5. Sjr   - Updated Sjr
# 6. Djr   - Updated Djr
##

def SeqNC(x_n,NC,c,v,n_ij,n_i,n,mu_ij,Sjr,Djr,d_j):

    # Update 'n_ij' of 'NC' - Total no of samples of class 'NC' per view
    for j in range(1, v + 1):
        n_ij[NC,j] = 1

    # Update 'n_i' of 'NC' - Total no of samples of class 'NC'
    n_i[NC] = v

    # Update 'n' - Total no of samples
    n = n + v


    # Update 'mu_ij' where 'i = NC' - Mean of 'NC' per view
    for j in range(1, v + 1):
        mu_ij[1:d_j[j]+1,NC,j] = x_n[1:d_j[j],j]


    # Update 'Sjr' - Within-class Scatter
    for j in range(1, v + 1):
        if j == 1:
            temp_dj = 1
        else:
            temp_dj = temp_dj + d_j[j - 1]
        for r in range(1, v + 1):
            if r == 1:
                temp_dr = 1
            else:
                temp_dr = temp_dr + d_j[r - 1]
            if j == r:
                Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r]) - 1] = Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] + (x_n[1:d_j[j],j] * numpy.transpose(x_n[1:d_j[j],j]) - (1 / v) * ((x_n[1:d_j[j],j] * numpy.transpose(x_n[1:d_j[j],j]))))
            else:
                Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r]) - 1] = Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] - (1 / v) * (x_n[1:d_j[j],j] * numpy.transpose(x_n[1:d_j[r],r]))


    # Update 'Djr' - Between-class Scatter
    for j in range(1, v + 1):
        if j == 1:
            temp_dj = 1
        else:
            temp_dj = temp_dj + d_j[j - 1]
        for r in range(1, v + 1):
            Temp_Djr1 = numpy.zeros[d_j[j],d_j[r]]
            Temp_Djr2 = numpy.zeros[d_j[j],1]
            Temp_Djr3 = numpy.zeros[d_j[r],1]
            for i in range(1, c + 1):
                Temp_Djr1 = Temp_Djr1 + ((n_ij[i,j] * n_ij[i,r]) / n_i[i]) * (x_n[1:d_j[j],j] * numpy.transpose(x_n[1:d_j[r],r]))
            for i in range(1, c + 1):
                Temp_Djr2 = Temp_Djr2 + (n_ij[i,j] * x_n[1:d_j[j],j])
            for i in range(1, c + 1):
                Temp_Djr3 = Temp_Djr3 + (n_ij[i,r] * x_n[1:d_j[r],r])
            if r == 1:
                temp_dr = 1
            else:
                temp_dr = temp_dr + d_j[r - 1]
            Djr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] = Temp_Djr1 - (1 / n) * (Temp_Djr2 * numpy.transpose(Temp_Djr3))

    return [n_ij,n_i,n,mu_ij,Sjr,Djr];









## 
# Sequential Increment - New sample belongs to an existing class 'EC'
# Here, existing data samples means the data samples before the increment.
# This function takes as input-
#  1. x_n   - The set of new data samples
#  2. EC    - Class label (number) of new data sample
#  3. c     - Number of existing classes
#  4. v     - Number of views
#  5. n_ij  - A (c,v) matrix of number of data samples per class per view
#  6. n_i   - A (c,1) vector of number of data samples per class
#  7. n     - Total number of data samples
#  8. mu_ij - A (max(d_j),c,v) matrix of per class per view means of existing samples
#  9. Sjr   - A (sum(d_j),sum(d_j)) within-class scatter matrix of existing samples
# 10. Djr   - A (sum(d_j),sum(d_j)) between-class scatter matrix of existing samples
# 11. d_j   - A (v,1) vector of the number of dimensions of each view

# Output arguments of this function are-
# 1. n_ij  - Updated n_ij
# 2. n_i   - Updated n_i
# 3. n     - Updated n
# 4. mu_ij - Updated mu_ij
# 5. Sjr   - Updated Sjr
# 6. Djr   - Updated Djr
##

def SeqEC(x_n,EC,c,v,n_ij,n_i,n,mu_ij,Sjr,Djr,d_j):

    # Update 'Sjr' - Scatter within class 
    for j in range(1, v + 1):
        if j == 1:
            temp_dj = 1
        else:
            temp_dj = temp_dj + d_j[j - 1]
        for r in range(1, v + 1):
            if r == 1:
                temp_dr = 1
            else:
                temp_dr = temp_dr + d_j[r - 1]
            if j == r:
                Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] = Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] + (n_i[EC] / (v * (n_i[EC] + v))) * ((x_n[1:d_j[j],j] - mu_ij[1:d_j[j],EC,j]) * numpy.transpose(x_n[1:d_j[j],j] - mu_ij[1:d_j[j],EC,j])) + (x_n[1:d_j[j],j] * numpy.transpose(x_n[1:d_j[j],j]) - (1 / v) * ((x_n[1:d_j[j],j] * numpy.transpose(x_n[1:d_j[j],j]))))
            else:
                Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] = Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] + (n_i[EC] / (v * (n_i[EC] + v))) * ((x_n[1:d_j[j],j] - mu_ij[1:d_j[j],EC,j]) * numpy.transpose(x_n[1:d_j[r],r] - mu_ij[1:d_j[r],EC,r])) - ((1 / v) * ((x_n[1:d_j[j],j] * numpy.transpose(x_n[1:d_j[r],r]))))


    # Update 'mu_ij' where 'i = EC' - Mean of class 'EC' per view
    for j in range(1, v + 1):
        mu_ij[1:d_j[j],EC,j] = ((n_ij[EC,j] * mu_ij[1:d_j[j],EC,j]) + x_n[1:d_j[j],j]) / (n_ij[EC,j] + 1)
    end


    # Update 'n_ij' of 'EC' - Total no of samples of 'EC' per view
    for j in range(1, v + 1):
        n_ij[EC,j] = n_ij[EC,j] + 1
    end

    # Update 'n_i' of 'EC' - Total no of samples of class 'EC'
    n_i[EC] = n_i[EC] + v

    # Update 'n' - Total no of samples
    n = n + v


    # Update 'Djr' - Scatter between class
    for j in range(1, v + 1):
        if j == 1:
            temp_dj = 1
        else:
            temp_dj = temp_dj + d_j[j - 1]
        for r in range(1, v + 1):
            Temp_Djr1 = numpy.zeros(d_j[j],d_j[r])
            Temp_Djr2 = numpy.zeros(d_j[j],1)
            Temp_Djr3 = numpy.zeros(d_j[r],1)
            for i in range(1, c + 1):
                Temp_Djr1 = Temp_Djr1 + ((n_ij[i,j] * n_ij[i,r]) / n_i[i]) * (mu_ij[1:d_j[j],i,j] * numpy.transpose(mu_ij[1:d_j[r],i,r]))
            for i in range(1, c + 1):
                Temp_Djr2 = Temp_Djr2 + (n_ij[i,j] * mu_ij[1:d_j[j],i,j])
            for i in range(1, c + 1):
                Temp_Djr3 = Temp_Djr3 + (n_ij[i,r] * mu_ij[1:d_j[r],i,r])
            if r == 1:
                temp_dr = 1
            else:
                temp_dr = temp_dr + d_j[r - 1]
            Djr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] = Temp_Djr1 - (1 / n) * (Temp_Djr2 * numpy.transpose(Temp_Djr3))

    return [n_ij,n_i,n,mu_ij,Sjr,Djr];









##
# Chunk Increment - Some of the new samples belong to a new class 'NC'
# Here, existing data samples means the data samples before the increment.
# This function takes as input-
#  1. x_n   - The set of new data samples
#  2. c     - Number of classes after increment
#  3. c_old - Number of classes before increment
#  4. v     - Number of views
#  5. n_ij  - A (c_old,v) matrix of number of data samples per class per view
#  6. n_i   - A (c_old,1) vector of number of data samples per class
#  7. n     - Total number of data samples
#  8. L_ij  - A (c,v) matrix of number of data samples per class per view of new data samples
#  9. L_i   - A (c,1) vector of number of data samples per class of new data samples
# 10. L     - Total number of new data samples
# 11. mu_ij - A (max(d_j),c,v) matrix of per class per view means of existing samples
# 12. Sjr   - A (sum(d_j),sum(d_j)) within-class scatter matrix of existing samples
# 13. Djr   - A (sum(d_j),sum(d_j)) between-class scatter matrix of existing samples
# 14. d_j   - A (v,1) vector of the number of dimensions of each view

# Output arguments of this function are-
# 1. n_ij  - Updated n_ij
# 2. n_i   - Updated n_i
# 3. n     - Updated n
# 4. mu_ij - Updated mu_ij
# 5. Sjr   - Updated Sjr
# 6. Djr   - Updated Djr

def ChNC(x_n,c,c_old,v,n_ij,n_i,n,L_ij,L_i,L,mu_ij,Sjr,Djr,d_j):

    # Compute 'mu_x_ij' - Mean per class per view of new samples
    mu_x_ij = numpy.zeros(max(d_j),c,v)
    for i in range(1, c + 1):
        for j in range(1, v + 1):
            if L_ij[i,j] != 0:
                for k in range(1, L_ij[i,j]):
                    mu_x_ij[1:d_j[j],i,j] = mu_x_ij[1:d_j[j],i,j] + x_n[1:d_j[j],i,j,k]
                mu_x_ij[1:d_j[j],i,j] = mu_x_ij[1:d_j[j],i,j] / L_ij[i,j]
                

    # Update 'Sjr' - Within-class Scatter
    for j in range(1, v + 1):
        if j == 1:
            temp_dj = 1
        else:
            temp_dj = temp_dj + d_j[j - 1]
        for r in range(1, v + 1):
            if r == 1:
                temp_dr = 1
            else:
                temp_dr = temp_dr + d_j[r - 1]
            
            Temp_Sjr = numpy.zeros(d_j[j],d_j[r])
            for i in range(1, c + 1):
                if L_ij[i,j] != 0:
                    Temp_Sjr1 = numpy.zeros(d_j[j],d_j[r])
                    if  j == r:
                        for k in range(1, L_ij[i,j]):
                            Temp_Sjr1 = Temp_Sjr1 + x_n[1:d_j[j],i,j,k] * numpy.transpose(x_n[1:d_j[j],i,j,k])
                    Temp_Sjr2 = ((L_ij[i,j] * L_ij[i,r]) / L_i[i]) * (mu_x_ij[1:d_j[j],i,j] * numpy.transpose(mu_x_ij[1:d_j[r],i,r]))
                    if i <= c_old:
                        Temp_Sjr3 = (1/((n_i[i] + L_i[i]) * n_i[i] * L_i[i])) * ((L_ij[i,j] * n_i[i] * mu_x_ij[1:d_j[j],i,j] - (L_i[i] * n_ij[i,j] * mu_ij[1:d_j[j],i,j])) * numpy.transpose((L_ij[i,r] * n_i[i] * mu_x_ij[1:d_j[r],i,r]) - (L_i[i] * n_ij[i,r] * mu_ij[1:d_j[r],i,r])))
                        Temp_Sjr = Temp_Sjr + Temp_Sjr1 - Temp_Sjr2 + Temp_Sjr3
                    else:
                        Temp_Sjr = Temp_Sjr + Temp_Sjr1 - Temp_Sjr2
            
            Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] = Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] + Temp_Sjr

    # Update 'mu_ij'
    for i in range(1, c + 1):
        for j in range(1, v + 1):
            if L_ij[i,j] != 0:
                if i == c:
                    mu_ij[1:d_j[j],i,j] = mu_x_ij[1:d_j[j],i,j]
                else:
                    mu_ij[1:d_j[j],i,j] = ((n_ij[i,j] * mu_ij[1:d_j[j],i,j]) + (L_ij[i,j] * mu_x_ij[1:d_j[j],i,j])) / (n_ij[i,j] + L_ij[i,j])

    # Update 'n_ij'
    for i in range(1, c + 1):
        for j in range(1, v + 1):
            if L_ij[i,j] != 0:
                if i == c:
                    n_ij[i,j] = L_ij[i,j]
                else:
                    n_ij[i,j] = n_ij[i,j] + L_ij[i,j]

    # Update n_i
    for i in range(1, c + 1):
        if L_i[i] != 0:
            if i == c:
                n_i[i] = L_i[i]
            else:
                n_i[i] = n_i[i] + L_i[i]

    # Update 'n'
    n = n + L

    # Update 'Djr' - Between-class Scatter
    for j in range(1, v + 1):
        if j == 1:
            temp_dj = 1
        else:
            temp_dj = temp_dj + d_j[j - 1]
        for r in range(1, v + 1):
            Temp_Djr1 = numpy.zeros(d_j[j],d_j[r])
            Temp_Djr2 = numpy.zeros(d_j[j],1)
            Temp_Djr3 = numpy.zeros(d_j[r],1)
            for i in range(1, c + 1):
                Temp_Djr1 = Temp_Djr1 + ((n_ij[i,j] * n_ij[i,r]) / n_i[i]) * (mu_ij[1:d_j[j],i,j] * numpy.transpose(mu_ij[1:d_j[r],i,r]))
            for i in range(1, c + 1):
                Temp_Djr2 = Temp_Djr2 + (n_ij[i,j] * mu_ij[1:d_j[j],i,j])
            for i in range(1, c + 1):
                Temp_Djr3 = Temp_Djr3 + (n_ij[i,r] * mu_ij[1:d_j[r],i,r])
            if r == 1:
                temp_dr = 1
            else:
                temp_dr = temp_dr + d_j[r - 1]
            Djr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] = Temp_Djr1 - (1 / n) * (Temp_Djr2 * numpy.transpose(Temp_Djr3))

    return [n_ij,n_i,n,mu_ij,Sjr,Djr];









##
# Chunk Increment - New samples belong to some of the existing classes
# Here, existing data samples means the data samples before the increment.
# This function takes as input-
#  1. x_n   - The set of new data samples
#  2. c     - Number of classes after increment
#  3. v     - Number of views
#  4. n_ij  - A (c_old,v) matrix of number of data samples per class per view
#  5. n_i   - A (c_old,1) vector of number of data samples per class
#  6. n     - Total number of data samples
#  7. L_ij  - A (c,v) matrix of number of data samples per class per view of new data samples
#  8. L_i   - A (c,1) vector of number of data samples per class of new data samples
#  9. L     - Total number of new data samples
# 10. mu_ij - A (max(d_j),c,v) matrix of per class per view means of existing samples
# 11. Sjr   - A (sum(d_j),sum(d_j)) within-class scatter matrix of existing samples
# 12. Djr   - A (sum(d_j),sum(d_j)) between-class scatter matrix of existing samples
# 13. d_j   - A (v,1) vector of the number of dimensions of each view

# Output arguments of this function are-
# 1. n_ij  - Updated n_ij
# 2. n_i   - Updated n_i
# 3. n     - Updated n
# 4. mu_ij - Updated mu_ij
# 5. Sjr   - Updated Sjr
# 6. Djr   - Updated Djr

def ChEC(x_n,c,v,n_ij,n_i,n,L_ij,L_i,L,mu_ij,Sjr,Djr,d_j):

    # Compute 'mu_x_ij' - Mean per class per view of new samples
    mu_x_ij = zeros(max(d_j),c,v)
    for i in range(1, c + 1):
        for j in range(1, v + 1):
            if L_ij[i,j] != 0:
                for k in range(1, L_ij[i,j]):
                    mu_x_ij[1:d_j[j],i,j] = mu_x_ij[1:d_j[j],i,j] + x_n[1:d_j[j],i,j,k]
                mu_x_ij[1:d_j[j],i,j] = mu_x_ij[1:d_j[j],i,j] / L_ij[i,j]

    # Update 'Sjr' - Within-class Scatter
    for j in range(1, v + 1):
        if j == 1:
            temp_dj = 1
        else:
            temp_dj = temp_dj + d_j[j - 1]
        for r in range(1, v + 1):
            if r == 1:
                temp_dr = 1
            else:
                temp_dr = temp_dr + d_j[r - 1]
            
            Temp_Sjr = zeros(d_j[j],d_j[r])
            for i in range(1, c + 1):
                if L_ij[i,j] != 0:
                    Temp_Sjr1 = zeros(d_j[j],d_j[r])
                    if j == r:
                        for k in range (1, L_ij[i,j]):
                            Temp_Sjr1 = Temp_Sjr1 + x_n[1:d_j[j],i,j,k] * numpy.transpose(x_n[1:d_j[j],i,j,k])
                    Temp_Sjr2 = ((L_ij[i,j] * L_ij[i,r]) / L_i[i]) * (mu_x_ij[1:d_j[j],i,j] * numpy.transpose(mu_x_ij[1:d_j[r],i,r]))
                    Temp_Sjr3 = (1/((n_i[i] + L_i[i]) * n_i[i] * L_i[i])) * ((L_ij[i,j] * n_i[i] * mu_x_ij[1:d_j[j],i,j] - (L_i[i] * n_ij[i,j] * mu_ij[1:d_j[j],i,j])) * numpy.transpose((L_ij[i,r] * n_i[i] * mu_x_ij[1:d_j[r],i,r]) - (L_i[i] * n_ij[i,r] * mu_ij[1:d_j[r],i,r])))
                    Temp_Sjr = Temp_Sjr + Temp_Sjr1 - Temp_Sjr2 + Temp_Sjr3
            Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] = Sjr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] + Temp_Sjr

    # Update 'mu_ij'
    for i in range(1, c + 1):
        for j in range(1, v + 1):
            if L_ij[i,j] != 0:
                mu_ij[1:d_j[j],i,j] = ((n_ij[i,j] * mu_ij[1:d_j[j],i,j]) + (L_ij[i,j] * mu_x_ij[1:d_j[j],i,j])) / (n_ij[i,j] + L_ij[i,j])

    # Update 'n_ij'
    for i in range(1, c + 1):
        for j in range(1, v + 1):
            if L_ij[i,j] != 0:
                n_ij[i,j] = n_ij[i,j] + L_ij[i,j]

    # Update 'n_i'
    for i in range(1, c + 1):
        if L_i[i] != 0:
            n_i[i] = n_i[i] + L_i[i]

    # Update 'n'
    n = n + L

    # Update 'Djr' - Between-class Scatter
    for j in range(1, v + 1):
        if j == 1:
            temp_dj = 1
        else:
            temp_dj = temp_dj + d_j[j - 1]
        for r in range(1, v + 1):
            Temp_Djr1 = zeros(d_j[j],d_j[r])
            Temp_Djr2 = zeros(d_j[j],1)
            Temp_Djr3 = zeros(d_j[r],1)
            for i in range(1, c + 1):
                Temp_Djr1 = Temp_Djr1 + ((n_ij[i,j] * n_ij[i,r]) / n_i[i]) * (mu_ij[1:d_j[j],i,j] * numpy.transpose(mu_ij[1:d_j[r],i,r]))
            for i in range(1, c + 1):
                Temp_Djr2 = Temp_Djr2 + (n_ij[i,j] * mu_ij[1:d_j[j],i,j])
            for i in range(1, c + 1):
                Temp_Djr3 = Temp_Djr3 + (n_ij[i,r] * mu_ij[1:d_j[r],i,r])
            if r == 1:
                temp_dr = 1
            else:
                temp_dr = temp_dr + d_j[r - 1]
            Djr[temp_dj:(temp_dj + d_j[j] - 1),temp_dr:(temp_dr + d_j[r] - 1)] = Temp_Djr1 - (1 / n) * (Temp_Djr2 * numpy.transpose(Temp_Djr3))

    return [n_ij,n_i,n,mu_ij,Sjr,Djr];
                                                                                                    
