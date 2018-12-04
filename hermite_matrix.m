function H = hermite_matrix(x,y,z)
for i=0:2
    for j=0:2
        for k=0:2
            H(i+1,j+1,k+1) = hermite(i,x)*hermite(j,y)*hermite(k,z);
        end
    end
end
