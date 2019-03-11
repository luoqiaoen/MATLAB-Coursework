function Y = projection(X,imagetype)
%% apply image constraint, such as real positivity (imagetype==1) and
%% complex positivity (imagetype==2)
[N1 , N2] = size(X);

if imagetype == 1
    Y  = real(X);
    for m = 1 : N1
        for l = 1 : N2
            if Y(m,l) < 0
                Y(m,l) = 0;
            end
        end
    end
elseif imagetype == 2
    newreal  = real(X);
    newimag  = imag(X);
    for m = 1 : N1
        for l = 1 : N2
            if newreal(m,l) < 0
                newreal(m,l) = 0;
            end
            if newimag(m,l) < 0
                newimag(m,l) = 0;
            end
        end
    end
    Y  = newreal + sqrt(-1) * newimag;
elseif imagetype == -2
    Y = X;
end