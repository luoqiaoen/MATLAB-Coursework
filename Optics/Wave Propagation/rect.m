function zz = rect(x, y, D)
if nargin == 2, 
    D = 1;
end
x = abs(x);
xx = double(x<D/2);
xx(x == D/2) = 0.5;
y = abs(y);
yy = double(y<D/2);
yy(y == D/2) = 0.5;
zz = xx.*yy;