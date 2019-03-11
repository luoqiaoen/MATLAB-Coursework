function [u]=fromcirc2rect(uc,n)
% convert a circular signal of length n to rectangular signal of size 2n
n2=2*n;
u=zeros(n2,n2);
for i=1:n
    x=n-i+1;
    for j=1:n
        y=n-j+1;
        r1=sqrt(x*x+y*y);
        r=fix(r1+.5)+1;
        if r <=512
            ucr=uc(r);
            u(i,j)=ucr; 
            u(n2-i++1,j)=ucr;
            u(n2-i+1,n2-j+1)=ucr;
            u(i,n2-j+1)=ucr;
        end
    end
end
