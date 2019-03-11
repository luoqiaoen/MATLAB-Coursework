clear all ; clc
load('data.mat');

M = 2;
K = 3;
pik = [1/3, 1/3, 1/3];

R = zeros(2,2,3);
R(:,:,1) = [1 0 ; 0 1];
R(:,:,2) = [1 0 ; 0 1];
R(:,:,3) = [1 0 ; 0 1];

u = [2,2 ; -2,-2 ;5.5, 2];

Size = size(x);
length = Size(1);

pn = zeros(1,3);
iterations = 20;

for iter = 1: iterations
    
Nk = zeros(1,3);

t1k = zeros(2,3);

t2k = zeros(2,2,3);
    for n = 1:1:length
        for i = 1:1:K
            pn(i) =pik(i)*(det(R(:,:,i))^(-0.5)).*exp(-0.5*(x(n,:)-u(i,:))*inv(R(:,:,i))*(x(n,:)-u(i,:))')./(2*pi)^(M/2); %#ok<MINV>
        end
        pd = sum(pn);
        
        Nk = Nk + pn/pd;
        
        t1k = t1k + ((pn/pd)'*x(n,:))';
        
        for i = 1:1:K
            t2k(:,:,i) = t2k(:,:,i)  + (pn(i)/pd)'*(x(n,:)'*x(n,:));
        end
    end
    fprintf('Iteration = %d \n',iter)
    
    for i = 1:1:K
        u(i,:) = (t1k(:,i)./Nk(i))'
        
        R(:,:,i) = t2k(:,:,i)./Nk(i)-t1k(:,i)*t1k(:,i)'/(Nk(i).^2)
        
        pik(i) = Nk(i)./length
    end
end
