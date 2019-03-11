clear all ; clc
load('data.mat');

M = 2;
K = 9;

for MDL = K:-1:1
    pii0 = 1/MDL;
    pik = zeros(1,MDL);
    for i = 1:1:MDL
        pik(i) = 1/MDL;
    end
    
    R = zeros(2,2,MDL);
    for i = 1:MDL
        R(:,:,i) = [1 0 ; 0 1];
    end
    
    u = zeros(MDL,2);
    for i = 1:MDL
        u(i,:) = x(i,:);
    end
    
    Size = size(x);
    length = Size(1);
    
    pn = zeros(1,MDL);
    iterations = 20;
    
    for iter = 1: iterations
        
        Nk = zeros(1,MDL);
        
        t1k = zeros(2,MDL);
        
        t2k = zeros(2,2,MDL);
        for n = 1:1:length
            for i = 1:1:MDL
                pn(i) =pik(i)*(det(R(:,:,i))^(-0.5)).*exp(-0.5*(x(n,:)-u(i,:))*inv(R(:,:,i))*(x(n,:)-u(i,:))')./(2*pi)^(M/2); %#ok<MINV>
            end
            pd = sum(pn);
            
            Nk = Nk + pn/pd;
            
            t1k = t1k + ((pn/pd)'*x(n,:))';
            
            for i = 1:1:MDL
                t2k(:,:,i) = t2k(:,:,i)  + (pn(i)/pd)'*(x(n,:)'*x(n,:));
            end
        end
        fprintf('Iteration = %d \n',iter)
        
        for i = 1:1:MDL
            u(i,:) = (t1k(:,i)./Nk(i))'
            
            R(:,:,i) = t2k(:,:,i)./Nk(i)-t1k(:,i)*t1k(:,i)'/(Nk(i).^2)
            
            pik(i) = Nk(i)./length
        end
        kk =1;
        
        resul(MDL,iter) = fun_MDL(MDL,x,pik,u,R);
    end
    resul_K(MDL) = min(resul(:,iter));
    [l, m, pilm, ulm, Rlm] = fun_d(length, MDL, pik, u, R);
    pik(l) = pilm;
    pik(m) = pik(MDL);
    u(l,:) = ulm;
    u(m,:) = u(MDL,:);
    R(:,:,l) = Rlm;
    R(:,:,m) = R(:,:,MDL);
    if iter < 20
        resul(MDL,iter) = fun_MDL(MDL, x, pik, u, R);
    end
end
figure(1)
plot(resul,'o-'); grid on;
xlabel('K'); ylabel('MDL');

figure(2)
surf(resul); grid on;
xlabel('Iteration');ylabel('K'); zlabel('MDL');
axis([1 iterations 1 K  min(min(resul)) max(max(resul))])
