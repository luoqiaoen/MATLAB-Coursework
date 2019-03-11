%% PCA decomposition of alphabet using different Fonts (Section 1 to 5)
clc;
% clear;
close ALL;
%% Section 1

mu = [0 0];
UNITY = [1 0; 0 1];
SIGMA = [2 -1.2; -1.2 1];

w = mvnrnd(mu,UNITY,1000);
w = w.';
figure(1)
plot(w(1,:),w(2,:),'+')
axis equal


[V,D] = eig(SIGMA);
x_bar = (D.^0.5)*w;

figure(2)
plot(x_bar(1,:),x_bar(2,:),'+')
axis equal

x = V*x_bar;
figure(3)
plot(x(1,:),x(2,:),'.')
axis equal
% Section 2

mu_2 = mean(x,2);
xx(1,:) = x(1,:)- mu_2(1);
xx(2,:) = x(2,:)- mu_2(2);

cov_2 = xx*(xx.')/999;
[V_2,D_2] = eig(cov_2);
xx_bar = inv(V_2)*x;
w_2 = inv(D_2.^0.5)*xx_bar;
figure(4)
plot(xx_bar(1,:),xx_bar(2,:),'+')
axis equal
figure(5)
plot(w_2(1,:),w_2(2,:),'+')
axis equal

mu_w2 = mean(w_2,2);
ww(1,:) = w_2(1,:)- mu_w2(1);
ww(2,:) = w_2(2,:)- mu_w2(2);

cov_w2 = ww*(ww.')/999;

%% Section 4
%run read_data.m in training_data

mu_3 = mean(X,2);
mu_3 = repmat(mu_3,1,312);
Z_3 = (X-mu_3)./64;
[U_3,S_3,V_3] = svd(Z_3,0);

sorting=reshape(S_3,[1,97344]);
sorting= sort(sorting,'descend');

twelve(1:12)= sorting(1:12);
[row,col] = find(S_3>=sorting(12));

for k=1:12
    img=reshape(U_3(:,k),[64,64]);
    figure(6); subplot(3,4,k); imagesc(img);
    axis('image'); colormap(gray(256));
end

Y = U_3.'*(X-mu_3);



figure(7)
plot(1:10,Y(1:10,1),'r',1:10,Y(1:10,2),'b',1:10,Y(1:10,3),'g',1:10,Y(1:10,4),'m');
legend('first 10 projection coefficients for a','first 10 projection coefficients for b','first 10 projection coefficients for c','first 10 projection coefficients for d')

figure(8)
img=reshape(X(:,1),[64,64]);
imagesc(img);
colormap(gray(256));

k=1;
for m=[1,5:5:20,30];
    U_zero = zeros(size(U_3));
    U_zero(:,1:m)=U_3(:,1:m);
    X_bar = U_zero*Y;
    figure(9);subplot(2,3,k);
    img=reshape(X_bar(:,1),[64,64]);
    imagesc(img);
    axis('image'); colormap(gray(256));
    k=k+1;
end


%% Using Eigendecomposition to classify test data (not very good)
if 1
    % % Section 5
    % run read_data.m in test_data
    Y_5=U_3(:,1:10).'*(X-mu_3);
    
    empty_cell=cell(26,2);
    params=cell2struct(empty_cell,{'M','R'},2);
    
    
    i = 1;
    
    
    while i < 27
        params(i).M = mean(Y_5(:,i:26:286+i),2);
        params(i).R = (Y_5(:,i:26:286+i)-repmat(mean(Y_5(:,i:26:286+i),2),1,12))*(Y_5(:,i:26:286+i)-repmat(mean(Y_5(:,i:26:286+i),2),1,12)).'/9;
        i = i+1;
    end
    
    YY = U_3(:,1:10).'*(XX-repmat(mean(X,2),1,26));
    m=zeros(26,1);
    n=zeros(26,1);
    
    for j=1:26
        YY_zero=zeros(size(YY));
        YY_zero(:,j)=YY(:,j);
        argmin = ((YY_zero-repmat(params(j).M,1,26)).'*inv(params(j).R)*(YY_zero-repmat(params(j).M,1,26)))+log((((YY_zero-repmat(params(j).M,1,26)).'*inv(params(j).R)*(YY_zero-repmat(params(j).M,1,26)))*((YY_zero-repmat(params(j).M,1,26)).'*inv(params(j).R)*(YY_zero-repmat(params(j).M,1,26))).').^(1/2));
        [maxA,ind] = min(argmin(:));
        [m(j,1),n(j,1)] = ind2sub(size(argmin),ind);
    end
    
    
    for j=1:26
        [V,D] = eig(params(j).R);
        YY_zero=zeros(size(YY));
        YY_zero(:,j)=YY(:,j);
        argmin = ((YY_zero-repmat(params(j).M,1,26)).'*inv(D)*(YY_zero-repmat(params(j).M,1,26)))++log((((YY_zero-repmat(params(j).M,1,26)).'*inv(D)*(YY_zero-repmat(params(j).M,1,26)))*((YY_zero-repmat(params(j).M,1,26)).'*inv(D)*(YY_zero-repmat(params(j).M,1,26))).').^(1/2));
        [maxA,ind] = min(argmin(:));
        [m(j,1),n(j,1)] = ind2sub(size(argmin),ind);
    end
    
    
    
    
    j=1;
    adding=zeros(10,10);
    while j<27
        adding=adding+params(j).R;
        j=j+1;
    end
    
    for j=1:26
        YY_zero=zeros(size(YY));
        YY_zero(:,j)=YY(:,j);
        argmin = ((YY_zero-repmat(params(j).M,1,26)).'*inv(adding./26)*(YY_zero-repmat(params(j).M,1,26)))++log((((YY_zero-repmat(params(j).M,1,26)).'*inv(adding./26)*(YY_zero-repmat(params(j).M,1,26)))*((YY_zero-repmat(params(j).M,1,26)).'*inv(adding./26)*(YY_zero-repmat(params(j).M,1,26))).').^(1/2));
        
        [maxA,ind] = min(argmin(:));
        [m(j,1),n(j,1)] = ind2sub(size(argmin),ind);
    end
    
    for j=1:26
        YY_zero=zeros(size(YY));
        YY_zero(:,j)=YY(:,j);
        [V,D] = eig(adding./26);
        argmin = ((YY_zero-repmat(params(j).M,1,26)).'*inv(D)*(YY_zero-repmat(params(j).M,1,26)))++log((((YY_zero-repmat(params(j).M,1,26)).'*inv(D)*(YY_zero-repmat(params(j).M,1,26)))*((YY_zero-repmat(params(j).M,1,26)).'*inv(D)*(YY_zero-repmat(params(j).M,1,26))).').^(1/2));
        
        [maxA,ind] = min(argmin(:));
        [m(j,1),n(j,1)] = ind2sub(size(argmin),ind);
    end
    
    for j=1:26
        YY_zero=zeros(size(YY));
        YY_zero(:,j)=YY(:,j);
        argmin = ((YY_zero-repmat(params(j).M,1,26)).'*inv(eye(10))*(YY_zero-repmat(params(j).M,1,26)))+log((((YY_zero-repmat(params(j).M,1,26)).'*inv(eye(10))*(YY_zero-repmat(params(j).M,1,26)))*((YY_zero-repmat(params(j).M,1,26)).'*inv(eye(10))*(YY_zero-repmat(params(j).M,1,26))).').^(1/2));
        [maxA,ind] = min(argmin(:));
        [m(j,1),n(j,1)] = ind2sub(size(argmin),ind);
    end
end