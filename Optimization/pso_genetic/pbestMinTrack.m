function [ppbest] = getPBestMin(xx,iter)
[~,b,~] = size(xx);
pp = zeros(2,b);
for i = 1:b
        pp(:,i) = rastrigins(squeeze(xx(:,i,iter:iter+1))');
end

if pp(1,:) < pp(2,:)
    ppbest = xx(:,:,iter);
else ppbest = xx(:,:,iter+1);
end
end
