function [ppbest] = getPBestMax(xx,iter)
[~,b,~] = size(xx);
pp = zeros(2,b);
ppbest = zeros(2,b);
for i = 1:b
        pp(:,i) = rastrigins(squeeze(xx(:,i,iter:iter+1))');
end

for i = 1:b
    
if pp(1,i) > pp(2,i)
    ppbest(:,i) = pbestMinTrack(:,:,iter);
else ppbest(:,i) = xx(:,i,iter+1);
end
end
end

