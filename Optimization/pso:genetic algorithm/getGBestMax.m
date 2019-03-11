function [gbestcor,gbest,gworstcor,gworst,...
    gaveragecor,gaverage] = getGBestMax(xx)
[~,b] = size(xx);
gg = zeros(1,b);
for i = 1:b
    gg(1,i) = rastrigins(xx(:,i)');
end

gworstcor = xx(:,find(gg==min(gg),1));
gworst =  min(gg);
gbestcor = xx(:,find(gg==max(gg),1));
gbest = max(gg);
gaveragecor = xx(:,find(gg==...
    min(gg-mean(gg))+mean(gg),1));
gaverage =  mean(gg);
end

