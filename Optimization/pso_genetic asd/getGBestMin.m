function [gbestcor,gbest,gworstcor,gworst,...
    gaveragecor,gaverage] = getGBestMin(xx)
[~,b] = size(xx);
gg = zeros(1,b);
for i = 1:b
    gg(1,i) = rastrigins(xx(:,i)');
end

gbestcor = xx(:,find(gg==min(gg),1));
gbest =  min(gg);
gworstcor = xx(:,find(gg==max(gg),1));
gworst =  max(gg);
gaveragecor = xx(:,find(gg==...
    min(gg-mean(gg))+mean(gg),1));
gaverage =  mean(gg);
end

