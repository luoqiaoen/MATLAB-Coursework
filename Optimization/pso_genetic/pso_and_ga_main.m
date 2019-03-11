clear
format long

psomin = 1;
psomax = 0;
cga = 0;
nga = 0;
tsp = 0;
%% plot surf and contour
x1 = linspace(-10,10,200);
x2 = linspace(-10,10,200);
[XX1,XX2]=meshgrid(x1,x2);
Z = 20+(XX1./10).^2+(XX2./10).^2 - 10*(cos(2*pi*XX1/10) +cos(2*pi*XX2/10));
figure(1)
zmin = floor(min(Z(:)));
zmax = ceil(max(Z(:)));
zinc = (zmax - zmin) / 40;
zlevs = zmin:zinc:zmax;
contour(XX1,XX2,Z,zlevs,'ShowText','on')
xlabel('X')
ylabel('Y')
title('global best evolution')
axis image
colorbar


%% PSO Minimum
% generate initial guess
if psomin
    c1 = 2.01;
    c2 = 2.09;
    w = 0.729; %inertial constant
    vmax = 4.0;
    iter = 150;
    n = 50;
    xx = 20*(rand(2,n)-0.5);
    pp = xx;
    vv = 2*(rand(2,n)-0.5);
    
    [gbestcor,gbest,gworstcor,gworst,...
        gaveragecor,gaverage]  = getGBestMin(xx);
    
    xxTrack = zeros(2,n,iter+1);
    xxTrack(:,:,1) = xx;
    vvTrack = zeros(2,n,iter+1);
    vvTrack(:,:,1) = vv;
    gbestMinTrack = zeros(1,iter+1);
    gbestMinTrack(1)=gbest;
    gbestCorMinTrack = zeros(2,iter+1);
    gbestCorMinTrack(:,1) = gbestcor;
    gaverageMinTrack = zeros(1,iter+1);
    gaverageMinTrack(1)=gaverage;
    gworstMinTrack = zeros(2,iter+1);
    gworstMinTrack(1)=gworst;
    pbestMinTrack = zeros(2,n,iter+1);
    pbestMinTrack(:,:,1)= pp;
    pptemp = zeros(2,n);
    rr = rand(2,n);
    ss = rand(2,n);
    
    for i = 1:iter-1
        rr = rand(2,n);
        ss = rand(2,n);
        vvTrack(:,:,i+1) = w*(vvTrack(:,:,i) + c1*rr.*(pbestMinTrack(:,:,i)...
            -xxTrack(:,:,i))) + c2*ss.*(repmat(gbestCorMinTrack(:,i),1,n)-xxTrack(:,:,i));
        vvTrack(1,:,i+1) = min(vmax,max(-vmax,vvTrack(1,:,i+1)));
        vvTrack(2,:,i+1) = min(vmax,max(-vmax,vvTrack(2,:,i+1)));
        
        xxTrack(:,:,i+1) = xxTrack(:,:,i)+vvTrack(:,:,i+1);
        %        xxTrack(1,:,i+1) = min(xmax,max(-xmax,xxTrack(1,:,i+1)));
        %        xxTrack(2,:,i+1) = min(xmax,max(-xmax,xxTrack(2,:,i+1)));
        for jj = 1:n
            if rastrigins(xxTrack(:,jj,i+1)')< rastrigins(xxTrack(:,jj,i)')
                pbestMinTrack(:,jj,i+1) = xxTrack(:,jj,i+1);
            else pbestMinTrack(:,jj,i+1) = pbestMinTrack(:,jj,i);
            end
        end
        
%         pbestMinTrack(:,:,i+1) = getPBestMin(xxTrack,i,pbestMinTrack);
        
        [gbestcor,gbest,gworstcor,gworst,...
            gaveragecor,gaverage]  = getGBestMin(xxTrack(:,:,i+1));
        if gbest < rastrigins(gbestCorMinTrack(:,i)')
            gbestMinTrack(i+1) = gbest;
            gbestCorMinTrack(:,i+1) = gbestcor;
            gworstMinTrack(i+1)= gworst;
            gaverageMinTrack(i+1) = gaverage;
        else
            gbestMinTrack(i+1) = gbestMinTrack(i);
            gbestCorMinTrack(:,i+1) =  gbestCorMinTrack(:,i);
            gworstMinTrack(i+1)= gworstMinTrack(i);
            gaverageMinTrack(i+1) = gaverageMinTrack(i);
            
        end
    end
    
    figure(1)
    hold on
    plot(gbestCorMinTrack(1,:),gbestCorMinTrack(2,:),'k+-')
    hold off
    
    
    finish = rastrigins(gbestMinTrack(:,end)');
    figure(4)
    plot(1:i+1, gbestMinTrack(1:i+1),'k+-',1:i+1, gworstMinTrack(1:i+1),'bo-',1:i+1, gaverageMinTrack(1:i+1),'g*-')
    axis image
    ylabel('Object Function Value');
    xlabel('Number of Iteration');
    legend('global best','global worst','global average');
end


%% PSO Maximum
% generate initial guess
if psomax
    c1 = 2.01;
    c2 = 2.09;
    w = 0.729; %inertial constant
    vmax = 4.0;
    iter = 150;
    n = 50;
    xx = 20*(rand(2,n)-0.5);
    pp = xx;
    vv = 0.5*(rand(2,n)-0.5);
    [gbestcor,gbest,gworstcor,gworst,...
        gaveragecor,gaverage]  = getGBestMax(xx);
    
    xxTrack = zeros(2,n,iter+1);
    xxTrack(:,:,1) = xx;
    vvTrack = zeros(2,n,iter+1);
    vvTrack(:,:,1) = vv;
    gbestMaxTrack = zeros(1,iter+1);
    gbestMaxTrack(1)=gbest;
    gbestCorMaxTrack = zeros(2,iter+1);
    gbestCorMaxTrack(:,1) = gbestcor;
    gaverageMaxTrack = zeros(1,iter+1);
    gaverageMaxTrack(1)=gaverage;
    gworstMaxTrack = zeros(2,iter+1);
    gworstMaxTrack(1)=gworst;
    pbestMaxTrack(:,:,1)= pp;
    pptemp = zeros(2,n);
    
    for i = 1:iter-1
        rr = rand(2,n);
        ss = rand(2,n);
        vvTrack(:,:,i+1) = w*(vvTrack(:,:,i) + c1*rr.*(pbestMaxTrack(:,:,i)...
            -xxTrack(:,:,i))) + c2*ss.*(repmat(gbestCorMaxTrack(:,i),1,n)-xxTrack(:,:,i));
        vvTrack(1,:,i+1) = min(vmax,max(-vmax,vvTrack(1,:,i+1)));
        vvTrack(2,:,i+1) = min(vmax,max(-vmax,vvTrack(2,:,i+1)));
        
        xxTrack(:,:,i+1) = xxTrack(:,:,i)+vvTrack(:,:,i+1);
        for jj = 1:n
            if rastrigins(xxTrack(:,jj,i+1)')> rastrigins(xxTrack(:,jj,i)')
                pbestMaxTrack(:,jj,i+1) = xxTrack(:,jj,i+1);
            else pbestMaxTrack(:,jj,i+1) = pbestMaxTrack(:,jj,i);
            end
        end
       % pbestMaxTrack(:,:,i+1) = getPBestMax(xxTrack,i);
        
        [gbestcor,gbest,gworstcor,gworst,...
            gaveragecor,gaverage]  = getGBestMax(xxTrack(:,:,i+1));
        if gbest > rastrigins(gbestCorMaxTrack(:,i)')
            gbestMaxTrack(i+1) = gbest;
            gbestCorMaxTrack(:,i+1) = gbestcor;
            gworstMaxTrack(i+1)= gworst;
            gaverageMaxTrack(i+1) = gaverage;
        else
            gbestMaxTrack(i+1) = gbestMaxTrack(i);
            gbestCorMaxTrack(:,i+1) =  gbestCorMaxTrack(:,i);
            gworstMaxTrack(i+1)= gworstMaxTrack(i);
            gaverageMaxTrack(i+1) = gaverageMaxTrack(i);
        end
    end
    
    figure(1)
    hold on
    plot(gbestCorMaxTrack(1,:),gbestCorMaxTrack(2,:),'k+-')
    hold off
    
    finish = rastrigins(gbestMaxTrack(:,end)');
    figure(4)
    plot(1:i+1, gbestMaxTrack(1:i+1),'k+-',1:i+1, gworstMaxTrack(1:i+1),'bo-',1:i+1, gaverageMaxTrack(1:i+1),'g*-')
    axis image
    ylabel('Object Function Value');
    xlabel('Number of Iteration');
    legend('global best','global worst','global average')
end

%% Canonical GA
% GA parameters
if cga
    iter = 150;
    n = 200;
    bitnumber = 16;
    maxi = 2^(16-1);
    xx_bin = randi(maxi,2,n);
    xx_dec = mapping(xx_bin, 16);
    pc = 0.85;
    pm = 0.15;
    count = 0;
    gbestTrack = zeros(1,iter+1);
    gworstTrack = zeros(1,iter+1);
    gaverageTrack = zeros(1,iter+1);
    [gbestcor,gbest,gworstcor,gworst,...
        gaveragecor,gaverage] = getGBestMin(xx_dec);
    gbestTrack(1) = gbest;
    gworstTrack(1) = gworst;
    gaverageTrack(1) = gaverage;
    for counter = 1:iter
        %% Selection and crossover
        prod = rastrigins(xx_dec');
        base = max(prod);
        fit = base-prod;
        if max(fit)<10^(-10);
            break
        else fit = fit/sum(fit);
            cumSum = cumsum(fit);
        end
        
        index1 = zeros(1,n/2);
        index2 = zeros(1,n/2);
        xx_dec_new1 = zeros(2,n/2);
        xx_dec_new2 = zeros(2,n/2);
        store1 = zeros(2,n/2,bitnumber);
        store2 = zeros(2,n/2,bitnumber);
        decode1 = zeros(2,n/2);
        decode2 = zeros(2,n/2);
        for i = 1:n/2
            index1(i) = find(cumSum-rand>0, 1 );
            index2(i) = find(cumSum-rand>0, 1 );
            xx_dec_new1(:,i) = xx_dec(:,index1(i));
            xx_dec_new2(:,i) = xx_dec(:,index2(i));
            parent1(1,:) = num2bit(xx_dec_new1(1,i),16,-10,10);
            parent1(2,:) = num2bit(xx_dec_new1(2,i),16,-10,10);
            parent2(1,:) = num2bit(xx_dec_new2(1,i),16,-10,10);
            parent2(2,:) = num2bit(xx_dec_new2(2,i),16,-10,10);
            % cross over
            if  rand < pc
                crosspoint = ceil(rand*(bitnumber-1));
                child1 = [parent1(:,1:crosspoint),...
                    parent2(:,crosspoint+1:end)];
                child2 = [parent2(:,1:crosspoint),...
                    parent1(:,crosspoint+1:end)];
            else
                child1 = parent1;
                child2 = parent2;
            end
            store1(:,i,:) = child1;
            store2(:,i,:) = child2;
            % mutation operator
            if rand < pm
                store1(:,i,:) = abs(store1(:,i,:)-1);
                store1(:,i,:) = abs(store1(:,i,:)-1);
            end
        end
        
        
        for i = 1:n/2
            decode1(1,i) = bit2num(store1(1,i,:),16,-10,10);
            decode1(2,i) = bit2num(store1(2,i,:),16,-10,10);
            decode2(1,i) = bit2num(store2(1,i,:),16,-10,10);
            decode2(2,i) = bit2num(store2(2,i,:),16,-10,10);
        end
        xx_dec(:,1:n/2)=decode1;
        xx_dec(:,n/2+1:end)=decode2;
        [gbestcor,gbest,gworstcor,gworst,...
            gaveragecor,gaverage] = getGBestMin(xx_dec);
        gbestTrack(counter+1) = gbest;
        gworstTrack(counter+1) = gworst;
        gaverageTrack(counter+1) = gaverage;
        count = count +1;
    end
    figure(1)
    hold on
    plot(xx_dec(1,:),xx_dec(2,:),'g*-')
    hold off
    figure(2)
    plot(1:counter+1, gbestTrack(1:counter+1),'k+-',1:counter+1, gworstTrack(1:counter+1),'bo-',1:counter+1, gaverageTrack(1:counter+1),'g*-')
    axis image
    ylabel('Object Function Value');
    xlabel('Number of Iteration');
    legend('global best','global worst','global average')
end


%% Numerical GA
% GA parameters
if nga
    iter = 100;
    n = 200;
    bitnumber = 16;
    maxi = 2^(16-1);
    xx_bin = randi(maxi,2,n);
    xx_dec = mapping(xx_bin, 16);
    pc = 0.75;
    pm = 0.1;
    count = 0;
    gbestTrack = zeros(1,iter+1);
    gworstTrack = zeros(1,iter+1);
    gaverageTrack = zeros(1,iter+1);
    [gbestcor,gbest,gworstcor,gworst,...
        gaveragecor,gaverage] = getGBestMin(xx_dec);
    gbestTrack(1) = gbest;
    gworstTrack(1) = gworst;
    gaverageTrack(1) = gaverage;
    for counter = 1:iter
        %% Selection and crossover
        prod = rastrigins(xx_dec');
        base = max(prod);
        fit = base-prod;
        if max(fit)<10^(-10);
            break
        else fit = fit/sum(fit);
            cumSum = cumsum(fit);
        end
        %     index1 = zeros(1,n/2);
        %     index2 = zeros(1,n/2);
        %     xx_dec_new1 = zeros(2,n/2);
        %     xx_dec_new2 = zeros(2,n/2);
        %     store1 = zeros(2,n/2);
        %     store2 = zeros(2,n/2);
        index1 = zeros(1,n);
        index2 = zeros(1,n);
        xx_dec_new1 = zeros(2,n/2);
        xx_dec_new2 = zeros(2,n/2);
        store1 = zeros(2,n);
        store2 = zeros(2,n/2);
        
        for i = 1:n
            index1(i) = find(cumSum-rand>0, 1 );
            index2(i) = find(cumSum-rand>0, 1 );
            xx_dec_new1(:,i) = xx_dec(:,index1(i));
            xx_dec_new2(:,i) = xx_dec(:,index2(i));
            if rand < pc
                parent1 = xx_dec_new1(:,i);
                parent2 = xx_dec_new2(:,i);
                % cross over
                child1 = (parent1+parent2)/2 + ...
                    normrnd(0, 0.01);
                %         child2 = (parent1+parent2)/2 + ...
                %             normrnd(0, 0.2)*0.7^(i);
                store1(:,i) = child1;
                %store2(:,i) = child2;
            end
            % mutation operator
            if rand < pm
                store1(:,i) = store1(:,i) +normrnd(0,0.01);
                %store2(:,i) = store2(:,i) +normrnd(0,0.1)*0.9^(i);
            end
        end
        
        
        %     xx_dec(:,1:n/2)=store1;
        %     xx_dec(:,n/2+1:end)=store2;
        xx_dec = store1;
        [gbestcor,gbest,gworstcor,gworst,...
            gaveragecor,gaverage] = getGBestMin(xx_dec);
        gbestTrack(counter+1) = gbest;
        gworstTrack(counter+1) = gworst;
        gaverageTrack(counter+1) = gaverage;
        count = count +1;
    end
    figure(1)
    hold on
    plot(xx_dec(1,:),xx_dec(2,:),'g*-')
    hold off
    figure(2)
    plot(1:counter+1, gbestTrack(1:counter+1),'k+-',1:counter+1, gworstTrack(1:counter+1),'bo-',1:counter+1, gaverageTrack(1:counter+1),'g*-')
    axis image
    ylabel('Object Function Value');
    xlabel('Number of Iteration');
    legend('global best','global worst','global average')
end


%% Traveling salesman problem
if tsp
    k = 200;
    d = 200;
    pc = .9;%prob_crossing
    makesure =5;
    p= [7.7176 0.8955
        9.8397 9.8352
        1.3001 7.7070
        4.6141 0.0784
        9.0204 4.8838
        1.5651 3.3703
        6.4311 6.1604
        8.0298 5.5397
        8.2941 1.2147
        0.6621 2.5383
        9.6546 0.5923
        6.8045 9.2022
        8.0196 2.8386
        8.1590 7.1041
        8.0531 1.8498];
    cities = transpose(p);
    cities_n = size(cities,2);
    order = zeros(d, cities_n);
    m = zeros(d,cities_n);
    fit = zeros(d,1);
    
    for i = 1:d
        order(i,:) = randperm(cities_n);
        fit(i) = getDist(cities, order(i,:));
    end
    
    
    for kk = 1:k
        [value, index] = sort(fit(:),'descend');
        for i = 1:makesure
            
            m(i,:) = order(index(i),:);
        end
        
        for i = makesure + 1:d
            j = find(cumsum(fit)/sum(fit) - rand()>0,1);
            m(i,:) = order(j,:);
            
            if(rand() <pc)
                randA = randi([1 cities_n]);
                randB = randi([1, cities_n]);
                
                temp = m(i,randA);
                m(i,randA) = m(i,randB);
                m(i,randB) = temp;
            end
        end
        order = m;
        for i = 1:d
            fit(i) = getDist(cities, order(i,:));
        end
    end
    
    
    bestRoute = order(find(fit ==max(fit),1),:);
    bextDist  = 1/getDist(cities, bestRoute);
    x = cities(1,bestRoute);
    x(cities_n +1) = cities(1,bestRoute(1));
    y = cities(2,bestRoute);
    y(cities_n +1) = cities(2,bestRoute(1));
    figure
    plot(x,y,'b-o')
end

