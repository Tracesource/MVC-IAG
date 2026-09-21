function [y_z,A,S,iter,obj] = ARMVC_new1(X,numanchor,beta,gamma)

%% initialize           
maxIter = 50 ; % the number of iterations
IterMax = 50;

m = numanchor;
numview = length(X);
numsample = size(X{1},1);

XX = [];
for p = 1 : numview
    X{p} = mapstd(X{p}',0,1);
    XX = [XX;X{p}];
end

rand('twister',12);
[IDX,AC] = kmeans(XX',m, 'MaxIter',100,'Replicates',10);

sigma = 0.3*((m-1)*0.5-1)+0.8;
S =zeros(m,m);
Gall=zeros(m,m);
count = 1;
for i = 1:numview
    di = size(X{i},1); 
    Z{i} =zeros(m,numsample);
    A{i} = AC(:,count:count+di-1)';
    XA{i}=X{i}'*A{i};
    AA{i} = A{i}'*A{i};
    count = count+di;
    D = eucl_dist2(A{i}',A{i}');
    G{i} = exp(-D/(2*sigma^2));
    G{i} = G{i}-eye(m);
    G{i} = G{i}./sum(G{i},2)';
    Gall = Gall+G{i};
end
S= S/numview;
Ls=(diag(sum((S+S')/2))-(S+S')/2);

flag = 1;
iter = 0;
%%
while flag
    iter = iter + 1;

    %% optimize Zv
    options = optimset( 'Algorithm','interior-point-convex','Display','off');
    parfor iv=1:numview
        Hv = AA{iv}+beta*Ls;
        for ii=1:numsample
            Z{iv}(:,ii) = quadprog(Hv,-XA{iv}(ii,:)',[],[],ones(1,m),1,zeros(m,1),ones(m,1),[],options);
        end
        ZZ{iv} = Z{iv}*Z{iv}';
    end
    
    
    %% optimize Si
    theta=0;
    for iv=1:numview
        theta = theta-beta/(2*gamma).*(diag(diag(ZZ{iv}))*ones(m,m)+ones(m,m)*diag(diag(ZZ{iv}))-2*ZZ{iv});
    end
    theta = (Gall+theta)/numview;
    for ii=1:m
        S([1:ii-1,ii+1:end],ii) = EProjSimplex_new(theta([1:ii-1,ii+1:end],ii));
    end
    Ls=(diag(sum((S+S')/2))-(S+S')/2);
    
    term1 = 0;
    term2 = 0;
    term3 = 0;
    for iv = 1:numview
        term1 = term1 + norm(X{iv}-A{iv}*Z{iv},'fro')^2;
        term2 = term2 + trace(Ls*ZZ{iv}');
        term3 = term3 + norm(S-G{iv},'fro')^2;
    end
    
    obj(iter) = term1+beta*term2+gamma*term3;
    
	if (iter>1) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-3 || iter>maxIter || obj(iter) < 1e-10)
        flag = 0;
        S = (S+S')/2;
        GSall = 0;
        for iv=1:numview
           G{iv}=(G{iv}+G{iv}')/2;
           aa = sort(G{iv}(:),'ascend');
           Gv10 = mean(aa(m+1:m+floor(0.9*(m^2-m))));
           Gvsort = (G{iv}<Gv10);
           GS{iv}= ones(m,m);
           GS{iv}(Gvsort)=0;
           GSall = GSall+GS{iv};
        end
        Gall = Gall/(numview);
        aa = sort(Gall(:),'ascend');
        average_top10 = mean(aa(m+1:m+floor(0.9*(m^2-m))));
        Ssort = (S<average_top10);
        S(Ssort)=0;
        S = S.*GSall;
        [~,y_s]=conncomp(sparse(S));
        
        
        Zlabel = full(sparse(1:numsample, IDX, 1, numsample, m));
        y_z = y_s*Zlabel';
    end
end
         
         
    
