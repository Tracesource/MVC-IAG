
clear;
clc;
warning off;
addpath(genpath('./'));

%% dataset
ds = {'Dermatology'};
dsPath = '.\datasets\';
resPath = './res-lmd/';


for dsi =1:1:length(ds)
    dataName = ds{dsi}; disp(dataName);
    load(strcat(dsPath,dataName));
    n = length(Y);
   realcluster = length(unique(Y));
   
   txtpathmean = strcat(resPath,strcat(dataName,'_mean.txt'));
    dlmwrite(txtpathmean, strcat('Dataset:',cellstr(dataName), '  Date:',datestr(now)),'-append','delimiter','','newline','pc');
    
    anchor = [10 20 30 40 50];
    beta = 10.^[-2 0 2 4];
    gamma = 10.^[-2 0 2 4];
    
    %%
    allresult = [];
    for ichor = 1:length(anchor)
        for ib = 1:length(beta)
            for ig = 1:length(gamma)
                tic;
                [y1,A,S,iter,obj] = MVC-IAG(X,anchor(ichor),beta(ib),gamma(ig));
                testcluster = length(unique(y1));
                res = Clustering8Measure(y1, Y);
                timer(ichor,ib,ig)  = toc;
                fprintf('Anchor:%d \t Beta:%d\t Gamma:%d\t Cluster:%d\t %d\t Res:%12.6f %12.6f %12.6f %12.6f \tTime:%12.6f \n',[anchor(ichor) beta(ib) gamma(ig) realcluster testcluster res(1) res(2) res(3) res(4) timer(ichor,ib,ig)]);
                dlmwrite(txtpathmean,  [anchor(ichor) beta(ib) gamma(ig) realcluster testcluster res timer(ichor,ib,ig)],'-append','delimiter','\t','newline','pc');
                allresult = [allresult;res realcluster testcluster  timer(ichor,ib,ig)];
            end
        end
    end
    [c,d] = max(allresult(:,1));
    maxresult = allresult(d,:);
    dlmwrite('./totalResults/AllDatasetResult.txt',char(dataName),'-append','delimiter','\t','newline','pc');
    dlmwrite('./totalResults/AllDatasetResult.txt',maxresult,'-append','delimiter','\t','newline','pc');
end


