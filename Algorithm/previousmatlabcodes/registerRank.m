function [rankedsoln] = registerRank(soln) 
  
%% we work mainly with the absolute value of the registers now because we are 
%mainly checking for the presence of parts and don't mind their directions. 
absoln = abs(soln); 
  
%% we want to minimize the number of elements that have promoter readthrough.  
prtVec = sum(ismember(absoln,[9 10 11 12 22 23 24]) + 2*ismember(absoln,[13 25]),2); 
[prtSort,prtSortIdx]=sort(prtVec); 
absoln = absoln(prtSortIdx,:); 
soln = soln(prtSortIdx,:); 
  
%% perform other subrankings by priority  
for i = unique(prtSort)'; 
    %minimize terminator readthrough 
    idxi = find(prtSort==i); 
    pVec = sum(ismember(absoln(idxi,:),7),2); 
    [pSort,pSortIdx]=sort(pVec); 
  
    absoln(idxi,:) = absoln(idxi(pSortIdx),:); 
    soln(idxi,:) = soln(idxi(pSortIdx),:); 
    for j = unique(pSort)'; 
        %maximize number of blank spaces 
        idxj = find(pSort==j); 
        eVec = sum(ismember(absoln(idxi(idxj),:),5),2); 
      
        [eSort,eSortIdx]=sort(eVec,'descend'); 
        absoln(idxi(idxj),:)=absoln(idxi(idxj(eSortIdx)),:); 
        soln(idxi(idxj),:)=soln(idxi(idxj(eSortIdx)),:); 
        for k = unique(eSort)'; 
            %minimize number of promoters 
            idxk = find(eSort==k); 
            tVec = sum(ismember(absoln(idxi(idxj(idxk)),:),[2 3 9 10 12 14 16 19 21]) + 2*ismember(absoln(idxi(idxj(idxk)),:),[6 11 13 17 18 22 24]) + 3*ismember(absoln(idxi(idxj(idxk)),:),[20 23 25]),2); 
             
            [tSort,tSortIdx]=sort(tVec); 
            absoln(idxi(idxj(idxk)),:)=absoln(idxi(idxj(idxk(tSortIdx))),:); 
            soln(idxi(idxj(idxk)),:)=soln(idxi(idxj(idxk(tSortIdx))),:); 
        end 
    end 
end 
  
rankedsoln=soln; 