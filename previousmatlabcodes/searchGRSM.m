function [registers] = searchGRSM(grp) 
  
%% load the database 
load('grsmDB.mat'); 
  
%number of states in the GRSM 
states = 5; 
  
%maximimum number of genes in any GRSM from the database 
grsmMax = 14; 
  
%% ensure that the number of states in the grp is compatible with the database 
if size(grp,1) ~= states 
    disp('error: number of states in grp incompatible with database'); 
    return; 
end 
  
%% if the number of genes in the grp exceeds the max in the database, return an error 
if size(grp,2) > grsmMax 
    disp('error: number of genes in grp exceeds max amount in the database'); 
    return; 
end 
  
%% if the number of genes in the grp is less than grsmMax 
%then add the extra genes to the grp (but specify 
%them as OFF (0) in every state) 
if size(grp,2) < grsmMax 
    grp = [grp zeros(states,grsmMax - size(grp,2))]; 
end 
  
%% find all solutions that match the grp: 
%put grp into same format as in the database (i.e. a 70-element vector); 
grp = sortrows(grp'); 
grp = grp(:); 
  
%search for grp, store indeces of matches 
[~,grpidx]=ismember(grp',grpArray','rows'); 
  
if grpidx==0 %if no match,return message and empty registers 
    disp('no matching registers for input gene regulation program') 
    registers = []; 
     
else %if there are matches, then pull the list of registers 
    registers = registerArray(register2grp==grpidx,:); 
    
     
    %rank the registers 
    registers = registerRank(registers);    
end