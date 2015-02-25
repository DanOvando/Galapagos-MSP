%% mjoin
% Dan Ovando
% 11/15/2011

%mjoin takes matrix b and merges it to matrix a. b must have unique entries
%for every key value of a. a and b should be organized
%such that data are stored in columns, entries in rows
%akey marks the column location of
%the key variables in a, bkey marks the location in b. get is a vector of
%column indices that show the data that should be taken from b and added
%to a. the final product tacks the merged data from b on to the end of the
%matrix a. Key valus in a with no matches in b receive an NaN in the final
%merged matrix




function [merged newdata]=mjoin(a,b,akey,bkey,get)

% load method
% load interpcount
% 
% a=method
% b=interpcount
% 
% b(:,5)=b(:,4)>.1
% 
% 
% akey=1
% bkey=1
% get=5

tempmerge=nan(size(a,1),length(get));

avals=unique(a(:,akey));

for i=1:length(avals)
    
    awhere=a(:,akey)==avals(i);
 
    bwhere=b(:,bkey)==avals(i);
    
    if sum(bwhere)~=0
        
        bdata=repmat(b(bwhere,get),sum(awhere),1);
    else
        
        blank=nan(1,length(get));

        bdata=repmat(blank,sum(awhere),1);
    end
   
    
    tempmerge(awhere,:)=bdata;
    

    
end

newdata=tempmerge;

merged=[a tempmerge];