function [PredictedHabitat,Model]=PredictHabitat(cpue,species,mstack,datalayers,continuous,data,BottomTypes)
global maps output figurefolder
% species='lobster';
% BottomTypes=[1,2,15,18,31,34,54,60];
% vars=RAWmstack;
% cpue=maps.lobcpue;
% names={'fondos','sstmean','upwell','chloro'};
% continuous=[0,1,0,1];
Dimension=size(cpue,1).*size(cpue,2);
RegStack=reshape(cpue,Dimension,1);

for i=1:length(datalayers)
    
    where=strcmp(data,datalayers{i});
    
    temp=reshape(mstack(:,:,where),Dimension,1);

    if strcmp(datalayers{i},'fondos')
        WrongBottom= findmatches(BottomTypes,temp,'number')==0;
        HasData=isnan(temp)==0;
        where=WrongBottom ==1 & HasData'==1; 
        temp(where)=0;
    end

    if continuous(i)==0
        temp=dummyvar(nominal(temp));
        temp=temp(:,2:end);
    end
    RegStack=[RegStack temp];
    
end

% cpuetest=isnan(maps.pepcpue);
% fondotest=isnan(maps.fondos);
% huh=cpuetest==0 & fondotest==0;
% sum(sum(huh))

RegStack=[RegStack ones(Dimension,1)];
RegStack2=RegStack(isnan(RegStack(:,1))==0,:);
[b,bint,r,rint,stats]=regress(RegStack2(:,1),RegStack2(:,2:end));
% s=regstats(RegStack2(:,1),RegStack2(:,2:(end-1)),'linear');
Model.b=b;
Model.bint=bint;
Model.r=r;
Model.rint=rint;
Model.stats=stats;
PredictedHabitat=reshape(RegStack(:,2:end)*b,size(maps.fondos,1),size(maps.fondos,2));

MakeMspMaps(PredictedHabitat,strcat(species,'PredictedHabitat'),output,figurefolder,1)

end