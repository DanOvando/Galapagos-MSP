function CriticalHabitat= FindHabitat(BottomTypes,maps)

if length(BottomTypes)==1
    CriticalHabitat=maps.fondos==BottomTypes;
else
    temp=[];
    FondoWeights=maps.fondos;
    FondoWeights(isnan(FondoWeights))=0;
    BottomTypes(:,2)=BottomTypes(:,2)./sum(BottomTypes(:,2));
    for x=1:length(BottomTypes)
        
        where=maps.fondos==BottomTypes(x,1);
        FondoWeights(where)=BottomTypes(x,2);
        
        if x<length(BottomTypes)
            temp=strcat(temp,'maps.fondos==',num2str(BottomTypes(x,1)),'|');
        else
            temp=strcat(temp,'maps.fondos==',num2str(BottomTypes(x,1)));
        end
        
    end
    
    CriticalHabitat=eval(temp);
    
    CriticalHabitat=CriticalHabitat.*FondoWeights;
    
    CriticalHabitat(CriticalHabitat==0)=NaN;
    
    CriticalHabitat=CriticalHabitat./nansum(nansum(CriticalHabitat));

end