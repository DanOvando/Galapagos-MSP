%% Marine Spatial Planning for Galapagos
% Dan Ovando
% Constructed to score the performance of marine spatial planning
% alternatives for the Galapagos Marine Reserve
%% Controlfile
% Set model parameters
clear all
close all
global output figurefolder maps

tic
set(0,'DefaultAxesFontName', 'Arial')
datafolder='DATAFOLDER';
propfolder='PROPOSALS';
output='MSP RESULTS Test';
figurefolder='MSP FIGURES Test';
mkdir(output)
mkdir(datafolder)
mkdir(figurefolder)
mkdir(propfolder)

scale=500; %each cell equals XXmeters
propnums=1; %which proposals to evaluate
noisl=1; %Set this to 1 to knock out all islands
setdistance=.08; %set this to the distances from shore that you want to include in the model, set to -1 to include all

mapfolder=['MAPFOLDER ', num2str(scale)];



%% Data
%Use this section to pull in and name data
% 1. Baselayer
% 2. Zoning
% 3. All kinds of data

%Use this section to individually pull in and name each layer. Pain in the
%ass but clearer
maps.blayer=flipud(csvread(strcat(mapfolder,'/basemap.csv'))); %baselayer of islands and oceans
maps.isls=flipud(csvread(strcat(mapfolder,'/islands.csv'))); %baselayer of islands and oceans
maps.bzones=flipud(csvread(strcat(mapfolder,'/rmg_zones_2006.csv'))); %current zoning sceme
maps.tvisit=flipud(csvread(strcat(mapfolder,'/tourism_visits.csv'))); %tourist vists
maps.upwell=flipud(csvread(strcat(mapfolder,'/upwelling.csv'))); %upwelling
maps.sstmin=flipud(csvread(strcat(mapfolder,'/sst_2003_2008_min.csv'))); %min sst
maps.sstmax=flipud(csvread(strcat(mapfolder,'/sst_2003_2008_max.csv'))); %max sst
maps.sstmean=flipud(csvread(strcat(mapfolder,'/sst_2003_2008_avg.csv'))); %mean sst
maps.sstsdv=flipud(csvread(strcat(mapfolder,'/sst_2003_2008_stdev.csv'))); %stdev of sst
maps.lobland=flipud(csvread(strcat(mapfolder,'/lobster_2008_catch_kriging.csv'))); %lobster catch (#s)
maps.lobfish=flipud(csvread(strcat(mapfolder,'/lobster_2008_effort_kriging.csv'))); %lobster effort (diver days)
maps.lobcpue=flipud(csvread(strcat(mapfolder,'/lobster_2008_mdcpue_kriging.csv'))); %lobster CPUE (lobster per diver hours)
% maps.lobcpue=flipud(csvread(strcat(mapfolder,'/lobster2008mediancpue.csv'))); %lobster CPUE (lobster per diver hours)
%maps.lobland=flipud(csvread(strcat(mapfolder,'/lobster2008catch.csv'))); %lobster catch (#s)
%maps.lobfish=flipud(csvread(strcat(mapfolder,'/lobster2008effort.csv'))); %lobster effort (diver days)
maps.pepland=flipud(csvread(strcat(mapfolder,'/pepino_2008_catch_kriging.csv'))); %pepino catch (#s)
maps.pepfish=flipud(csvread(strcat(mapfolder,'/pepino_2008_effort_kriging.csv'))); %pepino diver-hours
maps.pepcpue=flipud(csvread(strcat(mapfolder,'/pepino_2008_mdcpue_kriging.csv'))); %median pepino/diver hour
% maps.pepland=flipud(csvread(strcat(mapfolder,'/pepino2011catch.csv'))); %pepino catch (#s)
% maps.pepfish=flipud(csvread(strcat(mapfolder,'/pepino2008effort.csv'))); %pepino diver-hours
% maps.pepcpue=flipud(csvread(strcat(mapfolder,'/pepino2008mediancpue.csv'))); %median pepino/diver hour

maps.turtles=flipud(csvread(strcat(mapfolder,'/turtle_nests.csv'))); %sea turtle nesting sites
maps.fondos=flipud(csvread(strcat(mapfolder,'/fondos.csv'))); %bottom habitat types
maps.chloro=flipud(csvread(strcat(mapfolder,'/chl_a_avg_krig.csv'))); %Average chlorophyl production
maps.endemic=flipud(csvread(strcat(mapfolder,'/endemic_galapagos_species_richness.csv'))); %Endemic species richness
maps.indo_species=flipud(csvread(strcat(mapfolder,'/indo_pacific_species_richness.csv'))); % Indo pacific species richness
maps.norsubtrop_species=flipud(csvread(strcat(mapfolder,'/northern_sub_tropical_species_richness.csv'))); % northern subtropical species richness
maps.nortemp_species=flipud(csvread(strcat(mapfolder,'/northern_temperate_species_richness.csv'))); % northern subtropical species richness
maps.sousubtrop_species=flipud(csvread(strcat(mapfolder,'/southern_sub_tropical_species_richness.csv'))); % northern subtropical species richness
maps.soutemp_species=flipud(csvread(strcat(mapfolder,'/southern_temperate_species_richness.csv'))); % southern temperate species richness
maps.all_species=flipud(csvread(strcat(mapfolder,'/total_species_richness.csv'))); % all species richness
maps.tropic_species=flipud(csvread(strcat(mapfolder,'/tropical_species_richness.csv'))); % tropical species richness
maps.distancetoshore=flipud(csvread(strcat(mapfolder,'/distance_land.csv'))); %Distance to shore
maps.KrigedTourism=flipud(csvread(strcat(mapfolder,'/tourism_visits_kriging.csv'))); %Kriged tourism data
maps.MacroZones=flipud(csvread(strcat(mapfolder,'/macrozones.csv'))); %Kriged tourism data
data=fieldnames(maps); %pull out names of each data layer


IslandLocations=maps.isls;
IslandLocations(IslandLocations==-9999)=NaN;
IslandLocations(isnan(IslandLocations)==0)=0;

mstack=NaN(size(maps.blayer,1),size(maps.blayer,2),length(data)); %create empty array to stack maps for easy manilulation
%Fill in mstack witha approproiate data
for d=1:length(data)
    eval(strcat(['maps.',data{d},'(maps.',data{d},'==-9999)=NaN;']))
    eval(strcat(['mstack(:,:,',num2str(d),')=','maps.',data{d},';']))
    
    if noisl==1
        temp=mstack(:,:,d);
        temp(maps.blayer==1)=NaN;
        mstack(:,:,d)=temp;
    end
    
    
    if setdistance~=-1
        temp=mstack(:,:,d);
        temp(maps.distancetoshore>=setdistance)=NaN;
        mstack(:,:,d)=temp;
    end
    
    figure
    pcolor(mstack(:,:,d))
    shading flat
    title(data{d})
    print(gcf,'-dtiff',strcat(figurefolder,'/',data{d},'.tiff'))
    close
end

%% Load in zoning proposals and rules

%Once you get some actual proposals here, automate a system to pull in each
%proposal as a part of a layer by number
PropNames={'A','B','Base'};
Props=NaN(size(maps.blayer,1),size(maps.blayer,2),length(PropNames)); %create empty array to stack maps for easy manilulation
for p=1:length(PropNames)
    
    PropTemp=flipud(csvread(strcat(propfolder,'/Prop',PropNames{p},'.csv')));
    PropTemp(PropTemp==-9999)=NaN;
    PropTemp(IslandLocations==0)=IslandLocations(IslandLocations==0);
    
    Props(:,:,p)=PropTemp;
    
    figure
    pcolor(Props(:,:,p))
    shading flat
    title(['Proposal ',PropNames{p}])
    print(gcf,'-dtiff',strcat(figurefolder,'/Prop',PropNames{p},'.tiff'))
    close  
end

props=importdata(strcat(propfolder,'/PROPOSALS.csv'),',',1);
weights=importdata(strcat(propfolder,'/WEIGHTS.csv'),',',1);
weights.data(:,3:end)=weights.data(:,3:end)./repmat(sum(weights.data(:,3:end),2),1,size(weights.data(:,3:end),2)); %Scale the weights

RAWmstack=mstack;


%% Predict habitat quality 
[LobsterPredictedHabitat, LobsterHabitatModel]=PredictHabitat(maps.lobcpue,'lobster',RAWmstack,{'fondos','sstmean','upwell','chloro'},[0,1,0,1],data,[1,2,15,18,31,34,54,60]);
mstack(:,:,end+1)=LobsterPredictedHabitat;
data{end+1}='LobsterPredictedHabitat';

LobsterHabitatModel.r

[PepinoPredictedHabitat, PepinoHabitatModel]=PredictHabitat(maps.pepcpue,'pepino',RAWmstack,{'fondos','sstmean','upwell','chloro'},[0,1,0,1],data,[1,5,18,19,31,54,60]);
mstack(:,:,end+1)=PepinoPredictedHabitat;
data{end+1}='PepinoPredictedHabitat';


%% Process cardinal data
% Use this section to process all that that have a clear numerical ranking

cardinal={'lobland','lobfish','lobcpue','chloro','pepland','pepfish','pepcpue','tvisit','endemic','all_species','KrigedTourism','LobsterPredictedHabitat','PepinoPredictedHabitat'};

for c=1:length(cardinal)
    where=strcmp(data,cardinal{c});
    mstack(:,:,where)=mstack(:,:,where)/nansum(nansum(mstack(:,:,where)));
end


%% Process non-use data

%Scale SST
%How: bin into cold, temperate, hot. Create new layers for these

where=strcmp(data,'sstmean');
temptemp=mstack(:,:,where);
temptemp(maps.blayer==1)=NaN;
ssts=unique(temptemp);
[counts,mids]=hist(ssts,3);
distance=(mids(2)-mids(1))/2;
breaks(1)=mids(1)+distance;
breaks(2)=mids(2)+distance;

cold=temptemp;
cold(temptemp<breaks(1))=1;
cold(temptemp>=breaks(1))=NaN;
where=isnan(cold);
cold(where==0)=1./sum(sum(where==0));

temperate=temptemp;
temperate(temptemp>=breaks(1) & temptemp<breaks(2))=2;
temperate(temptemp<breaks(1) | temptemp>=breaks(2))=NaN;
where=isnan(temperate);
temperate(where==0)=1./sum(sum(where==0));

hot=temptemp;
hot(temptemp>=breaks(2))=3;
hot(temptemp<breaks(2))=NaN;
where=isnan(hot);
hot(where==0)=1./sum(sum(where==0));

data(end+(1:3))={'cold','temperate','hot'};
mstack(:,:,end+1)=cold;
mstack(:,:,end+1)=temperate;
mstack(:,:,end+1)=hot;

%Scale upwelling
where=strcmp(data,'upwell');
uptemp=mstack(:,:,where);
where2=uptemp==1;
uptemp(where2)=1/sum(sum(where2));
uptemp(where2==0)=NaN;
mstack(:,:,where)=uptemp;

 %% Add in Lobster Habitat %%
 
 BottomTypes=[1,2,13,15,18,31,60;1,1.5,2,1.5,1.5,1.5,1]';
 
 CriticalHabitat=FindHabitat(BottomTypes,maps);
 
 MakeMspMaps(CriticalHabitat,'Lobster Habitat',output,figurefolder,1)
 
 mstack(:,:,end+1)=CriticalHabitat;
 
 data{end+1}='lobhabitat';
 
 

 %% Add in Pepino Habitat %%
 
 BottomTypes=[1,2,5,18,19,31;2,1,1.5,1,1.5,1]';
 
 CriticalHabitat=FindHabitat(BottomTypes,maps);
 
 MakeMspMaps(CriticalHabitat,'Pepino Habitat',output,figurefolder,1)
 
 mstack(:,:,end+1)=CriticalHabitat;
 
 data{end+1}='pephabitat';
 
%% Put data into matrix for algorithm

propdata=NaN(size(props.data,1),size(props.data,2)); %make blank matrix to store algorithm data
PropNums=unique(props.data(:,1));
numprops=length(PropNums); %find number of proposals being evaluated
vars=deblank(props.colheaders); %find the variables being used in this analysis
% propmaps.testzone=maps.bzones;
% propmaps.testzone(propmaps.testzone==3)=5; %create a fake 'zone 5' in the curren mixed use areas
c=0;

for p=1:numprops % loop over each proposal
    TempMap=Props(:,:,p);
    wherep=props.data(:,1)==PropNums(p);
    ptemp=props.data(wherep,:); %pull out proposal p
    tzones=ptemp(:,2); %find zones inside that proposal
    for z=1:length(tzones) %loop over zones in proposal p
        wherez=TempMap==tzones(z); %find locations of that zone in proposal p
        c=c+1;
        propdata(c,[1,2])=[PropNums(p),tzones(z)];
        for v=3:length(vars) %loop over variables you want to measure
            whereget=strcmp(vars{v},data); %find the right data layer
            tempdat=mstack(:,:,whereget); %store that data
            tempscore=nansum(nansum(tempdat(wherez))); %find the score of that zone
            whereput=strcmp(vars{v},vars);
            propdata(c,whereput)=tempscore;
        end %close loop over variables
    end %close loop over zones
end %close loop over proposals

 
%% Score and Save Zoning

% [r,c,m]=shpath(maps.blayer,200,180,260,359);

propscore=propdata(:,1:2); %Place to store final results
rules=propdata(:,3:end).*props.data(:,3:end);%Apply rules
scores=100.*(rules.*weights.data(:,3:end)); %Apply weights to rules

%Break into categories
fishing_data={'lobland','lobfish','lobcpue','pepland','pepfish','pepcpue'};
pairs=findmatches(fishing_data,vars(3:end),'character');
fishing=sum(scores(:,pairs),2);
TempWeights=NaN(size(mstack(:,:,1),1),size(mstack(:,:,1),2),length(fishing_data));
ReScaleWeights=sum(weights.data(1,findmatches(fishing_data,vars,'character')));
DataPairs=findmatches(fishing_data,data,'character');
OrderdData= data(DataPairs);
for t=1:length(OrderdData)
    WhereWeight=strcmp(OrderdData(t),vars);
    TempWeight=weights.data(1,WhereWeight)./ReScaleWeights;
    TempWeights(:,:,t)=TempWeight;
end
TempScore=mstack(:,:,DataPairs).*TempWeights;
FishingLayer=nansum(TempScore,3);
FishingLayer(FishingLayer==0)=NaN;
csvwrite(strcat(output,'/FishingLayer.csv'),FishingLayer)
figure
pcolor(FishingLayer)
shading flat
title('Fishing')
colorbar
print(gcf,'-dtiff',strcat(figurefolder,'/FishingLayer.tiff'))
close

tourism_data={'endemic','KrigedTourism'};
pairs=findmatches(tourism_data,vars(3:end),'character');
tourism=sum(scores(:,pairs),2);
TempWeights=NaN(size(mstack(:,:,1),1),size(mstack(:,:,1),2),length(tourism_data));
ReScaleWeights=sum(weights.data(1,findmatches(tourism_data,vars,'character')));
DataPairs=findmatches(tourism_data,data,'character');
OrderdData= data(DataPairs);
for t=1:length(OrderdData)
    WhereWeight=strcmp(OrderdData(t),vars);
    TempWeight=weights.data(1,WhereWeight)./ReScaleWeights;
    TempWeights(:,:,t)=TempWeight;
end

TempScore=mstack(:,:,DataPairs).*TempWeights;
TourismLayer=nansum(TempScore,3);
TourismLayer(TourismLayer==0)=NaN;
csvwrite(strcat(output,'/TourismLayer.csv'),TourismLayer)
figure
pcolor(TourismLayer)
shading flat
title('Tourism')
colorbar
print(gcf,'-dtiff',strcat(figurefolder,'/TourismLayer.tiff'))
close


ecosystem_data={'endemic','all_species','chloro','hot','cold','temperate'};
pairs=findmatches(ecosystem_data,vars(3:end),'character');
ecosystem=sum(scores(:,pairs),2);
TempWeights=NaN(size(mstack(:,:,1),1),size(mstack(:,:,1),2),length(ecosystem_data));
ReScaleWeights=sum(weights.data(1,findmatches(ecosystem_data,vars,'character')));
DataPairs=findmatches(ecosystem_data,data,'character');
OrderdData= data(DataPairs);
for t=1:length(OrderdData)
    WhereWeight=strcmp(OrderdData(t),vars);
    TempWeight=weights.data(1,WhereWeight)./ReScaleWeights;
    TempWeights(:,:,t)=TempWeight;
end
DataPairs=findmatches(ecosystem_data,data,'character');
TempScore=mstack(:,:,DataPairs).*TempWeights; %There's a problem here, the weights and the data aren't in the right order!sky
EcosystemLayer=nansum(TempScore,3);
EcosystemLayer(EcosystemLayer==0)=NaN;
% EcosystemLayer(EcosystemLayer>.5*max(max(EcosystemLayer)))=.5*max(max(EcosystemLayer));
csvwrite(strcat(output,'/EcosystemLayer.csv'),EcosystemLayer)
figure
pcolor(EcosystemLayer)
shading flat
title('Ecosystem')
colorbar
print(gcf,'-dtiff',strcat(figurefolder,'/EcosystemLayer.tiff'))
close

%% Perform Conflict and Low Fruit Analysis
cap=0.5;
a=isnan(EcosystemLayer)==0;
b=isnan(FishingLayer)==0;
c=isnan(TourismLayer)==0;
UsesLayer= a+b+c;
MakeMspMaps(UsesLayer,'UsesLayer',output,figurefolder,cap)
EcoVsFishing=EcosystemLayer.*FishingLayer;
MakeMspMaps(EcoVsFishing,'EcoVsFishing',output,figurefolder,cap)
FishVsTour=FishingLayer.*TourismLayer;
MakeMspMaps(FishVsTour,'FishVsTour',output,figurefolder,cap)
EcoVsTour=EcosystemLayer.*TourismLayer;
MakeMspMaps(EcoVsTour,'EcoVsTour',output,figurefolder,cap)
ConflictLayer= UsesLayer.*(EcosystemLayer.*FishingLayer+FishingLayer.*TourismLayer+EcosystemLayer.*TourismLayer);
MakeMspMaps(ConflictLayer,'ConflictLayer',output,figurefolder,cap)

% Produce Low Hanging Fruit 

EcoWiZeros=EcosystemLayer;
EcoWiZeros(isnan(EcosystemLayer) & maps.distancetoshore<setdistance)=1e-5;
FishWiZeros=FishingLayer;
FishWiZeros(isnan(FishingLayer)& maps.distancetoshore<setdistance)=1e-5;
TourWiZeros=TourismLayer;
TourWiZeros(isnan(TourismLayer)& maps.distancetoshore<setdistance)=1e-5;


MakeMspMaps(EcoWiZeros,'EcoWiZeros',output,figurefolder,cap)
MakeMspMaps(FishWiZeros,'FishWiZeros',output,figurefolder,cap)
MakeMspMaps(TourWiZeros,'TourWiZeros',output,figurefolder,cap)

EasyEco=EcosystemLayer./(FishWiZeros+TourWiZeros);
MakeMspMaps(EasyEco,'EasyEco',output,figurefolder,cap)
EasyFish=FishingLayer./(EcoWiZeros+TourWiZeros);
MakeMspMaps(EasyFish,'EasyFish',output,figurefolder,cap)
EasyTour=TourismLayer./(FishWiZeros+EcoWiZeros);
MakeMspMaps(EasyTour,'EasyTour',output,figurefolder,cap)%Write layer and csv that has the final score for each proposal and subzone

%% Score individual zones %%

%Analyze reductions in effort and habitat protection

WherePepinoFishing=findmatches({'pepfish','PepinoPredictedHabitat'},vars(3:end),'character');
PepinoFishingChange=[propdata(:,1:2) rules(:,WherePepinoFishing)];

WhereLobsterFishing=findmatches({'lobfish','LobsterPredictedHabitat'},vars(3:end),'character');
LobsterFishingChange=[propdata(:,1:2) rules(:,WhereLobsterFishing)];


% propscore=[propscore scores fishing tourism ecosystem sum(scores,2)]; %compile data

propscore=[propscore scores fishing tourism ecosystem fishing+tourism+ecosystem]; %compile data

for p=1:length(PropNums)
   where=propscore(:,1)==PropNums(p);
   IndPropScores(p,:)=[PropNums(p),1,sum(propscore(where,3:end),1)];
end

ScoreNames=[props.colheaders,'Fishing','Tourism','Ecosystem','Total'];

figure
plot3(IndPropScores(:,strcmp('Fishing',ScoreNames)),IndPropScores(:,strcmp('Tourism',ScoreNames)),IndPropScores(:,strcmp('Ecosystem',ScoreNames)),'o','markersize',12)
xlabel('Fishing')
ylabel('Tourism')
zlabel('Ecosystem')
grid on
print(gcf,'-dtiff',strcat(figurefolder,'/3dtradeoff.tiff'))
close

figure
plot(IndPropScores(:,strcmp('Fishing',ScoreNames)),IndPropScores(:,strcmp('Ecosystem',ScoreNames)),'o','markersize',12)
xlabel('Fishing')
ylabel('Ecosystem')
print(gcf,'-dtiff',strcat(figurefolder,'/2dtradeoff.tiff'))
close


figure
barh([IndPropScores(:,strcmp('Fishing',ScoreNames)), IndPropScores(:,strcmp('Tourism',ScoreNames)), IndPropScores(:,strcmp('Ecosystem',ScoreNames)),IndPropScores(:,strcmp('Total',ScoreNames))]')
set(gca,'YTickLabel',[])
set(gca,'YTickLabel',{'Fishing','Tourism','Ecosystem','Total'})
colormap summer
legend(strcat('Proposal',PropNames),'Location','NorthWest')
print(gcf,'-dtiff',strcat(figurefolder,'/bargraph.tiff'))

csvwrite_with_headers(strcat(output,'/Proposal Scores.csv'),propscore,[props.colheaders,'Fishing','Tourism','Ecosystem','Total']) %Spit out final results as a named .csv

%Plot scores of each zone in each proposal

PlotData={'Fishing','Tourism','Ecosystem','Total','endemic'};
for d=1:length(PlotData) %loop over data to be plotted/stored
    
    for p=1:length(PropNums) %loop over proposals
        where=propscore(:,1)==PropNums(p);
        TempProp=propscore(where,:); %Pull out proposal p
        TempZones=unique(TempProp(:,2));
        TempMap=Props(:,:,p);
        for z=1:length(TempZones)
            WhereZone1=TempMap==TempZones(z); %find location of zone z in proposal p
            WhereZone2=TempProp(:,2)==TempZones(z); %find data for zone z in proposal p
            WhereScore=strcmp(PlotData(d),ScoreNames); %find the right data
            TempMap(WhereZone1)=TempProp(WhereZone2,WhereScore);  %assign zone z in proposal p the score for data d
        end
        PlotName=strcat('Proposal',num2str(PropNums(p)),PlotData{d});
        csvwrite(strcat(output,'/',PlotName,'.csv'),TempMap)
        
        figure
        pcolor(TempMap)
        shading flat
        title(PlotName)
        colorbar
        print(gcf,'-dtiff',strcat(figurefolder,'/',PlotName,'.tiff'))
        close
        
    end
    
    
end






