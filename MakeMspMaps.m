function MakeMspMaps(Layer,LayerName,output,figurefolder,cap)

% Layer=ConflictLayer
% LayerName='ConflictLayer'

Layer(Layer==0)=NaN;
csvwrite(strcat(output,'/',LayerName,'.csv'),Layer)
Layer(Layer>cap*max(max(Layer)))=cap*max(max(Layer));
figure
pcolor(Layer)
shading flat
title(LayerName)
colorbar
print(gcf,'-dtiff',strcat(figurefolder,'/',LayerName,'.tiff'))
close

end