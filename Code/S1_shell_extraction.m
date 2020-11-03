clear;

dataDir='V:\david m\export040920\spherod_segmentation[2104]\calcein_penetration_phnix2_01[4067]\2020-01-24T134730-0500[4123]\Code\Data\';

%all c1
c1DataDir=[dataDir, 'c1\'];

%c1 at working plane
c1WorkingDir=[dataDir, 'c1_working_plane\'];

%spherod segmentation mask of c1 at working plane (by maskRCNN)
c1MaskDir=[dataDir, 'c1_working_plane_mask\'];

%normalized c1 at working plane
c1NormalizeDir=[dataDir, 'c1_working_plane_normalize\'];
if ~exist(c1NormalizeDir, 'dir')
    mkdir(c1NormalizeDir);
end

%shells of c1
c1ShellDir=[dataDir, 'c1_working_plane_shell\'];
if ~exist(c1ShellDir, 'dir')
    mkdir(c1ShellDir);
end

%write tiff options
options.color     = false;
options.compress  = 'no';
options.message   = true;
options.append    = false;
options.overwrite = true;
options.big       = false;

%list of all wells
rawfiles=dir(fullfile([c1DataDir  '*.tif']));
for ii=1:length(rawfiles)
    raw_name=rawfiles(ii).name;
    index1=strfind(raw_name,'_');
    well_all{ii}=raw_name(1:index1(1)-1);
    
end
well_list=unique(well_all);

for mm=1:length(well_list)
    wellname=well_list{mm}
    
    %Step 1: determine bottom/working plane for each well
    bottom_plane=bottom_plane_computation(c1DataDir, wellname);
    working_plane=bottom_plane+5;
    
    %     % copy working wells to a seperate dir
    %     workingFileName=[c1DataDir wellname '_f1_z' num2str(working_plane) '_c1.tif'];
    %     copyfile(workingFileName,c1WorkingDir);
     
    %Step 2: spherod segmentation by maskRCNN (here we assume we already
    %had all spherod segmentations)
    
    %c1 of working plane
    c1_name=[wellname '_f1_z' num2str(working_plane) '_c1.tif'];
    %spheroid segmentation by maskRCNN
    c1_mask_name=strrep(c1_name,'.tif','_nuclei_post.tif');
    
    c1Data = loadtiff([c1DataDir c1_name]);
    [ny, nx, nz] = size(c1Data);
    
    c1MaskData = loadtiff([c1MaskDir c1_mask_name]);
    
    %Step 3: normalization c1
    csvFile=[c1NormalizeDir c1_name(1:end-4) '.csv'];
    c1NormData=normalize_c1(c1Data,c1MaskData,bottom_plane,working_plane,csvFile,wellname);
    saveastiff(uint16(c1NormData),[c1NormalizeDir c1_name(1:end-4) '_norm.tif'],options);
    
    %Step 4: Shell extraction and compute mean FU of each shell
    csvFile=[c1ShellDir c1_name(1:end-4) '.csv'];
    c1ShellData=shellExtraction(uint16(c1NormData),c1MaskData,csvFile);
    saveastiff(uint16(c1ShellData), [c1ShellDir c1_name(1:end-4) '_shells.tif'], options);
    
    disp('pause for next well ...\n');
    pause
    close all;
    
end

function bottom_plane=bottom_plane_computation(c1DataDir,wellname)

rawfiles=dir(fullfile([c1DataDir  [wellname '_f1_*.tif']]));

for ii=1:length(rawfiles)
    raw_name=rawfiles(ii).name;
    index1=strfind(raw_name,'_');
    well=raw_name(1:index1(1)-1);
    slice=str2num(raw_name(index1(2)+2:index1(3)-1));
    sliceTiff=loadtiff([c1DataDir,raw_name]);
    channel1(:,:,slice)=sliceTiff(:,:,30);
end

%plot of mean attenuation to determine bottom plane
mean_att=[];
zplane=[];
for zz=1:size(channel1,3)
    
    sliceI=channel1(:,:,zz);
    mean_att=[mean_att,mean(nonzeros(sliceI(:)))];
    zplane=[zplane,zz];
    
end

plot(mean_att);
title('plot of mean FU attenuation to determine bottom/working plane');

indx=find(mean_att>=max(mean_att)*0.8);
bottom_plane=zplane(indx(1));
% working_plane=bottom_plane+5;

end

function c1_norm=normalize_c1(c1_data,c1_maskData,bottom_plane,working_plane,csvFile,wellname)

FU_all=[];
FU=[];

vox_x=1.19595216191;
vox_z=12;

[ny, nx, nz] = size(c1_data);


%from mean curve to determine cutoff_distance to avoid artifically
%increasing in the center dark region
min_dist2center=[];
for zz=1:nz
    cell1=c1_maskData(:,:,zz);
    cell_dist=bwdist(1-cell1);
    min_dist2center=[min_dist2center max(cell_dist(:))];
end
min_dist2center=min(nonzeros(min_dist2center));

for zz=1:nz
    
    rawI=c1_data(:,:,zz);
    cell1=c1_maskData(:,:,zz);
    cell_dist=bwdist(1-cell1);
    
    FU=[];
    dist2surface=[];
    
    for nn=1:5:round(min_dist2center)
        
        boundaryI=(cell_dist > (nn-1) & cell_dist <=nn);
        boundaryI=rawI.*uint16(boundaryI);
        FU=[FU; mean(nonzeros(boundaryI))];
        dist2surface=[dist2surface;nn];
    end
    
    dist2center=min_dist2center-dist2surface;
    FU_all(zz,:)=FU;
end

FU_mean=mean(FU_all);

%     figure;
%     plot(dist2center,FU_mean,'b');

figure;
FU_diff=abs(diff(FU_mean));
cutoff=max(FU_diff(:))*0.1;
indx=find(FU_diff>=cutoff);
cutoff_indx=indx(end);
cutoff_dist=dist2center(cutoff_indx);


fid3 = fopen(csvFile, 'w');
fprintf(fid3, '%s,%s,%s,%s\n','time point','distance (um from center)','normalized_attenuation','original_attenuation');

c1_norm=zeros(size(c1_data));

for zz=1:nz
    
    rawI=c1_data(:,:,zz);
    cell1=c1_maskData(:,:,zz);
    cell_dist=bwdist(1-cell1);
    
    % %     %y = 1.1051205 * e^(-0.0207123*x) - 0.0183193
    depth=(working_plane-bottom_plane)*vox_z*vox_x*cell_dist./(vox_x*max(cell_dist(:)));
    depth_cutoff=(working_plane-bottom_plane)*vox_z*vox_x*cutoff_dist/(vox_x*max(cell_dist(:)));
    depth(depth>depth_cutoff)=depth_cutoff;
    
    att_decay=exp(log(1.1051205)-0.0207123.*depth)- 0.0183193;
    att_decay(att_decay>1)=1;
    
    norm_rawI=(1-att_decay)*abs(FU_all(zz,end)-FU_all(zz,1))+double(rawI);
    norm_rawI=uint16(norm_rawI);
    c1_norm(:,:,zz)=norm_rawI;
    
    %     norm_rawI=double(rawI)./att_decay;
    
    FU=[];
    dist2surface=[];
    
    FU_n=[];
    
    for nn=0:2:round(max(cell_dist(:)))-2
        
        boundaryI=(cell_dist > (nn) & cell_dist <=nn+2);
        boundaryI_n=norm_rawI.*uint16(boundaryI);
        boundaryI_o=rawI.*uint16(boundaryI);
        
        FU=[FU; mean(nonzeros(boundaryI_o))];
        dist2surface=[dist2surface;nn];
        
        FU_n=[FU_n; mean(nonzeros(boundaryI_n))];
        
    end
    
    dist2center=max(cell_dist(:))-dist2surface;
    plot(vox_x*dist2center,FU,'r'); hold on;
    plot(vox_x*dist2center,FU_n,'b');
    
    for nn=1:length(dist2center)
        fprintf(fid3, '%d,%d,%f,%f\n', zz, round(vox_x*dist2center(nn)), FU_n(nn), FU(nn));
    end
    
    
end

fclose(fid3);

title(wellname);
xlabel('distance to center (um)');
ylabel('FU intensity');

end

function c1_shells=shellExtraction(c1_normalized,c1_mask,csvFile)

[ny, nx, nz] = size(c1_normalized);

fid3 = fopen(csvFile, 'w');
fprintf(fid3, '%s,%s,%s\n','time point','shell','mean_attenuation');

c1_shells=zeros(size(c1_normalized));

for zz=1:nz
    
    rawI=c1_normalized(:,:,zz);
    
    cell1=c1_mask(:,:,zz);
    cell1=logical(cell1);
    
    props1 = regionprops(cell1, 'Image', 'Centroid','Area','BoundingBox');
    largestPart=props1(1).Image;
    %     imtool(largestPart);
    %     imtool close all;
    
    factor=sqrt(0.75);
    tform = affine2d([factor 0 0; 0 factor 0; 0 0 1]);
    ref = imref2d(size(cell1)); %relate intrinsic and world coordinates
    cell75_s = imwarp(largestPart,tform,'OutputView',ref);
    cell75_s = bwareafilt(cell75_s, 1);
    props75 = regionprops(cell75_s, 'Centroid','Area','BoundingBox');
    
    transx=abs(props1(1).Centroid(1)-props75.Centroid(1));
    transy=abs(props1(1).Centroid(2)-props75.Centroid(2));
    
    tform = affine2d([1 0 0; 0 1 0;  transx transy 1]);
    cell75_t = imwarp(cell75_s,tform,'OutputView',ref);
    shell75=(cell1-cell75_t);
    
    shell75_crop=uint16(shell75).*rawI;
    att_shell75=mean(nonzeros(shell75_crop));
    fprintf(fid3, '%d,%d,%f\n', zz, 4, att_shell75);
    
    factor=sqrt(0.5);
    sform = affine2d([factor 0 0; 0 factor 0; 0 0 1]);
    ref = imref2d(size(cell1)); %relate intrinsic and world coordinates
    cell50_s = imwarp(largestPart,sform,'OutputView',ref);
    cell50_s = bwareafilt(cell50_s, 1);
    props50 = regionprops(cell50_s, 'Centroid');
    
    transx=abs(props1(1).Centroid(1)-props50.Centroid(1));
    transy=abs(props1(1).Centroid(2)-props50.Centroid(2));
    
    tform = affine2d([1 0 0; 0 1 0;  transx transy 1]);
    cell50_t = imwarp(cell50_s,tform,'OutputView',ref);
    shell50=(cell75_t-cell50_t);
    
    shell50_crop=uint16(shell50).*rawI;
    att_shell50=mean(nonzeros(shell50_crop));
    fprintf(fid3, '%d,%d,%f\n', zz, 3, att_shell50);
    
    
    factor=sqrt(0.25);
    sform = affine2d([factor 0 0; 0 factor 0; 0 0 1]);
    ref = imref2d(size(cell1)); %relate intrinsic and world coordinates
    cell25_s = imwarp(largestPart,sform,'OutputView',ref);
    cell25_s = bwareafilt(cell25_s, 1);
    props25 = regionprops(cell25_s, 'Centroid');
    
    transx=abs(props1(1).Centroid(1)-props25.Centroid(1));
    transy=abs(props1(1).Centroid(2)-props25.Centroid(2));
    
    tform = affine2d([1 0 0; 0 1 0;  transx transy 1]);
    cell25_t = imwarp(cell25_s,tform,'OutputView',ref);
    shell25=(cell50_t-cell25_t);
    
    shell25_crop=uint16(shell25).*rawI;
    att_shell25=mean(nonzeros(shell25_crop));
    fprintf(fid3, '%d,%d,%f\n', zz, 2, att_shell25);
    
    shell10=double(cell1)-shell75-shell25-shell50;
    
    shell10_crop=uint16(shell10).*rawI;
    att_shell10=mean(nonzeros(shell10_crop));
    fprintf(fid3, '%d,%d,%f\n', zz, 1, att_shell10);
    
    c1_shells(:,:,zz)=10*shell10+75*shell75+50*shell50+25*shell25;
    
end

fclose(fid3);


end

