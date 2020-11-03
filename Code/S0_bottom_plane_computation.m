clear;

dataDir='V:\david m\export040920\spherod_segmentation[2104]\calcein_penetration_phnix2_01[4067]\2020-01-24T134730-0500[4123]\Code\Data\';

%all c1
c1DataDir=[dataDir, 'c1\'];

%c1 at working plane
c1WorkingDir=[dataDir, 'c1_working_plane\'];
if ~exist(c1WorkingDir, 'dir')
    mkdir(c1WorkingDir);
end

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
    
    % copy working wells to a seperate dir
    workingFileName=[c1DataDir wellname '_f1_z' num2str(working_plane) '_c1.tif'];
    copyfile(workingFileName,c1WorkingDir);
    
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
title(wellname);

indx=find(mean_att>=max(mean_att)*0.8);
bottom_plane=zplane(indx(1));

end

