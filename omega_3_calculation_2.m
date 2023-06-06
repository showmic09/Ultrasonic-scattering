%%%load voxelID before running the code

dim=500-6*18; %%%dimension of the block 
%cd(directory);

    voxelID=reshape(voxelID,dim*[1 1 1]);
    grainid=unique(voxelID);
count=0;
for ii=1:numel(grainid)
    vox_vol(ii)=numel(find(voxelID==grainid(ii)));
end

%disp(vox_vol);
tic;
%for ii=1:numel(grainid)
for ii=2140:2140 
    count=0;
    id_no=ii;
    x_sum=0;
    y_sum=0;
    z_sum=0;
    mu_2=[0 0 0];
    mu_110=0;
    mu_101=0;
    mu_011=0;
    [x_val,y_val,z_val]=ind2sub(size(voxelID),find(voxelID==grainid(ii)));
    height=unique(z_val);
    %%%%part for calculating aspect ratio
    for jj=1:numel(height)
        [y_plane,x_plane]=find(voxelID(:,:,height(jj))==grainid(ii));
        xlen(jj)=numel(unique(x_plane));
        ylen(jj)=numel(unique(y_plane));
    end
    zlen=numel(height);
    sort_val=[max(xlen) max(ylen) zlen];
    ratio_len{ii}=sort(sort_val);
    abyb(ii)=ratio_len{ii}(3)/ratio_len{ii}(2);
    bbyc(ii)=ratio_len{ii}(2)/ratio_len{ii}(1);
    %%%%%%%
    x_dist=y_val;
    x_dist=x_dist-0.5;
    y_dist=x_val;
    y_dist=y_dist-0.5;
    x_sum=x_sum+ sum(x_dist)
    y_sum=y_sum+ sum(y_dist)
    z_dist=z_val-0.5;
    z_sum=z_sum+sum(z_dist)
    origin=[x_sum/vox_vol(ii) y_sum/vox_vol(ii) z_sum/vox_vol(ii)] %%%location of centroid of the grain
    mu_2=mu_2+[sum((x_dist-origin(1)).^2) sum((y_dist-origin(2)).^2) sum((z_dist-origin(3)).^2)];
    mu_110=mu_110+sum((x_dist-origin(1)).*(y_dist-origin(2)));
    mu_101=mu_101+sum((x_dist-origin(1)).*(z_dist-origin(3)));
    mu_011=mu_011+sum((z_dist-origin(3)).*(y_dist-origin(2)));
    O_3(ii)=mu_2(1)*mu_2(2)*mu_2(3)+2*mu_110*mu_101*mu_011-mu_2(1)*mu_011^2-mu_2(2)*mu_101^2-mu_2(3)*mu_110^2;
    omega_3(ii)=(vox_vol(ii))^5/O_3(ii)
    omega_3_bar(ii)=omega_3(ii)/(2000*pi^2/9);
    clear xlen ylen zlen
end
toc;
save('omega3_sd_15_normal','omega_3_bar'); %%omega_3_bar is the omega3 that is used in DREAM.3D
save('ratio_sd15','ratio_len','abyb','bbyc');