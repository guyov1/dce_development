function [Rval,  Pval] = CorrNiiMaps( Path_Map1, Path_Map2, SliceNum, MinRange, MaxRange, Mask1, Mask2, Mask3, Mask4, Use_Mask, Mask_Names, Remove_Extreme_Vals)


NiiMap1   = loadniidata(Path_Map1);
NiiMap2   = loadniidata(Path_Map2);

NiiMap1   = NiiMap1(:, :, SliceNum);
NiiMap2   = NiiMap2(:, :, SliceNum);

Mask = uint8(Mask1) | uint8(Mask2) | uint8(Mask3) | uint8(Mask4);

NiiMap1_Masked = zeros(size(NiiMap1));
NiiMap2_Masked = zeros(size(NiiMap2));

if Use_Mask
    NiiMap1_Masked(Mask(:, :, SliceNum)==1)   = NiiMap1(Mask(:, :, SliceNum)==1);
    NiiMap2_Masked(Mask(:, :, SliceNum)==1)   = NiiMap2(Mask(:, :, SliceNum)==1);
    
    figure;
    subplot(2,2,1);
    imshow(NiiMap1);
    title(['DCE Map. Slice: ' num2str(SliceNum)]);
    subplot(2,2,2);
    imshow(NiiMap1_Masked,'Colormap',jet(255));
    subplot(2,2,3);
    imshow(NiiMap2);
    title(['DSC Map. Slice: ' num2str(SliceNum)]);
    subplot(2,2,4);
    imshow(NiiMap2_Masked,'Colormap',jet(255));
else
    NiiMap1_Masked = NiiMap1;
    NiiMap2_Masked = NiiMap2;
end

maxVal1   = max(max(max(NiiMap1_Masked)));
maxVal2   = max(max(max(NiiMap2_Masked)));
minVal1   = min(min(min(NiiMap1_Masked)));
minVal2   = min(min(min(NiiMap2_Masked)));

idx1      = find( NiiMap1_Masked > (MinRange*minVal1) & NiiMap1_Masked < (MaxRange*maxVal1) );
idx2      = find( NiiMap2_Masked > (MinRange*minVal2) & NiiMap2_Masked < (MaxRange*maxVal2) );
final_idx = intersect(idx1, idx2);

msk1_idx  = find(Mask1(:, :, SliceNum)~=0);
msk2_idx  = find(Mask2(:, :, SliceNum)~=0);
msk3_idx  = find(Mask3(:, :, SliceNum)~=0);
msk4_idx  = find(Mask4(:, :, SliceNum)~=0);

if Remove_Extreme_Vals
   not_max_idx_1 = find(NiiMap1_Masked ~= max(max(NiiMap1_Masked)));
   not_max_idx_2 = find(NiiMap2_Masked ~= max(max(NiiMap2_Masked)));
   not_min_idx_1 = find(NiiMap1_Masked ~= min(min(NiiMap1_Masked)));
   not_min_idx_2 = find(NiiMap2_Masked ~= min(min(NiiMap2_Masked)));
   
   new_idx       = intersect(not_max_idx_1, not_max_idx_2);
   new_idx       = intersect(new_idx, not_min_idx_1);
   new_idx       = intersect(new_idx, not_min_idx_2);
   final_idx     = intersect(new_idx,final_idx);
else
   new_idx       = final_idx;
end
   
vector1   = NiiMap1_Masked( final_idx);
vector2   = NiiMap2_Masked( final_idx);

%[R, PValue] = corrplot([vector1  vector2], 'testR','on');
[Rval, Pval] = corrplot([vector1  vector2], 'testR','on');

display(['-I- Displaying ' num2str(length(vector1)) ' voxels.']);

figure;
idx1 = intersect(new_idx,msk1_idx);
idx1 = intersect(idx1,final_idx);
idx2 = intersect(new_idx,msk2_idx);
idx2 = intersect(idx2,final_idx);
idx3 = intersect(new_idx,msk3_idx);
idx3 = intersect(idx3,final_idx);
idx4 = intersect(new_idx,msk4_idx);
idx4 = intersect(idx4,final_idx);

h1 = scatter(NiiMap1_Masked( idx1 ), NiiMap2_Masked( idx1 ), 'r');
hold on;

h2 = scatter(NiiMap1_Masked( idx2 ), NiiMap2_Masked( idx2 ), 'g');
h3 = scatter(NiiMap1_Masked( idx3 ), NiiMap2_Masked( idx3 ), 'b');
h4 = scatter(NiiMap1_Masked( idx4 ), NiiMap2_Masked( idx4 ), 'c');
hold off;
%SliceNum = 4 ; Path1 ='\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\Perfusion_DCE\coregFlow\DSC2DCE\Flow_Larsson_Relative_WM_30_6_brain_Thresholded_200_Normalized_0_1.nii' ; Path2 ='\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\Perfusion_DCE\coregFlow\DSC2DCE\rdsc_oCBFlr_Thresholded_200_Normalized_0_1.nii' ;CorrNiiMaps(Path1, Path2,4);
legend( [h1 h2 h3 h4], char(Mask_Names{1}), char(Mask_Names{2}), char(Mask_Names{3}), char(Mask_Names{4}));


end



