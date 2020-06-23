function NoBoneMask = ...
    RemoveBoneFromThighVolume...
    (OneLegVolume, threshold, strelSizeFactor )
% Remove the cortical bone.
midSlice = round(size(OneLegVolume, 3)/2);
strelSize = round(sqrt(numel(OneLegVolume(:,:,1))) / strelSizeFactor);
se = strel('disk', strelSize);

MuscleandFatVolume = OneLegVolume>threshold;
% Use Otsu thresholding to remove cortical bone.
% MuscleandFatVolume = matitk('SOT', [512] , single(OneLegVolume));
% MuscleandFatVolume = MuscleandFatVolume > 0;

NoBoneMask = zeros(size(MuscleandFatVolume));
NoBoneMask = MuscleandFatVolume;


% Morphological opening followed by 
% boundary tracing to keep the parent object
% thus removing the bone marrow.
% for i=1:size(OneLegVolume, 3)
% %     if i==midSlice
% %         figure, subplot(131), imagesc(OneLegVolume(:,:,i)),
% %         colormap(gray), title('Original'), colorbar, axis image
% %     end
% 
%     % Morphological opening.
%     MuscleandFatSlice = imopen(MuscleandFatVolume(:,:,i), se);
%     
%     % Conn. components.
%     [~, ObjectMapSlice, nObjects] = bwboundaries(MuscleandFatSlice, 'noholes');
%   
%     NoBoneMask(:,:,i) = ObjectMapSlice==1;
% end
% 
NoBoneMask = int16(NoBoneMask);

end
