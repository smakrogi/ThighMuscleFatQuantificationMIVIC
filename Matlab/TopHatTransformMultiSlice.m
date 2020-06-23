function [EnhancedNonSuppressed, ...
    EnhancedFatSuppressed, ...
    EnhancedWaterSuppressed] = ...
    TopHatTransformMultiSlice...
    (NonSuppressed, ...
    FatSuppressed, ...
    WaterSuppressed, ...
    ExperimentInfo, ...
    strelSizeFactor, ...
    allVolumesFlag)

midSlice = round(size(NonSuppressed, 3)/2);
strelSize = round(sqrt(prod(size(NonSuppressed(:,:,1)))) / strelSizeFactor);
se = strel('disk', strelSize);

EnhancedNonSuppressed = zeros(size(NonSuppressed));
EnhancedWaterSuppressed = zeros(size(WaterSuppressed));
EnhancedFatSuppressed = zeros(size(FatSuppressed));

for i=1:size(NonSuppressed, 3)
    if i==midSlice
        figure, subplot(131), imagesc(NonSuppressed(:,:,i)),
        colormap(gray), title('Original'), colorbar, axis image
    end
    EnhancedNonSuppressed(:,:,i) = imtophat( NonSuppressed(:,:,i), se );
    if i==midSlice
        subplot(132), imagesc(imopen( NonSuppressed(:,:,i), se)),
        colormap(gray), title('After opening'), colorbar, axis image
        subplot(133), imagesc(EnhancedNonSuppressed(:,:,i)),
        colormap(gray), title('After top hat transform'), colorbar, axis image
        saveas(gcf, [ExperimentInfo.SubjectID '_top_hat_transform.png']);
    end
    if( allVolumesFlag )
        EnhancedFatSuppressed(:,:,i) = imtophat( FatSuppressed(:,:,i), se );
        EnhancedWaterSuppressed(:,:,i) = imtophat(WaterSuppressed(:,:,i), se);
    end
end

end

