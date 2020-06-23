function OutputVolume = TransposeVolume(InputVolume)
% Used after reading and writing routines.

% OutputVolume = zeros([size(InputVolume,2),size(InputVolume,1),size(InputVolume,3)]);
OutputVolume = zeros(size(InputVolume));

for i=1:size(InputVolume, 3)
    OutputVolume(:,:,i) = InputVolume(:,:,i)';
end
OutputVolume = int16(OutputVolume);
end