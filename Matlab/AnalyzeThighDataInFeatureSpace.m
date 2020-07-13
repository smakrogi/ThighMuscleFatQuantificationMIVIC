function [tissueIdx, OutDataMatrix, DistributionModels] = ...
    AnalyzeThighDataInFeatureSpace( subjectNumber, dataMatrix, ExperimentInfo )
% syntax: [tissueIdx, OutDataMatrix] = ...
% AnalyzeThighDataInFeatureSpace( subjectNumber, dataMatrix, ExperimentInfo )
% Feature analysis, clustering and plotting of the input 
% samples.

% Retrieve path to stpr toolbox.
% STPR_PATH = getenv('STPR_PATH');
STPR_PATH = 'C:\Users\smakrogiannis\Documents\Codes\src\thalesrepos\m-files\stprtool';
% DistanceLearningMethods = {'', 'PCA', 'ISOMAP', 'Laplacian', 'LLE', 'LTSA', 'MDS'};
% distancemetriclearningMethod = 'PCA';
distancemetriclearningMethod = ExperimentInfo.DistanceLearningMethod;
nNeighbors = 8;
outputDimensionality = ExperimentInfo.OutputDimensionality;
DistributionModels = [];

% Feature extraction and dimensionality reduction.
switch lower(distancemetriclearningMethod)
    case{''}
    case{'pca'}
        % Run PCA.
%         % Matlab's stats toolbox implementation.
%         [eigVector, scoreMatrix, eigValues] = princomp(dataMatrix);
%         dataMatrix = scoreMatrix;
        % UIUC's implementation.
        optionsPCA.ReducedDim = outputDimensionality;
        [eigVector, eigValues] = PCA(dataMatrix, optionsPCA);
        scoreMatrix = dataMatrix * eigVector;
        dataMatrix = scoreMatrix;
    case{'isomap'} % very slow and memory demanding.
        options.dims = outputDimensionality;
%         DM = L2_distance(dataMatrix',dataMatrix',1);
        DM = L2_distance(dataMatrix',dataMatrix');
        [Y, R] = Isomap(DM, 'k', nNeighbors, options);
        dataMatrix = Y.coords{1}';
    case{'laplacian'} % not so good results on our data.
        [E,V] = leigs(dataMatrix, 'nn', nNeighbors, outputDimensionality+1);
        dataMatrix = E(:,1:outputDimensionality);
    case{'lle'} % very slow, unexpected results.
        dataMatrix = lle(dataMatrix, outputDimensionality, nNeighbors);
    case{'ltsa'}  % slow, unexpected results.
        dataMatrix = ltsa(dataMatrix, outputDimensionality, nNeighbors);
    case('mds') % very slow.
         DM = L2_distance(dataMatrix',dataMatrix');
         dataMatrix = mdsFast(DM, outputDimensionality);
%          dataMatrix = mds(dataMatrix, outputDimensionality);
     case('kpca') % very slow.
%          stprpath(STPR_PATH);
%          model = kpca( dataMatrix', struct('ker','rbf','arg',4,'new_dim',2) );
%          XR = kpcarec( dataMatrix', model );
%          dataMatrix = XR';
         dataMatrix = compute_mapping(dataMatrix, 'KernelPCA', outputDimensionality);
    case('ppca')
        dataMatrix = compute_mapping(dataMatrix, 'ProbPCA', outputDimensionality);
end

infoString = sprintf('Distance metric learning using %s, done. \n', distancemetriclearningMethod);
fprintf( infoString );

% Clustering.
[temp, D] = size(dataMatrix);
clear temp;
selectedClusteringMethod = ExperimentInfo.ClusteringMethod;

switch lower(selectedClusteringMethod)
    case{'k-means'}
        % Apply k-means.
%         seeds = [min(dataMatrix); mean(dataMatrix); max(dataMatrix)];
%         opts = statset('Display','iter');
%         % Use matlab's statistics toolbox.
%         %tissueIdx = kmeans( dataMatrix(:,:), ExperimentInfo.nClusters, 'Start', seeds(1:ExperimentInfo.nClusters,:), 'Options', opts );
%         tissueIdx = kmeans( dataMatrix(:,:), ExperimentInfo.nClusters, 'Options', opts );
        % Use Fuzzy Clustering Toolbox's k-means implementation.
        data.X = dataMatrix;
        param.c = ExperimentInfo.nClusters;
        param.vis = 0;
        result = Kmeans(data,param);
        [temp, tissueIdx] = max(result.data.f');
        clear temp;
        tissueIdx = tissueIdx';
    case{'fcm'}
        data.X = dataMatrix;
        param.c = ExperimentInfo.nClusters;
        result = FCMclust(data,param);
        [temp, tissueIdx] = max(result.data.f');
        clear temp;
        tissueIdx = tissueIdx';
    case('emgmmbaycls') %emgmm method with bayesian decision.
        stprpath(STPR_PATH);
        
        % Optional subtractive clustering for cluster validity.
        % subclustRadii = 0.25;
        % clusterCenters = subclust(dataMatrix, subclustRadii);
        % model = mlcgmm(dataMatrix', struct('ncomp', ExperimentInfo.nClusters, ...
        % 'verb',1, 'init', 'cmeans', 'cov_type', 'full'));
        
        % Build GMM model using EM.
        covarianceMatrixType = lower( ExperimentInfo.CovarianceMatrixType );
        model = emgmm(dataMatrix', struct('ncomp', ...
            ExperimentInfo.GMMnDistributions, 'verb', 0, 'init', ...
            'cmeans', 'cov_type', covarianceMatrixType));
        
        % Generate 2 or 3-class GMM model.
        if (ExperimentInfo.nClusters == 2)
            DistributionModels = Generate2ClassGMM(model, ...
                ExperimentInfo.GMMnDistributions);
        elseif  (ExperimentInfo.nClusters == 3)
                DistributionModels = Generate3ClassGMM(model, ...
                    ExperimentInfo.GMMnDistributions);
        else error('Sorry, unexpected number of clusters.');
        end
        
        % Classification using Bayesian inference.
        tissueIdx = bayescls(dataMatrix', DistributionModels);
        tissueIdx = tissueIdx';
    case('ncutclus') % memory intensive, needs undersampling.
        [W,temp] = compute_relation(single(dataMatrix'));
        clear temp;
        [NcutDiscrete, temp1, temp2] = ...
            ncutW(W, ExperimentInfo.nClusters);
        clear temp1; clear temp2;
        [temp, tissueIdx] = max(NcutDiscrete');
        clear temp;
        tissueIdx = tissueIdx';
        clear W;
    case{'meanshift'}
        bandwidth = ExperimentInfo.MeanShiftBandwidth;
        [temp, tissueIdx] = MeanShiftCluster(dataMatrix',bandwidth);
        clear temp;
        tissueIdx = tissueIdx';
        nDetectedClusters = length(unique(tissueIdx));
        infoString = sprintf('Mean shift detected %d clusters.\n', nDetectedClusters);
        fprintf( infoString );
        if( nDetectedClusters~=2 ) error('Sorry, unexpected number of clusters.'); end
    case('emgmmquad') %emgmm method with bayesian decision boundary and quadratic classifier
        %!! Use only when nClusters=2.
        if( ExperimentInfo.nClusters ~= 2 ) error('Sorry, valid selection for 2 clusters only.'); end
        stprpath(STPR_PATH);
        
        % Build GMM model using EM.
        covarianceMatrixType = lower( ExperimentInfo.CovarianceMatrixType );
        %model = mlcgmm(dataMatrix', struct('ncomp', ExperimentInfo.nClusters, 'verb',1, 'init', ...
        %'cmeans', 'cov_type', 'full'));
        model = emgmm(dataMatrix', struct('ncomp', ExperimentInfo.nClusters, 'verb', 0, ...
            'init', 'cmeans', 'cov_type', covarianceMatrixType));  % cov_type: full, diag, spherical.
        
        % Create a tissue-specific model to be used in quantification.
        DistributionModels = Generate2ClassGMM(model, ...
            ExperimentInfo.GMMnDistributions);
        
        % Create quadratic model.
        quad_model = bayesdf(model);
        tissueIdx = quadclass(dataMatrix',quad_model);
        tissueIdx = tissueIdx';
    case('mst')
%         [patterns, targets] = min_spanning_tree(dataMatrix', ...
%             ones(size(dataMatrix,1)), '[''NN'', 2]');
    case('parzen')
    case('hierarchical')
end

infoString = sprintf('Clustering using %s, done.\n', selectedClusteringMethod);
fprintf( infoString );

% Plot results in 2D and 3D, before and after clustering.
if (ExperimentInfo.commandlineMode)
    purple = [119/255 73/255 152/255];
    plot_labels = {'k', 'r', 'g', purple, 'm', 'y'};
    imagefileprefix = [ExperimentInfo.SubjectID{subjectNumber}, '_', ...
        ExperimentInfo.DistanceLearningMethod, ...
        '_', ExperimentInfo.ClusteringMethod];
    
    plotMarkerSize = 4;
    plotMarkerShape = 'o';
    
    if(D==2)
        % 2D.
        figure,
        subplot(221),
        plot(dataMatrix(:,1), dataMatrix(:,2), plotMarkerShape, 'MarkerSize', plotMarkerSize);
        hold on, xlabel('fat suppressed'); ylabel('water suppressed'); grid on, axis square;
        xlim = get(gca, 'xlim');
        ylim = get(gca, 'ylim');
        set(gca, 'xlim', xlim);
        set(gca, 'ylim', ylim);
%         figure,
        for count=1:ExperimentInfo.nClusters
            vectorIndex = find(tissueIdx==count);
%             subplot(222),
            plot(dataMatrix(vectorIndex,1), dataMatrix(vectorIndex,2), plotMarkerShape, ...
                'MarkerEdgeColor', plot_labels{count}, 'MarkerFaceColor', plot_labels{count}, ...
                'MarkerSize', plotMarkerSize);
            hold on, xlabel('fat suppressed'); ylabel('water suppressed'); grid on, axis square;
            set(gca, 'xlim', xlim);
            set(gca, 'ylim', ylim);
        end
        legend('muscle','interMF','intraMF','group4','group5', 'Location', 'NorthEast');
        switch lower(selectedClusteringMethod)
            case('emgmmbaycls')
%                 figure,
                subplot(223),
                axis square, grid on
                options.fill = 1;
                options.line_style = 'gs';
                pgauss(model, options);
                xlabel('fat suppressed'); ylabel('water suppressed');
                set(gca, 'xlim', xlim);
                set(gca, 'ylim', ylim);
%                 figure,
                subplot(224)
                options.fill = 1;
                options.line_style = 'gs';
                muscleIdx = 1;
                vectorIndex = find(tissueIdx==muscleIdx);
                plot(dataMatrix(vectorIndex,1), dataMatrix(vectorIndex,2), plotMarkerShape, ...
                'MarkerEdgeColor', plot_labels{muscleIdx}, 'MarkerFaceColor', plot_labels{muscleIdx}, ...
                'MarkerSize', plotMarkerSize);
                hold on,
                pgauss(DistributionModels.Pclass{1}, options);
                xlabel('fat suppressed'); ylabel('water suppressed');
                axis square, grid on
            case('emgmmquad')
%                 figure,
                subplot(223),
                axis square, grid on
                options.fill = 1;
                pgauss(model, options)
                options.fill = 0;
                options.line_style = 'rs';
                pboundary(quad_model, options),
                xlabel('fat suppressed'); ylabel('water suppressed');
                set(gca, 'xlim', xlim);
                set(gca, 'ylim', ylim);
        end
        saveas(gcf, [imagefileprefix, '_', 'fat_water_2D_plots.fig']);
        saveas(gcf, [imagefileprefix, '_', 'fat_water_2D_plots.png']);
    end
    
    if( D==3 )
        % 3D.
        azimuth = -37.5; elevation = 30;
        figure, subplot(121), title('3D plots of unlabeled data'), hold on,
        plot3(dataMatrix(:,1),dataMatrix(:,2),dataMatrix(:,3), plotMarkerShape); view(azimuth, elevation);
        xlabel('t1w'); ylabel('fat suppressed'); zlabel('water suppressed', ...
            'MarkerSize', plotMarkerSize);
        grid on, axis square;
        
        subplot(122), title('3D plots of clustered data'), hold on,
        for count=1:ExperimentInfo.nClusters
            vectorIndex = find(tissueIdx==count);
            plot3(dataMatrix(vectorIndex,1),dataMatrix(vectorIndex,2),...
                dataMatrix(vectorIndex,3), plotMarkerShape, 'MarkerEdgeColor', plot_labels{count}, ...
                'MarkerFaceColor', plot_labels{count}, 'MarkerSize', plotMarkerSize);
            hold on;
        end
        view(azimuth, elevation);
        
        %     legend('tissue1','tissue2','tissue3','tissue4','tissue5');
        xlabel('t1w'); ylabel('fat suppressed'); zlabel('water suppressed');
        grid on, axis square;
        saveas(gcf, [imagefileprefix, '_', 't1w_fat_water_3D_plots.fig']);
        saveas(gcf, [imagefileprefix, '_', 't1w_fat_water_3D_plots.png']);        
    end
end

% Return data matrix after dimensionality reduction.
OutDataMatrix = dataMatrix;

% Group samples by cluster number for plottng.
DistributionModels.Samples = CopySamplesByCluster(dataMatrix, tissueIdx, ExperimentInfo.nClusters);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bayesian_model = Generate2ClassGMM(model, nClusters)

classWSMeans = cat(1, model.Mean(2,:));
[temp,fatLabel] = max(classWSMeans);
clear temp;
count=1;
bayesian_model.Pclass = cell(1,2);
bayesian_model.Prior = zeros(1,2);
% Construct a 2-class model from the Gaussians.
for(i=1:nClusters)
    if(i==fatLabel) % interMF group.
        bayesian_model.Pclass{2}.Mean = model.Mean(:,i);
        bayesian_model.Pclass{2}.Cov = model.Cov(:,:,i);
        bayesian_model.Pclass{2}.fun = model.fun;
        bayesian_model.Pclass{2}.Prior = model.Prior(i);
        bayesian_model.Prior(2) = model.Prior(i);
    else % muscle group.
        bayesian_model.Pclass{1}.Mean(:,count) = model.Mean(:,i);
        bayesian_model.Pclass{1}.Cov(:,:,count) = model.Cov(:,:,i);
        bayesian_model.Pclass{1}.fun = model.fun;
        bayesian_model.Pclass{1}.Prior(count) = model.Prior(i);
        bayesian_model.Prior(1) = bayesian_model.Prior(1) + ...
            model.Prior(i);
        count = count + 1;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bayesian_model = Generate3ClassGMM(model, nClusters)

classWSMeans = cat(1, model.Mean(2,:));
[temp,fatLabel] = max(classWSMeans);
[temp,muscleLabel] = min(classWSMeans);

clear temp;
count=1;
bayesian_model.Pclass = cell(1,3);
bayesian_model.Prior = zeros(1,3);
% Construct a 3-class model from the Gaussians.
for(i=1:nClusters)
    if(i==fatLabel) % interMF group.
        bayesian_model.Pclass{2}.Mean = model.Mean(:,i);
        bayesian_model.Pclass{2}.Cov = model.Cov(:,:,i);
        bayesian_model.Pclass{2}.fun = model.fun;
        bayesian_model.Pclass{2}.Prior = model.Prior(i);
        bayesian_model.Prior(2) = model.Prior(i);
    elseif(i==muscleLabel) % muscle group.
        bayesian_model.Pclass{1}.Mean = model.Mean(:,i);
        bayesian_model.Pclass{1}.Cov = model.Cov(:,:,i);
        bayesian_model.Pclass{1}.fun = model.fun;
        bayesian_model.Pclass{1}.Prior = model.Prior(i);
        bayesian_model.Prior(1) = model.Prior(i);
    else % intraMF
        bayesian_model.Pclass{3}.Mean(:,count) = model.Mean(:,i);
        bayesian_model.Pclass{3}.Cov(:,:,count) = model.Cov(:,:,i);
        bayesian_model.Pclass{3}.fun = model.fun;
        bayesian_model.Pclass{3}.Prior(count) = model.Prior(i);
        bayesian_model.Prior(3) = bayesian_model.Prior(2) + ...
            model.Prior(i);
        count = count + 1;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Samples = CopySamplesByCluster(dataMatrix, tissueIdx, nClusters)

for count=1:nClusters
    vectorIndex = find(tissueIdx==count);
    Samples{count}.FS = dataMatrix(vectorIndex,1);
    Samples{count}.WS = dataMatrix(vectorIndex,2);
end

end
