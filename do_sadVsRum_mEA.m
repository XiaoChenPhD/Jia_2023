% do mixed effect analysis on degree centrality on sad vs. rum, as the
% reviewer suggested
%
% 
% Xiao Chen 230228
% chenxiaophd@gmail.com

%% initialization
clear; clc;
work_dir = 'working directory';
data_dir = 'DPABISurf Preprocessing Files';
TemplateDir = 'DPABISurf surf template directory';

load([work_dir, '/sub_info_clean.mat']);

HemisphereNameSet = {'lh','rh'};
HemisphereSet = {'LH','RH'};

MeasureSet = {'DegreeCentrality'};
PipesuffixSet = {'_FunSurfWCF'};
SuffixSet = {'_Bilateral_PositiveWeightedSumBrain'};

%% Organize metric files
% Rumination
% MDD
for i = 1:length(find(Dx == 1))
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Task/ResultsS/FunSurf',HemisphereSet{iHem},'/', ...
                              MeasureSet{iMeasure},PipesuffixSet{iMeasure}, ...
                              '/sz',MeasureSet{iMeasure},SuffixSet{iMeasure}, ...
                              '_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/rum/MDD_FunSurf',HemisphereSet{iHem}, ...
                              '/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
            if ~exist(TargetDir); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure}, ...
                            SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

%HC
for i = length(find(Dx == 1))+1:length(sub_list)
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Task/ResultsS/FunSurf',HemisphereSet{iHem},'/', ...
                                 MeasureSet{iMeasure},PipesuffixSet{iMeasure}, ...
                                 '/sz',MeasureSet{iMeasure},SuffixSet{iMeasure}, ...
                                 '_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/rum/HC_FunSurf',HemisphereSet{iHem}, ...
                                 '/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
            if ~exist(TargetDir); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure}, ...
                                  '_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

% Sad
%MDD
for i = 1:length(find(Dx == 1))
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Sad/ResultsS/FunSurf',HemisphereSet{iHem},'/', ...
                             MeasureSet{iMeasure},PipesuffixSet{iMeasure}, ...
                              '/sz',MeasureSet{iMeasure}, ...
                              SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/sad/MDD_FunSurf',HemisphereSet{iHem}, ...
                              '/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
            if ~exist(TargetDir); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure},SuffixSet{iMeasure}, ...
                            '_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

%HC
for i = length(find(Dx == 1))+1:length(sub_list)
    for iHem = 1:2
        for iMeasure = 1:length(MeasureSet)
            OriginFile = [data_dir,'/Sad/ResultsS/FunSurf',HemisphereSet{iHem}, ...
                                    '/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}, ...
                                    '/sz',MeasureSet{iMeasure},SuffixSet{iMeasure}, ...
                                    '_sub-',sub_list{i},'.func.gii'];
            TargetDir = [work_dir,'/sad/HC_FunSurf',HemisphereSet{iHem},'/', ...
                                     MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
            if ~exist(TargetDir); mkdir(TargetDir);end
            TargetFile = [TargetDir,'/sz',MeasureSet{iMeasure}, ...
                                   SuffixSet{iMeasure},'_sub-',sub_list{i},'.func.gii'];
            copyfile(OriginFile,TargetFile);
        end
    end
end

%% do mixed effect analysis
PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
PALMSettings.SavePermutations = 0;

% left hemisphere
PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_lh_white.surf.gii'];
PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_lh_white_avg.area.gii'];
MaskFile =  [TemplateDir,'/fsaverage5_lh_cortex.label.gii'];
for iMeasure = 1:length(MeasureSet)
        DependentDir{1,1} = [work_dir,'/rum/MDD_FunSurfLH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{2,1} = [work_dir,'/sad/MDD_FunSurfLH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{3,1} = [work_dir,'/rum/HC_FunSurfLH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{4,1} = [work_dir,'/sad/HC_FunSurfLH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];

        OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
        OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),3);
        OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
        OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),3);

        OutputDir = [work_dir,'/stats'];
        if ~exist(OutputDir); mkdir(OutputDir); end
        OutputName = [OutputDir,'/',MeasureSet{iMeasure},PipesuffixSet{iMeasure},'_lh'];
        y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);
end

% right hemisphere
PALMSettings.SurfFile = [TemplateDir,'/fsaverage5_rh_white.surf.gii'];
PALMSettings.SurfAreaFile = [TemplateDir,'/fsaverage5_rh_white_avg.area.gii'];
MaskFile =  [TemplateDir,'/fsaverage5_rh_cortex.label.gii'];
for iMeasure = 1:length(MeasureSet)
        DependentDir{1,1} = [work_dir,'/rum/MDD_FunSurfRH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{2,1} = [work_dir,'/sad/MDD_FunSurfRH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{3,1} = [work_dir,'/rum/HC_FunSurfRH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];
        DependentDir{4,1} = [work_dir,'/sad/HC_FunSurfRH/',MeasureSet{iMeasure},PipesuffixSet{iMeasure}];

        OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
        OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),3);
        OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
        OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),3);

        OutputDir = [work_dir,'/stats'];
        if ~exist(OutputDir); mkdir(OutputDir); end
        OutputName = [OutputDir,'/',MeasureSet{iMeasure},PipesuffixSet{iMeasure},'_rh'];
        y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);
end