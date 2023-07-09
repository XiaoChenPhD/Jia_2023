% As the reviewer requested, I analysed the sliding windowed FCs
% This script is for doing stats
%
% Xiao Chen 230706
% chenxiaophd@gmail.com

%% initialization
clear; clc;

work_dir = 'DynamicFC';
data_dir = 'DPABISurf Preprocessing Files';
template_dir = 'DPABISurf surf template directory';

hemisp_name_set = {'lh','rh'};
hemisp_set = {'LH','RH'};
measure_set = {'CV', 'Mean', 'Std'};
condition_set = {'rum', 'dis'};
condition_prefix = {'', 'S2_'};

% load clean subjects
load('sub_info_clean.mat');

%% Organize metric files
% surface
% MDD
for iCondition = 1:length(condition_set) %condition
    for iHem = 1:2 %hemisphere
        for iMeasure = 1:length(measure_set) %brain metrics
            for i = 1:length(find(Dx == 1)) %subject
                origin_file = [data_dir,'/',condition_prefix{iCondition}, ...
                          'Results/FunSurf',hemisp_set{iHem}, ...
                          '/TemporalDynamics/TemporalDynamicsMetrics/FC_FunSurfWCFS/', ...
                          measure_set{iMeasure}, ...
                              'zFC_ROI1_sub-',sub_list{i},'.gii'];
                target_dir = [work_dir,'/metric_maps/',condition_set{iCondition}, ...
                             '/MDD_FunSurf',hemisp_set{iHem}, ...
                             '/',measure_set{iMeasure}];
                if ~exist(target_dir, 'dir'); mkdir(target_dir);end
                target_file = [target_dir,'/',measure_set{iMeasure}, ...
                              'zFC_ROI1_sub-',sub_list{i},'.gii'];
                copyfile(origin_file,target_file);
            end
        end
    end
end

% HC
for iCondition = 1:length(condition_set) %condition
    for iHem = 1:2 %hemisphere
        for iMeasure = 1:length(measure_set) %brain metrics
            for i = length(find(Dx == 1))+1:length(sub_list) %subject
                origin_file = [data_dir,'/',condition_prefix{iCondition}, ...
                          'Results/FunSurf',hemisp_set{iHem}, ...
                          '/TemporalDynamics/TemporalDynamicsMetrics/FC_FunSurfWCFS/', ...
                          measure_set{iMeasure}, ...
                              'zFC_ROI1_sub-',sub_list{i},'.gii'];
                target_dir = [work_dir,'/metric_maps/',condition_set{iCondition}, ...
                             '/HC_FunSurf',hemisp_set{iHem}, ...
                             '/',measure_set{iMeasure}];
                if ~exist(target_dir, 'dir'); mkdir(target_dir);end
                target_file = [target_dir,'/',measure_set{iMeasure}, ...
                              'zFC_ROI1_sub-',sub_list{i},'.gii'];
                copyfile(origin_file,target_file);
            end
        end
    end
end

% subcortical
% MDD
for iCondition = 1:length(condition_set) %condition
    for iMeasure = 1:length(measure_set) %brain metrics
        for i = 1:length(find(Dx == 1)) %subject
            origin_file = [data_dir,'/',condition_prefix{iCondition}, ...
                      'Results/FunVolu', ...
                      '/TemporalDynamics/TemporalDynamicsMetrics/FC_FunVoluWCFS/', ...
                      measure_set{iMeasure}, ...
                          'zFC_ROI1_sub-',sub_list{i},'.nii'];
            target_dir = [work_dir,'/metric_maps/',condition_set{iCondition}, ...
                         '/MDD_FunVolu/', ...
                         measure_set{iMeasure}];
            if ~exist(target_dir, 'dir'); mkdir(target_dir);end
            target_file = [target_dir,'/',measure_set{iMeasure}, ...
                          'zFC_ROI1_sub-',sub_list{i},'.nii'];
            copyfile(origin_file,target_file);
        end
    end
end

%HC
for iCondition = 1:length(condition_set) %condition
    for iMeasure = 1:length(measure_set) %brain metrics
        for i = length(find(Dx == 1))+1:length(sub_list) %subject
            origin_file = [data_dir,'/',condition_prefix{iCondition}, ...
                      'Results/FunVolu', ...
                      '/TemporalDynamics/TemporalDynamicsMetrics/FC_FunVoluWCFS/', ...
                      measure_set{iMeasure}, ...
                          'zFC_ROI1_sub-',sub_list{i},'.nii'];
            target_dir = [work_dir,'/metric_maps/',condition_set{iCondition}, ...
                         '/HC_FunVolu/', ...
                         measure_set{iMeasure}];
            if ~exist(target_dir, 'dir'); mkdir(target_dir);end
            target_file = [target_dir,'/',measure_set{iMeasure}, ...
                          'zFC_ROI1_sub-',sub_list{i},'.nii'];
            copyfile(origin_file,target_file);
        end
    end
end

%% do mixed effect analysis, vertex-wised
PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
PALMSettings.SavePermutations = 0;

% left hemisphere
PALMSettings.SurfFile = [template_dir,'/fsaverage5_lh_white.surf.gii'];
PALMSettings.SurfAreaFile = [template_dir,'/fsaverage5_lh_white_avg.area.gii'];
MaskFile =  [template_dir,'/fsaverage5_lh_cortex.label.gii'];
for iMeasure = 1:length(measure_set)
        DependentDir{1,1} = [work_dir,'/metric_maps/rum/MDD_FunSurfLH/',measure_set{iMeasure}];                          
        DependentDir{2,1} = [work_dir,'/metric_maps/dis/MDD_FunSurfLH/', measure_set{iMeasure}];                           
        DependentDir{3,1} = [work_dir,'/metric_maps/rum/HC_FunSurfLH/', measure_set{iMeasure}];                          
        DependentDir{4,1} = [work_dir,'/metric_maps/dis/HC_FunSurfLH/', measure_set{iMeasure}];
                                    

        OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
        OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
        OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
        OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

        OutputDir = [work_dir,'/stats'];
        if ~exist(OutputDir, 'dir'); mkdir(OutputDir); end
        OutputName = [OutputDir,'/',measure_set{iMeasure},'_lh'];
        y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);
end

% right hemisphere
PALMSettings.SurfFile = [template_dir,'/fsaverage5_rh_white.surf.gii'];
PALMSettings.SurfAreaFile = [template_dir,'/fsaverage5_rh_white_avg.area.gii'];
MaskFile =  [template_dir,'/fsaverage5_rh_cortex.label.gii'];
for iMeasure = 1:length(measure_set)
        DependentDir{1,1} = [work_dir,'/metric_maps/rum/MDD_FunSurfRH/',measure_set{iMeasure}];
        DependentDir{2,1} = [work_dir,'/metric_maps/dis/MDD_FunSurfRH/',measure_set{iMeasure}];
        DependentDir{3,1} = [work_dir,'/metric_maps/rum/HC_FunSurfRH/',measure_set{iMeasure}];
        DependentDir{4,1} = [work_dir,'/metric_maps/dis/HC_FunSurfRH/',measure_set{iMeasure}];

        OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
        OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
        OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
        OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

        OutputDir = [work_dir,'/stats'];
        if ~exist(OutputDir, 'dir'); mkdir(OutputDir); end
        OutputName = [OutputDir,'/',measure_set{iMeasure},'_rh'];
        y_MixedEffectsAnalysis_Image(DependentDir,OutputName,MaskFile,[],OtherCovariates, PALMSettings);
end


%% do one-sample t test, vertex
PALMSettings.nPerm = 5000;
PALMSettings.ClusterInference=0;
PALMSettings.ClusterFormingThreshold=2.3;
PALMSettings.TFCE=1;
PALMSettings.FDR=0;
PALMSettings.TwoTailed=1;
PALMSettings.AccelerationMethod='NoAcceleration'; % or 'tail', 'gamma', 'negbin', 'lowrank', 'noperm'
PALMSettings.SavePermutations = 0;

group_set = {'MDD', 'HC'};
Base = 0;

for iGroup = 1:2
    OtherCovariates{1,1} = HeadMotion(Dx == iGroup);
    for iCondition = 1:length(condition_set)
        % left hemisphere
        PALMSettings.SurfFile = [template_dir,'/fsaverage5_lh_white.surf.gii'];
        PALMSettings.SurfAreaFile = [template_dir,'/fsaverage5_lh_white_avg.area.gii'];
        MaskFile =  [template_dir,'/fsaverage5_lh_cortex.label.gii'];
        DependentDirs{1,1} = [work_dir, '/metric_maps/', condition_set{iCondition}, ...
                             '/', group_set{iGroup},'_FunSurfLH/CV'];
        OutputName = [work_dir, '/stats_one_sample/', group_set{iGroup}, '_', ...
                                condition_set{iCondition}, '_lh'];
        y_TTest1_Image(DependentDirs,OutputName,MaskFile,[],OtherCovariates,Base,PALMSettings)

        % right hemisphere
        PALMSettings.SurfFile = [template_dir,'/fsaverage5_rh_white.surf.gii'];
        PALMSettings.SurfAreaFile = [template_dir,'/fsaverage5_rh_white_avg.area.gii'];
        MaskFile =  [template_dir,'/fsaverage5_rh_cortex.label.gii'];
        DependentDirs{1,1} = [work_dir, '/metric_maps/', condition_set{iCondition}, ...
                             '/', group_set{iGroup},'_FunSurfRH/CV'];
        OutputName = [work_dir, '/stats_one_sample/', group_set{iGroup}, '_', ...
                                condition_set{iCondition}, '_rh'];
        y_TTest1_Image(DependentDirs,OutputName,MaskFile,[],OtherCovariates,Base,PALMSettings)
    end
end

%% do network-wised analyses
network_temp_dir = 'directory storing Schaefer template';

%% extract signals
TopOutputDir = [work_dir, '/dynamic_network_extracted'];
if ~exist(TopOutputDir, 'dir'); mkdir(TopOutputDir); end
% surface
IsMultipleLabel = 1;
GHeader = [];
CUTNUMBER = 1;
ROIDef = {};
IsNeedDetrend = 0;
for iCondition = 1:length(condition_set)
    for iHem = 1:2
        for iMeasure = 1:length(measure_set)
            ROIDef{1} = [network_temp_dir, '/fsaverage5_', hemisp_name_set{iHem}, ...
                                                '_Schaefer2018_400Parcels_17Networks_order.label.gii'];
            % MDD
            AllVolume = [work_dir, '/metric_maps/', condition_set{iCondition}, ...
                         '/MDD_FunSurf', hemisp_set{iHem},'/', measure_set{iMeasure}];
            OutputName = [TopOutputDir, '/MDD_', condition_set{iCondition}, '_', ...
                          measure_set{iMeasure}, '_',hemisp_name_set{iHem}];
            AMaskFilename = [template_dir, '/fsaverage5_', hemisp_name_set{iHem}, '_cortex.label.gii'];
            [ROISignals] = y_ExtractROISignal_Surf(AllVolume, ROIDef, OutputName, AMaskFilename, ...
                                                                  IsMultipleLabel, GHeader, CUTNUMBER);
                                                              
            %HC
            AllVolume = [work_dir, '/metric_maps/', condition_set{iCondition}, ...
                         '/HC_FunSurf', hemisp_set{iHem},'/', measure_set{iMeasure}];
            OutputName = [TopOutputDir, '/HC_', condition_set{iCondition}, '_', ...
                          measure_set{iMeasure}, '_',hemisp_name_set{iHem}];
            AMaskFilename = [template_dir, '/fsaverage5_', hemisp_name_set{iHem}, '_cortex.label.gii'];
            [ROISignals] = y_ExtractROISignal_Surf(AllVolume, ROIDef, OutputName, AMaskFilename, ...
                                                                  IsMultipleLabel, GHeader, CUTNUMBER);
        end
    end
end

%%
% subcortical
IsNeedDetrend = 0;
TemporalMask = [];

for iCondition = 1:length(condition_set)
    for iMeasure = 1:length(measure_set)
        ROIDef{1} = [network_temp_dir, '/Tian_Subcortex_S4_3T_2009cAsym.nii'];
        % MDD
        AllVolume = [work_dir, '/metric_maps/', condition_set{iCondition}, ...
                     '/MDD_FunVolu/', measure_set{iMeasure}];
        OutputName = [TopOutputDir, '/MDD_', condition_set{iCondition}, '_', ...
                      measure_set{iMeasure}, '_subcortical'];
        MaskData = [data_dir,'/Masks/AllResampled_BrainMask_05_91x109x91.nii'];
        [ROISignals] = y_ExtractROISignal(AllVolume, ROIDef, OutputName, MaskData, ...
                                          IsMultipleLabel, [], IsNeedDetrend, [], ...
                                          TemporalMask, [], [], [], CUTNUMBER);
                                      
        % HC
        AllVolume = [work_dir, '/metric_maps/', condition_set{iCondition}, ...
                     '/HC_FunVolu/', measure_set{iMeasure}];
        OutputName = [TopOutputDir, '/HC_', condition_set{iCondition}, '_', ...
                      measure_set{iMeasure}, '_subcortical'];
        MaskData = [data_dir,'/Masks/AllResampled_BrainMask_05_91x109x91.nii'];
        [ROISignals] = y_ExtractROISignal(AllVolume, ROIDef, OutputName, MaskData, ...
                                          IsMultipleLabel, [], IsNeedDetrend, [], ...
                                          TemporalMask, [], [], [], CUTNUMBER);
    end
end

%% get network FC
%load template info file, could be found in ~/template
load('DPABISurf_Schaefer2018_400_17Networks_Tian2020_54_Info438_V3.mat');
conditon4loop_set = {'rum', 'rum', 'dis', 'dis'};
group4loop_set = {'MDD', 'MDD', 'HC', 'HC'};
            
% network_FC: 1 x 4 cells, 4 conditions in condition_set,
% each cell has n x 10 matrixs, 10 networks defined as the abovementrioned
% .mat file
for iMeasure = 1:length(measure_set)
    network_FC = {};
    for iCondition = 1:length(conditon4loop_set)
        load([TopOutputDir,'/ROISignals_', group4loop_set{iCondition}, '_', ...
                conditon4loop_set{iCondition}, '_', ...
                measure_set{iMeasure}, '_lh.mat']);
        ROISignals_lh = ROISignals;
        
        load([TopOutputDir,'/ROISignals_', group4loop_set{iCondition}, '_', ...
                conditon4loop_set{iCondition}, '_', ...
                measure_set{iMeasure}, '_rh.mat']);
        ROISignals_rh = ROISignals;
        
        load([TopOutputDir,'/ROISignals_', group4loop_set{iCondition}, '_', ...
                conditon4loop_set{iCondition}, '_', ...
                measure_set{iMeasure}, '_subcortical.mat']);
        ROISignals_subcortical = ROISignals;
        
        ROISignals_full{iCondition} = [ROISignals_lh,ROISignals_rh,ROISignals_subcortical];
        ROISignals_full_438{iCondition} = ROISignals_full{iCondition}(:, ROIIndex_181_Schafer_Tian);
        network_label = DPABISurf_Schaefer2018_400_Tian2020_54_YeoNetwork';
        for iNetwork = unique(network_label)
            network_FC{iCondition}(:,iNetwork) = mean(ROISignals_full_438{iCondition}(:, find(network_label == iNetwork)), 2);
        end
    end
    save([work_dir,'/', measure_set{iMeasure}, '_network_dynamic.mat'], 'network_FC');
end


%% do mixed effect analysis
iMeasure = 2;
load([work_dir,'/', measure_set{iMeasure}, '_network_dynamic.mat']);

MDD_Rum_n = size(network_FC{1},1);
MDD_Dis_n = size(network_FC{2},1);
HC_Rum_n = size(network_FC{3},1);
HC_Dis_n = size(network_FC{4},1);

BetweenSubjectFactor1 = ones(MDD_Rum_n+MDD_Dis_n,1);
BetweenSubjectFactor2 = -1*ones(HC_Rum_n+HC_Dis_n,1);
BetweenSubjectFactor = [BetweenSubjectFactor1;BetweenSubjectFactor2];
WithinSubjectFactor = [ones(MDD_Rum_n,1);-1*ones(MDD_Dis_n,1);ones(HC_Rum_n,1);-1*ones(HC_Dis_n,1)];
nSubject = [MDD_Rum_n MDD_Dis_n HC_Rum_n HC_Dis_n];

OtherCovariates{1,1} = HeadMotion(1:length(find(Dx == 1)),1);
OtherCovariates{2,1} = HeadMotion(1:length(find(Dx == 1)),2);
OtherCovariates{3,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),1);
OtherCovariates{4,1} = HeadMotion(length(find(Dx == 1))+1:length(sub_list),2);

SubjectRegressorsAll = [];
for i=1:4
    SubjectRegressorsAll=[SubjectRegressorsAll;(ceil(i/2)-1)*10000+([1:nSubject(i)]')];
end
SubIndex=unique(SubjectRegressorsAll);
nSub=length(SubIndex);
SubjectRegressors=[];
for i=1:nSub
    SubjectRegressors(:,i) = zeros(size(SubjectRegressorsAll));
    SubjectRegressors(SubjectRegressorsAll==SubIndex(i),i) = 1;
end
Interaction = WithinSubjectFactor.*BetweenSubjectFactor;
HeadMotion_Regressor = [OtherCovariates{1};OtherCovariates{2};OtherCovariates{3};OtherCovariates{4}];
AllCov = [WithinSubjectFactor,Interaction,SubjectRegressors,HeadMotion_Regressor];
Contrast = zeros(1,size(AllCov,2));
TF_Flag = 'T';

% condition effect
Contrast(1)=1;
for iNetwork = unique(network_label)
    y = [network_FC{1}(:,iNetwork);network_FC{2}(:,iNetwork);network_FC{3}(:,iNetwork);network_FC{4}(:,iNetwork)];
    [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,AllCov,Contrast,TF_Flag);
    Stats_Mixed_Condition{iNetwork,1} = YeoSCNetwork_Label{iNetwork,2};
    Stats_Mixed_Condition{iNetwork,2} = TF_ForContrast;
    Stats_Mixed_Condition{iNetwork,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
    Stats_Mixed_Condition{iNetwork,4} = 2*tcdf(-abs(TF_ForContrast),nSub-2)*9;
    Stats_Mixed_Condition{iNetwork,5} = Cohen_f2;
end

% interaction effect
Contrast(1)=0;
Contrast(2)=1;
for iNetwork = unique(network_label)
    y = [network_FC{1}(:,iNetwork);network_FC{2}(:,iNetwork);network_FC{3}(:,iNetwork);network_FC{4}(:,iNetwork)];
    [b,r,SSE,SSR, T, TF_ForContrast, Cohen_f2] = y_regress_ss(y,AllCov,Contrast,TF_Flag);
    Stats_Mixed_Interaction{iNetwork,1} = YeoSCNetwork_Label{iNetwork,2};
    Stats_Mixed_Interaction{iNetwork,2} = TF_ForContrast;
    Stats_Mixed_Interaction{iNetwork,3} = 2*tcdf(-abs(TF_ForContrast),nSub-2);
    Stats_Mixed_Interaction{iNetwork,4} = 2*tcdf(-abs(TF_ForContrast),nSub-2)*9;
    Stats_Mixed_Interaction{iNetwork,5} = Cohen_f2;
    p_raw(iNetwork) =  2*tcdf(-abs(TF_ForContrast),nSub-2);
end