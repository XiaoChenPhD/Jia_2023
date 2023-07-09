% As one reviewer suggested, I calculated the windowed FCs from the left
% SFG cluster.
%
% modified from DPABI_TDA_Surf_run.m
% Xiao Chen 230705
% chenxiaophd@gmail.com

%% initialization
clear; clc;

load('Parameters_of_TDA.mat');% Cfg file of TDA analysis GUI
Cfg.StartingDirName = 'FunSurfWCFS';
Cfg.StartingDirForDCetc = Cfg.StartingDirName;
Cfg.FunctionalSessionNumber = 2;
Cfg.CalFC.ROIDefSurfLH{1,1} = 'DC_interaction_Cluster_Mask.gii';
Cfg.StartingDirName_Volume = ['FunVolu',Cfg.StartingDirName(8:end)];
Cfg.SubjectNum=length(Cfg.SubjectID);

if ~isfield(Cfg.CalFC,'ROIDefVolu')
        Cfg.CalFC.ROIDefVolu = {};
end
if ~isfield(Cfg.CalFC,'ROISelectedIndexVolu')
    Cfg.CalFC.ROISelectedIndexVolu = [];
end
if ~isfield(Cfg.CalFC,'ROISelectedIndexSurfLH')
    Cfg.CalFC.ROISelectedIndexSurfLH = [];
end
if ~isfield(Cfg.CalFC,'ROISelectedIndexSurfRH')
    Cfg.CalFC.ROISelectedIndexSurfRH = [];
end
                           
FunSessionPrefixSet = {'', 'S2_'};


%% Extract ROI Signals for Dynamic Functional Connectivity Analysis
for iFunSession=1:Cfg.FunctionalSessionNumber
    mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'TemporalDynamics',filesep,'ROISignals_',Cfg.StartingDirName]);
    mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'TemporalDynamics',filesep,'ROISignals_',Cfg.StartingDirName]);
    if (Cfg.IsProcessVolumeSpace==1)
        mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'ROISignals_',Cfg.StartingDirName_Volume]);
    end
    mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName]);

    %Extract the ROI time courses
    parfor i=1:Cfg.SubjectNum
        ROISignalsSurfLH=[];
        ROISignalsSurfRH=[];
        ROISignalsVolu=[];
        % Left Hemi
        if ~isempty(Cfg.CalFC.ROIDefSurfLH)
            DirName=dir(fullfile(Cfg.WorkingDir,[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],Cfg.SubjectID{i},'*fsaverage5_hemi-L*.func.gii'));
            for iFile=1:length(DirName)
                FileName=DirName(iFile).name;
                [ROISignalsSurfLH] = y_ExtractROISignal_Surf(fullfile(Cfg.WorkingDir,[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],Cfg.SubjectID{i},FileName), ...
                    Cfg.CalFC.ROIDefSurfLH, ...
                    [Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'TemporalDynamics',filesep,'ROISignals_',Cfg.StartingDirName,filesep,Cfg.SubjectID{i}], ...
                    '', ... % Will not restrict into the brain mask in extracting ROI signals
                    Cfg.CalFC.IsMultipleLabel,Cfg.CalFC.ROISelectedIndexSurfLH);
            end
        end

        % Right Hemi
        if ~isempty(Cfg.CalFC.ROIDefSurfRH)
            DirName=dir(fullfile(Cfg.WorkingDir,[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],Cfg.SubjectID{i},'*fsaverage5_hemi-R*.func.gii'));
            for iFile=1:length(DirName)
                FileName=DirName(iFile).name;
                [ROISignalsSurfRH] = y_ExtractROISignal_Surf(fullfile(Cfg.WorkingDir,[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],Cfg.SubjectID{i},FileName), ...
                    Cfg.CalFC.ROIDefSurfRH, ...
                    [Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'TemporalDynamics',filesep,'ROISignals_',Cfg.StartingDirName,filesep,Cfg.SubjectID{i}], ...
                    '', ... % Will not restrict into the brain mask in extracting ROI signals
                    Cfg.CalFC.IsMultipleLabel,Cfg.CalFC.ROISelectedIndexSurfRH);
            end
        end

        % Volume
        if ~isempty(Cfg.CalFC.ROIDefVolu)
            [ROISignalsVolu] = y_ExtractROISignal([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName_Volume,filesep,Cfg.SubjectID{i}], ...
            Cfg.CalFC.ROIDefVolu, ...
            [Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'ROISignals_',Cfg.StartingDirName_Volume,filesep,Cfg.SubjectID{i}], ...
            '', ... % Will not restrict into the brain mask in extracting ROI signals
            Cfg.CalFC.IsMultipleLabel,Cfg.CalFC.ROISelectedIndexVolu);
        end

        ROISignals = [ROISignalsSurfLH, ROISignalsSurfRH, ROISignalsVolu];
        y_CallSave([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROISignals_',Cfg.SubjectID{i},'.mat'], ROISignals, '');
        y_CallSave([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROISignals_',Cfg.SubjectID{i},'.txt'], ROISignals, ' ''-ASCII'', ''-DOUBLE'',''-TABS''');
        ROICorrelation = corrcoef(ROISignals);
        y_CallSave([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_',Cfg.SubjectID{i},'.mat'], ROICorrelation, '');
        y_CallSave([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_',Cfg.SubjectID{i},'.txt'], ROICorrelation, ' ''-ASCII'', ''-DOUBLE'',''-TABS''');
        ROICorrelation_FisherZ = 0.5 * log((1 + ROICorrelation)./(1- ROICorrelation));
        y_CallSave([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_FisherZ_',Cfg.SubjectID{i},'.mat'], ROICorrelation_FisherZ, '');
        y_CallSave([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROICorrelation_FisherZ_',Cfg.SubjectID{i},'.txt'], ROICorrelation_FisherZ, ' ''-ASCII'', ''-DOUBLE'',''-TABS''');
    end
end


%% Dynamic FC

for iFunSession=1:Cfg.FunctionalSessionNumber
    mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'TemporalDynamics',filesep,'TemporalDynamics4D',filesep,'FC_',Cfg.StartingDirName]);
    mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName]);
    mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'TemporalDynamics',filesep,'TemporalDynamics4D',filesep,'FC_',Cfg.StartingDirName]);
    mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName]);

    if (Cfg.IsProcessVolumeSpace==1)
        mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'TemporalDynamics4D',filesep,'FC_',Cfg.StartingDirName_Volume]);
        mkdir([Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName_Volume]);
    end
    for i=1:Cfg.SubjectNum
        ROIDef = {[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'ROISignals_SurfLHSurfRHVolu_',Cfg.StartingDirName,filesep, 'ROISignals_',Cfg.SubjectID{i},'.txt']};
        IsMultipleLabel = 1;

        % Left Hemi
        DirName=dir(fullfile(Cfg.WorkingDir,[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],Cfg.SubjectID{i},'*fsaverage5_hemi-L*.func.gii'));
        for iFile=1:length(DirName)
            FileName=DirName(iFile).name;
            InFiles = fullfile(Cfg.WorkingDir,[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],Cfg.SubjectID{i},FileName);
            OutFile = [Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'TemporalDynamics',filesep,'TemporalDynamics4D',filesep,'FC_',Cfg.StartingDirName,filesep,'FC_',Cfg.SubjectID{i},'.func.gii'];

            %[FCBrain_AllWindow, zFCBrain_AllWindow, GHeader] = y_SCA_Surf_Window(WindowSize, WindowStep, WindowType, AllVolume, ROIDef, OutputName, AMaskFilename, IsMultipleLabel, IsNeedDetrend, GHeader, CUTNUMBER)
            [FCBrain_AllWindow, zFCBrain_AllWindow, GHeader] = y_SCA_Surf_Window(Cfg.WindowSize, Cfg.WindowStep, Cfg.WindowType, InFiles, ROIDef, OutFile, Cfg.MaskFileSurfLH, IsMultipleLabel, Cfg.IsDetrend);

            %Calculate mean and std
            for iROI=1:size(zFCBrain_AllWindow,3)
                y_Write(squeeze(mean(zFCBrain_AllWindow(:,:,iROI),2)),GHeader,[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName,filesep,'MeanzFC_','ROI',num2str(iROI),'_',Cfg.SubjectID{i}]);
                y_Write(squeeze(std(zFCBrain_AllWindow(:,:,iROI),0,2)),GHeader,[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName,filesep,'StdzFC_','ROI',num2str(iROI),'_',Cfg.SubjectID{i}]);
                Temp = squeeze(std(zFCBrain_AllWindow(:,:,iROI),0,2)) ./ squeeze(mean(zFCBrain_AllWindow(:,:,iROI),2));
                Temp(find(isnan(Temp)))=0;
                y_Write(Temp,GHeader,[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfLH',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName,filesep,'CVzFC_','ROI',num2str(iROI),'_',Cfg.SubjectID{i}]);
            end
        end

        % Right Hemi
        DirName=dir(fullfile(Cfg.WorkingDir,[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],Cfg.SubjectID{i},'*fsaverage5_hemi-R*.func.gii'));
        for iFile=1:length(DirName)
            FileName=DirName(iFile).name;
            InFiles = fullfile(Cfg.WorkingDir,[FunSessionPrefixSet{iFunSession},Cfg.StartingDirName],Cfg.SubjectID{i},FileName);
            OutFile = [Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'TemporalDynamics',filesep,'TemporalDynamics4D',filesep,'FC_',Cfg.StartingDirName,filesep,'FC_',Cfg.SubjectID{i},'.func.gii'];

            %[FCBrain_AllWindow, zFCBrain_AllWindow, GHeader] = y_SCA_Surf_Window(WindowSize, WindowStep, WindowType, AllVolume, ROIDef, OutputName, AMaskFilename, IsMultipleLabel, IsNeedDetrend, GHeader, CUTNUMBER)
            [FCBrain_AllWindow, zFCBrain_AllWindow, GHeader] = y_SCA_Surf_Window(Cfg.WindowSize, Cfg.WindowStep, Cfg.WindowType, InFiles, ROIDef, OutFile, Cfg.MaskFileSurfRH, IsMultipleLabel, Cfg.IsDetrend);

            %Calculate mean and std
            for iROI=1:size(zFCBrain_AllWindow,3)
                y_Write(squeeze(mean(zFCBrain_AllWindow(:,:,iROI),2)),GHeader,[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName,filesep,'MeanzFC_','ROI',num2str(iROI),'_',Cfg.SubjectID{i}]);
                y_Write(squeeze(std(zFCBrain_AllWindow(:,:,iROI),0,2)),GHeader,[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName,filesep,'StdzFC_','ROI',num2str(iROI),'_',Cfg.SubjectID{i}]);
                Temp = squeeze(std(zFCBrain_AllWindow(:,:,iROI),0,2)) ./ squeeze(mean(zFCBrain_AllWindow(:,:,iROI),2));
                Temp(find(isnan(Temp)))=0;
                y_Write(Temp,GHeader,[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunSurfRH',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName,filesep,'CVzFC_','ROI',num2str(iROI),'_',Cfg.SubjectID{i}]);
            end
        end

        % Volume
        if (Cfg.IsProcessVolumeSpace==1)
            InFiles = [Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},Cfg.StartingDirName_Volume,filesep,Cfg.SubjectID{i}];
            OutFile = [Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'TemporalDynamics4D',filesep,'FC_',Cfg.StartingDirName_Volume,filesep,'FC_',Cfg.SubjectID{i}];

            %[FCBrain_AllWindow, zFCBrain_AllWindow, Header] = y_SCA_Window(WindowSize, WindowStep, WindowType, AllVolume, ROIDef, OutputName, MaskData, IsMultipleLabel, IsNeedDetrend, Band, TR, TemporalMask, ScrubbingMethod, ScrubbingTiming, Header, CUTNUMBER)
            [FCBrain_AllWindow, zFCBrain_AllWindow, Header] = y_SCA_Window(Cfg.WindowSize, Cfg.WindowStep, Cfg.WindowType, InFiles, ROIDef, OutFile, Cfg.MaskFileVolu, IsMultipleLabel, Cfg.IsDetrend);

            %Calculate mean and std
            for iROI=1:size(zFCBrain_AllWindow,5)
                y_Write(squeeze(mean(zFCBrain_AllWindow(:,:,:,:,iROI),4)),Header,[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName_Volume,filesep,'MeanzFC_','ROI',num2str(iROI),'_',Cfg.SubjectID{i}]);
                y_Write(squeeze(std(zFCBrain_AllWindow(:,:,:,:,iROI),0,4)),Header,[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName_Volume,filesep,'StdzFC_','ROI',num2str(iROI),'_',Cfg.SubjectID{i}]);
                Temp = squeeze(std(zFCBrain_AllWindow(:,:,:,:,iROI),0,4)) ./ squeeze(mean(zFCBrain_AllWindow(:,:,:,:,iROI),4));
                Temp(find(isnan(Temp)))=0;
                y_Write(Temp,Header,[Cfg.WorkingDir,filesep,FunSessionPrefixSet{iFunSession},'Results',filesep,'FunVolu',filesep,'TemporalDynamics',filesep,'TemporalDynamicsMetrics',filesep,'FC_',Cfg.StartingDirName_Volume,filesep,'CVzFC_','ROI',num2str(iROI),'_',Cfg.SubjectID{i}]);
            end
        end
    end
end

