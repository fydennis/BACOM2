clc,close all;
% IO parameters
ori_BACOM_result_dir = [pwd,'\Samples\']; % where original BACOM1 results are stored
SampleListFileDir = [pwd,'\'];
SampleListFileName = 'SampleList.txt';
TxtSaveDir = '';
TxtSaveFileName = 'BACOM2_Results.txt';

% some platform para.
% ChrLoc = dlmread('\ChrMarkerLocation.txt'); % notice: Line# N+1 corresponds to Chr. N
ChrLoc = [0;
146401;
300064;
427830;
548126;
663798;
776623;
877619;
975896;
1058064;
1151656;
1241181;
1328502;
1394569;
1451672;
1505228;
1559410;
1606042;
1658135;
1688434;
1732062;
1757173;
1781657];
MarkerSize = ChrLoc(23);

% for debug use
NeedReadFile = 1;
NeedBuidDataMtx = 1;
NeedSegmentation = 1;    
NeedCorrectCrossTalk = 1;
 
% genotype portion setting
N_genotypeAB_portion = 0.22; %AB
N_genotypeAABB_portion = 0.66;% AA/BB
T_genotypeAB_portion = 0.15;
T_genotypeAABB_portion = 0.62;

% thresholds
AlleleLowThres = -1;
CNLowThres = -2;
LongABThres = 30;
SegConvMask = pdf('norm',-0.25:0.01:0.25,0,0.05);
AmpThres = 2.3;
DelThres = 1.7;
RobustStdMultiplier = 2.5;
SlideWindowSize = 201;
SlideWindowHalfSize = floor(SlideWindowSize/2);
CN_DiffThres = 1.5E4;
NeuLengthThres = 1E3;
LocNormalizeThres = 0.15;
LongSegThres = 500;

fCompOut = fopen([TxtSaveDir,TxtSaveFileName],'wb');
fprintf(fCompOut,'Sample_Name\tTumor_Purity\tAverage_Ploidy\n');
%% 
% list of tumor sample names
SampleList = textread([SampleListFileDir, SampleListFileName],'%s');
SampleSize = size(SampleList,1);
for SampleIdx = 1:SampleSize
SampleName = SampleList{SampleIdx}(1:end-4);
%% Read raw intensity data
if NeedReadFile > 0 
    IntnDataFin = fopen([ori_BACOM_result_dir,SampleName,'.CEL_outputdata.csv']);
    if IntnDataFin<0
        fprintf(fCompOut,'%s\tNA\tNA\n',SampleName);
        continue;
    else        
        cnMtx = textscan(IntnDataFin, '%s %s %u %f %f %f %f %u', 'delimiter', ',');
        fclose(IntnDataFin);
    end
    
    if ~exist('MarkerNameList','var')
        ChrList= cnMtx{1};
        MarkerNameList = cnMtx{2};
        MarkerPosList  = cnMtx{3};
        ChrList= ChrList(1:ChrLoc(23));   
        MarkerNameList = MarkerNameList(1:ChrLoc(23));    
        MarkerChrList = str2double(ChrList(1:ChrLoc(23)));
        MarkerPosList  = double(MarkerPosList(1:ChrLoc(23)));
        MarkerTypeList = true(MarkerSize,1); 
        for idx = 1:MarkerSize
            MarkerTypeList(idx) = (MarkerNameList{idx}(1)=='S'); % SNP or CN probes 
        end
        SMarkerSize = sum(MarkerTypeList);
        CMarkerSize = sum(~MarkerTypeList);
    end
    % exclude sex chromosomes
    ITAList = cnMtx{4}(1:ChrLoc(23)); %intensity, tumor sample, allele A
    ITBList = cnMtx{5}(1:ChrLoc(23)); %intensity, tumor sample, allele B
    INAList = cnMtx{6}(1:ChrLoc(23)); %intensity, normal sample, allele A
    INBList = cnMtx{7}(1:ChrLoc(23)); %intensity, normal sample, allele B 
    
    % only interest in SNP probes
    IntensityMtxS = cat(2,MarkerChrList(MarkerTypeList),MarkerPosList(MarkerTypeList),ITAList(MarkerTypeList),ITBList(MarkerTypeList),INAList(MarkerTypeList),INBList(MarkerTypeList));
    IntensityMtxC = cat(2,MarkerChrList(~MarkerTypeList),MarkerPosList(~MarkerTypeList),ITAList(~MarkerTypeList),INAList(~MarkerTypeList));
end

%% preprocessing
if NeedBuidDataMtx == 1        
    ITAListS = IntensityMtxS(:,3);
    ITBListS = IntensityMtxS(:,4);
    INAListS = IntensityMtxS(:,5);
    INBListS = IntensityMtxS(:,6);
    
    ITListC = IntensityMtxC(:,3);
    INListC = IntensityMtxC(:,4);
    
%% background noise offset removal
    ItensitySteps = 10:20:1000;
    Ht = hist([ITAListS;ITBListS],ItensitySteps);
    Hn = hist([INAListS;INBListS],ItensitySteps); 
    Esti_Nt = ItensitySteps(find(Ht(2:end)>5E2,1,'first')+2);
    Esti_Nn = ItensitySteps(find(Hn(2:end)>5E2,1,'first')+2);    

    ITAListS = ITAListS-Esti_Nt;
    ITBListS = ITBListS-Esti_Nt;
    INAListS = INAListS-Esti_Nn;
    INBListS = INBListS-Esti_Nn;
    ITListC = ITListC-Esti_Nt;
    INListC = INListC-Esti_Nn;
    
%% Genotyping    
    N_genotypeAB_quantity = floor(SMarkerSize*N_genotypeAB_portion);
    N_genotypeAABB_quantity = floor(SMarkerSize*N_genotypeAABB_portion);
    
    NGenotypeMtx = cat(2,(1:SMarkerSize)',abs(log(INAListS)-log(INBListS)),zeros(SMarkerSize,1)); % distance to 45 degree line
    NGenotypeMtx = sortrows(NGenotypeMtx,2);
    NGenotypeMtx(1:N_genotypeAB_quantity,3) = ones(N_genotypeAB_quantity,1); % assign genotype AB to those with smallest distances 
    NGenotypeMtx(end-N_genotypeAABB_quantity+1:end,3) = -ones(N_genotypeAABB_quantity,1); % assign genotype AABB to those with largest distances 
    NGenotypeMtx = sortrows(NGenotypeMtx,1);
    N_ABGenotype = (NGenotypeMtx(:,3)==1); % boolean indicater of whether AB genotype
    N_AABBGenotype = (NGenotypeMtx(:,3)==-1);
    
    T_genotypeAB_quantity = floor(SMarkerSize*T_genotypeAB_portion);
    T_genotypeAABB_quantity = floor(SMarkerSize*T_genotypeAABB_portion);
    TGenotypeMtx = cat(2,(1:SMarkerSize)',abs(log(ITAListS)-log(ITBListS)),zeros(SMarkerSize,1));
    TGenotypeMtx = sortrows(TGenotypeMtx,2);
    TGenotypeMtx(1:T_genotypeAB_quantity,3) = ones(T_genotypeAB_quantity,1);
    TGenotypeMtx(end-T_genotypeAABB_quantity+1:end,3) = -ones(T_genotypeAABB_quantity,1);
    TGenotypeMtx = sortrows(TGenotypeMtx,1);
    T_ABGenotype = (TGenotypeMtx(:,3)==1);
    T_AABBGenotype = (TGenotypeMtx(:,3)==-1);
     
%% Cross-talk estimation and correction
    N_CrosstalkList = min(INAListS(N_AABBGenotype),INBListS(N_AABBGenotype))./max(INAListS(N_AABBGenotype),INBListS(N_AABBGenotype));  
    N_Crosstalk = mean(N_CrosstalkList(~isnan(N_CrosstalkList) & ~isinf(N_CrosstalkList)));    
    N_Crosstalk = N_Crosstalk/(1+N_Crosstalk);

    T_CrosstalkList = min(ITAListS(T_AABBGenotype),ITBListS(T_AABBGenotype))./max(ITAListS(T_AABBGenotype),ITBListS(T_AABBGenotype));
    T_Crosstalk = mean(T_CrosstalkList(~isnan(T_CrosstalkList) & ~isinf(T_CrosstalkList)));
    T_Crosstalk = T_Crosstalk/(1+T_Crosstalk);

    ITS = 1/(1-T_Crosstalk*2)*[1-T_Crosstalk,-T_Crosstalk;-T_Crosstalk,1-T_Crosstalk]*[ITAListS';ITBListS'];
    INS = 1/(1-N_Crosstalk*2)*[1-N_Crosstalk,-N_Crosstalk;-N_Crosstalk,1-N_Crosstalk]*[INAListS';INBListS'];
    
    if NeedCorrectCrossTalk ~= 0
        ITAListS = ITS(1,:)';
        ITBListS = ITS(2,:)';
        INAListS = INS(1,:)';
        INBListS = INS(2,:)';
    end
    
%% generate allelic copy nubmer signals
    alleleAS = ITAListS./(INAListS+INBListS)*2;
    alleleBS = ITBListS./(INAListS+INBListS)*2;
    X_C = ITListC./ITListC*2;
%% attenuation correction
    IntensityOutlierThres = quantile([ITAListS;INBList],0.95);
    if IntensityOutlierThres<5500
        IntensityOutlierThres=5500;
    end    
    AttenuationM = max([alleleAS(ITAListS<IntensityOutlierThres);alleleBS(ITBListS<IntensityOutlierThres)]);
    if AttenuationM>8
        AttenuationM = 8;
    end
    alleleAS = alleleAS./(1-(alleleAS-1)./AttenuationM); % allele A copy number signal
    alleleBS = alleleBS./(1-(alleleBS-1)./AttenuationM); % allele B copy number signal
    XabS = alleleAS+alleleBS; % total copy number signal
    X_C = X_C./(1-(X_C-1)./AttenuationM);    
%% Outlier Removal
    % for SNP probes
    OutlierIdxS = ITAListS>IntensityOutlierThres | ITBListS>IntensityOutlierThres |...
        INAListS>IntensityOutlierThres | INBListS>IntensityOutlierThres |(INAListS+INBListS+ITAListS+ITBListS)<2 |...    
        alleleAS>7 | alleleBS>7 | isnan(alleleAS) | isnan(alleleBS) | isinf(alleleAS) | isinf(alleleBS) | ...
        XabS>10| XabS<CNLowThres | isnan(XabS) | isinf(XabS);
    % for CN probes
    X_C = ITListC./INListC*2;
    OutlierIdxC = ITListC>IntensityOutlierThres | INListC>IntensityOutlierThres |...
        INListC<0 | ITListC<0 | X_C>10 | X_C<CNLowThres | isnan(X_C) | isinf(X_C);

    Xab = zeros(MarkerSize,1);
    Xab(MarkerTypeList) = XabS;
    Xab(~MarkerTypeList) = X_C;
    OutlierIdx = true(MarkerSize,1);
    OutlierIdx(MarkerTypeList) = OutlierIdxS;
    OutlierIdx(~MarkerTypeList) = OutlierIdxC;
    
    alleleA = zeros(MarkerSize,1);
    alleleB = zeros(MarkerSize,1);
    alleleA(MarkerTypeList) = alleleAS;
    alleleB(MarkerTypeList) = alleleBS;  
    
    NGenoABList = false(MarkerSize,1);
    NGenoABList(MarkerTypeList) = N_ABGenotype;
    NGenoABList = NGenoABList & ~OutlierIdx;    
    NGenoAABBList = false(MarkerSize,1);
    NGenoAABBList(MarkerTypeList) = N_AABBGenotype;
    NGenoAABBList = NGenoAABBList & ~OutlierIdx;
end 

%% Segmentation
if NeedSegmentation == 1
    % read segmentation data from original BACOM's result
    segFin = fopen([ori_BACOM_result_dir,SampleName,'.CEL','_outputdetectionResult.csv']);
    segMtx = textscan(segFin, '%s %u %u %f', 'delimiter', ',');
    fclose(segFin);
    SegChrList = segMtx{1};
    SegHeadList = segMtx{2};
    SegTailList = segMtx{3};
    SegChrIdx = zeros(size(SegChrList,1),1);
    
    SegSize = 0;
    for SIdx = 1:size(SegChrList,1);
        CIdx = str2double(SegChrList{SIdx});
        if (CIdx<=22)
            SegChrIdx(SIdx) = CIdx;
            SegHeadList(SIdx) = SegHeadList(SIdx)+ChrLoc(CIdx)+1;
            SegTailList(SIdx) = SegTailList(SIdx)+ChrLoc(CIdx)+1;
        else
            SegSize = SIdx-1;
            break;
        end
    end    
    SegChrIdx = SegChrIdx(1:SegSize);
    SegHeadList = SegHeadList(1:SegSize);
    SegTailList = SegTailList(1:SegSize);

    SegLength = zeros(SegSize,1);
    SegABLength = zeros(SegSize,1);
    SegLOHLength = zeros(SegSize,1);
    SegMean = zeros(SegSize,1);
    SegCorruptIdx = false(SegSize,1);
    ProbeRhoList = -1.1*ones(MarkerSize,1);

    for SIdx = 1:SegSize 
        SegAlleleA = alleleA(SegHeadList(SIdx):SegTailList(SIdx));
        SegAlleleB = alleleB(SegHeadList(SIdx):SegTailList(SIdx));
        SegXab = Xab(SegHeadList(SIdx):SegTailList(SIdx));
        SegOutlierIdx = OutlierIdx(SegHeadList(SIdx):SegTailList(SIdx));
        SegXabMiu = mean(SegXab(~SegOutlierIdx)); % mean
        SegXabStd = std(SegXab(~SegOutlierIdx)); % std

        SegNormRobustCNIdx = SegXab<SegXabMiu+RobustStdMultiplier*SegXabStd & SegXab>SegXabMiu-RobustStdMultiplier*SegXabStd; % outliers in the segment
        OutlierIdx(SegHeadList(SIdx):SegTailList(SIdx)) = OutlierIdx(SegHeadList(SIdx):SegTailList(SIdx)) | ~SegNormRobustCNIdx;
        SegRobustCNIdx = ~OutlierIdx(SegHeadList(SIdx):SegTailList(SIdx));
        NGenoABList(SegHeadList(SIdx):SegTailList(SIdx)) = NGenoABList(SegHeadList(SIdx):SegTailList(SIdx)) & SegRobustCNIdx;    
        NGenoAABBList(SegHeadList(SIdx):SegTailList(SIdx)) = NGenoAABBList(SegHeadList(SIdx):SegTailList(SIdx)) & SegRobustCNIdx;    
        SegGenotypeABIdx = NGenoABList(SegHeadList(SIdx):SegTailList(SIdx));
        SegLength(SIdx) = sum(SegRobustCNIdx);       
        SegABLength(SIdx) = sum(SegGenotypeABIdx);

        % calculate allelic correlation coef. (rho) in sliding window        
        SegGenoAB_alleleA = SegAlleleA(SegGenotypeABIdx);
        SegGenoAB_alleleB = SegAlleleB(SegGenotypeABIdx);
        SegCNABRhoList = -ones(SegABLength(SIdx),1);
        SegCNRhoList = -ones(SegTailList(SIdx)-SegHeadList(SIdx)+1,1);
        if SegABLength(SIdx)>=SlideWindowSize
            SegCNABRhoList(1:SlideWindowHalfSize) = corr(SegGenoAB_alleleA(1:SlideWindowSize),SegGenoAB_alleleB(1:SlideWindowSize));
            SegCNABRhoList((end-SlideWindowHalfSize+1):end) = corr( SegGenoAB_alleleA((end-SlideWindowSize+1):end),SegGenoAB_alleleB((end-SlideWindowSize+1):end));            
            for PIdx = (SlideWindowHalfSize+1):(SegABLength(SIdx)-SlideWindowHalfSize)
                SegCNABRhoList(PIdx) = corr( SegGenoAB_alleleA((PIdx-SlideWindowHalfSize):(PIdx+SlideWindowHalfSize)),SegGenoAB_alleleB((PIdx-SlideWindowHalfSize):(PIdx+SlideWindowHalfSize)) );
            end
        end
        SegCNRhoList(SegGenotypeABIdx) = SegCNABRhoList;
        ProbeRhoList(SegHeadList(SIdx):SegTailList(SIdx)) = SegCNRhoList;

        if SegABLength(SIdx)>0
            SegMean(SIdx) = sum(SegAlleleA(SegGenotypeABIdx)+SegAlleleB(SegGenotypeABIdx))/ SegABLength(SIdx);
        end
        if SegABLength(SIdx)<SlideWindowSize        
            SegCorruptIdx(SIdx) = true;
        end
    end
    SegCorruptIdx = SegCorruptIdx | isnan(SegMean) | SegMean>12;
    
    % low correlations
    RhoEdges = -1.01:0.01:1;
    RhoHist = hist(ProbeRhoList(NGenoABList),RhoEdges);
    RhoHist(1) = 0; % remove cnt of rho<-1
    RhoHist = conv(RhoHist,ones(1,10)/10,'same');    
    RhoHistDiff = [RhoHist(3:end),0,0]-RhoHist;
    HistIdx = find(RhoHistDiff>200,1,'last');
    RhoHistDiff = RhoHistDiff(1:HistIdx-5);
    LOH_RhoThresIdx = find(RhoHistDiff<=0,1,'last');    
    LOH_RhoThres = RhoEdges(LOH_RhoThresIdx);    
    if isempty(LOH_RhoThres)
        LOH_RhoThres = -0.2;
    elseif LOH_RhoThres<-0.3
        LOH_RhoThres = -0.3;
    elseif LOH_RhoThres>0.3
        LOH_RhoThres = 0.3;
    end
    
    % LOH
    GenoABNonLOHList = NGenoABList & ProbeRhoList>=LOH_RhoThres;    
    for SIdx = 1:SegSize 
        SegLOHLength(SIdx) = sum(NGenoABList(SegHeadList(SIdx):SegTailList(SIdx)) & ProbeRhoList(SegHeadList(SIdx):SegTailList(SIdx))<LOH_RhoThres);
    end    
    SegLOHPortion = SegLOHLength./SegABLength;
    if sum(1-SegLOHPortion)<0.1 % suspicious low quality sample
        fprintf(fCompOut,'%s\tNF\tNF\n',SampleName);
        continue;
    end

end

%% Copy-neutral component identification
SegHistEdges = 0:0.01:10;
SegHistEdgeSize = size(SegHistEdges,2);
SegHistAllCount = zeros(SegHistEdgeSize,1);
SegHistABCount = zeros(SegHistEdgeSize,1);

for SIdx = 1:SegSize
    if ~SegCorruptIdx(SIdx) && SegLength(SIdx)>LongSegThres
        HistIdx = floor(SegMean(SIdx)/0.01)+1;
        if HistIdx>SegHistEdgeSize
            HistIdx = SegHistEdgeSize;
        end
        SegHistAllCount(HistIdx) = SegHistAllCount(HistIdx)+SegLength(SIdx);
        SegHistABCount(HistIdx) = SegHistABCount(HistIdx)+SegABLength(SIdx)-SegLOHLength(SIdx);        
    end
end

SegHistAllCount = conv(SegHistAllCount,SegConvMask,'same');
SegHistABCount = conv(SegHistABCount,SegConvMask,'same');

SegHistABCountDiff = [zeros(10,1);SegHistABCount(1:end-10)]-SegHistABCount;
DiffHistIdx = find(SegHistABCountDiff>CN_DiffThres,1,'first');
if isempty(DiffHistIdx)
    DiffHistIdx = 1000;
end
SegHistABCount_1stPeak = SegHistABCount(1:DiffHistIdx);
HistIdx = find(SegHistABCount_1stPeak==max(SegHistABCount_1stPeak),1,'first');
Baseline = SegHistEdges(HistIdx); % global normalization baseline

HistIdx = find(SegHistAllCount>0,1,'last');
if isempty(HistIdx) || HistIdx<2 % damaged sample
    fprintf(fCompOut,'%s\tNF\tNF\n',SampleName);
    continue;
end


%% Global normalization
CorrAlelleA = alleleA*2/Baseline; % normalized allele A copy number signal
CorrAlelleB = alleleB*2/Baseline;
CorrSegMean = SegMean*2/Baseline; % normalized segment mean

%% Local normalization

SegNeutralIdx = ~SegCorruptIdx &(abs(CorrSegMean-2)<=LocNormalizeThres & SegLength>LongSegThres);
LocCorrSegMean = zeros(SegSize,1);
for CIdx = 1:22
    LocalChrIdx = (SegChrIdx==CIdx);
    LocGenChrIdx = (MarkerChrList==CIdx);
    LocalNeuIdx = LocalChrIdx & SegNeutralIdx;
    NeuLength = sum(SegLength.* LocalNeuIdx); 
    if NeuLength>NeuLengthThres
        BaselineLoc = sum(CorrSegMean.*SegLength.*LocalNeuIdx)/NeuLength;
        LocCorrSegMean(LocalChrIdx) = CorrSegMean(LocalChrIdx)*2/BaselineLoc;
        CorrAlelleA(LocGenChrIdx) = CorrAlelleA(LocGenChrIdx)*2/BaselineLoc;
        CorrAlelleB(LocGenChrIdx) = CorrAlelleB(LocGenChrIdx)*2/BaselineLoc;
    else
        LocCorrSegMean(LocalChrIdx) = CorrSegMean(LocalChrIdx);
    end                
end

% copy-neutral probes
NeutralNonLOHABProbeIdx = false(MarkerSize,1);
for SIdx = 1:SegSize
    if SegNeutralIdx(SIdx)
        NeutralNonLOHABProbeIdx(SegHeadList(SIdx):SegTailList(SIdx)) = GenoABNonLOHList(SegHeadList(SIdx):SegTailList(SIdx));
    end
end
% global allelic correlation coefficient 
GlobalRho = corr(CorrAlelleA(NeutralNonLOHABProbeIdx),CorrAlelleB(NeutralNonLOHABProbeIdx),'type','Pearson');

%% deletion type classification
SegDelIdx = ~SegCorruptIdx & LocCorrSegMean<DelThres & SegABLength>LongABThres;
% deletion segments
DelSegSize  = sum(SegDelIdx);
DelSegPChiX2 = zeros(DelSegSize,1);
DelSegPNChiX2 = zeros(DelSegSize,1);
DelSegY = zeros(DelSegSize,1);
DelSegLength = zeros(DelSegSize,1);
DelSegABLength = zeros(DelSegSize,1);
DelSegLamda = zeros(DelSegSize,1);
DelSegMiu = zeros(DelSegSize,1);
DelSegVarXab =  zeros(DelSegSize,1);
DelSegRho = zeros(DelSegSize,1);
DelSegChr = zeros(DelSegSize,1);
DelHemiIdx = false(DelSegSize,1);
DelHomoIdx = false(DelSegSize,1);
DelNaNIdx  = false(DelSegSize,1);
DIdx = 0;

for SIdx = 1:SegSize
    SegGenotypeABIdx = NGenoABList(SegHeadList(SIdx):SegTailList(SIdx));    
    SegAlleleA = CorrAlelleA(SegHeadList(SIdx):SegTailList(SIdx));
    SegAlleleB = CorrAlelleB(SegHeadList(SIdx):SegTailList(SIdx));      
    SegAlleleA_AB = SegAlleleA(SegGenotypeABIdx); 
    SegAlleleB_AB = SegAlleleB(SegGenotypeABIdx);
    SegXab = SegAlleleA_AB+SegAlleleB_AB; 
    N = size(SegAlleleA_AB,1);

    if SegDelIdx(SIdx)
        DIdx = DIdx+1;                                  
        DelSegChr(DIdx) = SegChrIdx(SIdx);        
        DelSegLength(DIdx) = SegLength(SIdx);
        DelSegABLength(DIdx) = N;
        DelSegVarXab(DIdx) = var(SegXab);        
        DelSegMiu(DIdx) = mean(SegXab);
        SegASubBX2 = (SegAlleleA_AB-SegAlleleB_AB).^2; 
        DelSegY(DIdx) = sum(SegASubBX2)/(DelSegVarXab(DIdx)*(1-GlobalRho)/(1+GlobalRho));
        DelSegLamda(DIdx) = N*(2-DelSegMiu(DIdx))^2/(DelSegVarXab(DIdx)*(1-GlobalRho)/(1+GlobalRho));         
        DelSegPChiX2(DIdx) = chi2pdf(DelSegY(DIdx),N);
        DelSegPNChiX2(DIdx) = ncx2pdf(DelSegY(DIdx),N,DelSegLamda(DIdx)); 
        
        if(DelSegPChiX2(DIdx)+DelSegPNChiX2(DIdx)== 0)
            DelNaNIdx(DIdx) = true;
        elseif DelSegY(DIdx)<= N
            DelHomoIdx(DIdx) = true;
        elseif DelSegY(DIdx)>= N+DelSegLamda(DIdx)
            DelHemiIdx(DIdx) = true;
        elseif DelSegPChiX2(DIdx)>=DelSegPNChiX2(DIdx)
            DelHomoIdx(DIdx) = true;
        else 
            DelHemiIdx(DIdx) = true;
        end
    end
end

% use the mean of all hemi-deletion segments as the upper-bound for homo-deletion
if any(DelHemiIdx)
    SegMeanHemi = sum(DelSegMiu.*DelSegABLength.*double(DelHemiIdx))./sum(DelSegABLength.*double(DelHemiIdx));
else
    SegMeanHemi = 2;
end
DelHomoIdx = DelHomoIdx & DelSegMiu<=SegMeanHemi;

% normal tissue fraction (alpha) estimation
DelClassifiedIdx = DelHomoIdx | DelHemiIdx;
SegAlpha = zeros(DelSegSize,1);
SegAlpha(DelHemiIdx) = DelSegMiu(DelHemiIdx)-1;
SegAlpha(DelHomoIdx) = DelSegMiu(DelHomoIdx)/2;
Alpha = quantile(SegAlpha(DelClassifiedIdx),0.09);
if isnan(Alpha) || isinf(Alpha)
    Alpha = 0;
end

% purity-corrected segment means
PurityCorrSegMean = (LocCorrSegMean-2*Alpha)/(1-Alpha);
% average copy number/average ploidy
AvgPloidy = sum(PurityCorrSegMean.*SegLength)/sum(SegLength);

% output
fprintf(fCompOut,'%s\t%s\t%s\n',SampleName,num2str(1-Alpha,'%1.3f'),num2str(AvgPloidy,'%1.3f'));
end

fclose(fCompOut);