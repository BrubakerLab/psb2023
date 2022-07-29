clear
clc
close all force

%-------------------------------------------------------------------------%
% PARAMETERS
%-------------------------------------------------------------------------%
mkdir Plots
pctThresh = 67;                                 % percentile of gene variance to keep for models
LV_init   = 8;                                  % initial LV's for internal PLS CV
k         = 6;                                  % Cross validation folds for PLS-DA model
numLV     = 2;                                  % LV selected after cross-validation
sz1       = 600;                                % size of scatter plot dots for Case 2 and Case 3 Models
%-------------------------------------------------------------------------%
% READ DATA IN
%-------------------------------------------------------------------------%

gse_data  = readtable('GSE16879_UC_pre-treatment_expression.txt','Delimiter','\t','ReadRowNames',1);
hmp_gene  = readtable('ihmp_uc_rna.txt','Delimiter','\t','ReadRowNames',1);
hmp_metb  = readtable('ihmp_uc_metabolomics.txt','Delimiter','\t','ReadRowNames',1);
hmp_micr  = readtable('ihmp_uc_microbiota.txt','Delimiter','\t','ReadRowNames',1);
hmp_phen  = readtable('ihmp_uc_phenotypes.txt','Delimiter','\t','ReadRowNames',1);

%----------------------------------------------------------------------------------------------------------------------%
% PREPROCESSING
%----------------------------------------------------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% PROCESS GSE16879 DATA
%-------------------------------------------------------------------------%

Y_ifx           = gse_data{1,:};                % create infliximab response variable (1 = Responder, -1 = Resistant)
gse_data(1,:)   = [];                           % remove infliximab response row from table
geneVar         = var(gse_data{:,:}');          % calculate variances of genes
varThresh       = prctile(geneVar,pctThresh);   % gene variance threshold to drop
varIdx          = find(geneVar >= varThresh);   % get rows with ihg variance
gse_data        = gse_data(varIdx,:);           % make new gene data with high variance genes

%-------------------------------------------------------------------------%
% PROCESS iHMP Data - Filter for UC patients and High Variance GSE Genes
%-------------------------------------------------------------------------%
ihmp_idx = zeros(height(hmp_metb),1);

for i = 1:height(hmp_metb)
    if isequal(hmp_metb.DiseaseStatus(i),{'UC'})
        ihmp_idx(i) = 1;
    end
end

idx = find(ihmp_idx);
hmp_metb = hmp_metb(idx,:);
hmp_micr = hmp_micr(idx,:);
hmp_gene = hmp_gene(:,idx);
hmp_phen = hmp_phen(idx,:);

% Check that iHMP samples match
if  isequal(hmp_metb.Properties.RowNames,hmp_micr.Properties.RowNames) &...
        isequal(hmp_metb.Properties.RowNames,hmp_phen.Properties.RowNames) &...
        isequal(hmp_metb.Properties.RowNames,hmp_gene.Properties.VariableNames')

    disp('All iHMP Samples Match')

end

% Check that genes match
[~,a,b] = intersect(gse_data.Properties.RowNames,hmp_gene.Properties.RowNames);
gse_data = gse_data(a,:);
hmp_gene = hmp_gene(b,:);
if isequal(gse_data.Properties.RowNames,hmp_gene.Properties.RowNames)
    disp('GSE and iHMP Genes Match')
end

% Create normalized datasets
gse_norm        = zscore(gse_data{:,:}');       % Create normalized dataset
hmp_geneNorm    = zscore(hmp_gene{:,:}');
hmp_metabNorm   = zscore(hmp_metb{:,3:end});
hmp_microNorm   = hmp_micr{:,3:end};

%----------------------------------------------------------------------------------------------------------------------%
% GSE PLS-DA Model Construction, Validation, and Prediction of iHMP Data
%----------------------------------------------------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Construct GSE-PLSDA Model, identify VIP Genes for IFX, Extract LV's
%-------------------------------------------------------------------------%

[~, ~, ~, ~,~,pctvar, mse, ~]  = plsregress(gse_norm,Y_ifx',LV_init,'cv',k);

figure
plot(1:LV_init,cumsum(100*pctvar(2,:)),'-bo','LineWidth',2);
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in Y');
title('GSE16879 Training Data')
saveas(gcf,'Plots/Fig2_trainPLS_pctvar','epsc')

figure
plot(1:LV_init,mse(2,1:end-1),'-bo','LineWidth',2);
xlabel('Number of PLS components');
ylabel('Mean Squared Error in Y');
title('GSE16879 Training Data')
%ylim([0, 1])
saveas(gcf,'Plots/Fig2_trainPLS_mse','epsc')

[XL, YL, XS, ~,beta,pctvar, mse, stats]     = plsregress(gse_norm,Y_ifx',numLV,'cv',k);
gseVip                                      = pls_vip(stats, XS, XL, YL,[1:numLV]);
gseVipTbl                                   = table(gseVip,'RowNames',gse_data.Properties.RowNames);
gseVipsByLV                                 = zeros(size(XL));
for i = 1:2
    gseVipsByLV(:,i) = pls_vip(stats, XS(:,i), XL(:,i), YL,i);
end
gseVipsByLV        = splitvars(table(gseVipsByLV,'RowNames',gse_data.Properties.RowNames));
gseVipTbl          = [gseVipTbl gseVipsByLV];

% Outputs for VIP-IFX Score Calculation

% Plot Top VIPs for Case 3
s = gse_data(:,Y_ifx == 1);
r = gse_data(:,Y_ifx == -1);
gsePlotVIP = [gseVipTbl table(log2(mean(s{:,:}')./mean(r{:,:}'))')];

figure
scatter(gsePlotVIP{:,4},gsePlotVIP{:,1},40,'k','filled')
xlabel('Log2 Fold Change S/R')
ylabel('GSE VIP Score')
yline(2)
xline(0)
saveas(gcf,'Plots/Fig2_vip_volcano','epsc')

writetable(gseVipTbl,'Plots/GSE_PLSDA_IFX_VIP_table.txt','WriteRowNames',1)
[x] = find(gseVipTbl{:,1} > 2);                                                 % 70 gene Infliximab Resistance signature
tData= gse_data(x,:);
S = tData(:,Y_ifx == 1);
R = tData(:,Y_ifx == -1);
FC = log2(mean(S{:,:}')./mean(R{:,:}'));
[y] = find(FC < 0);                                                            
iData = hmp_gene(x,:);
writetable(iData(y,:),'Plots/InfliximabScoreGenes.txt','WriteRowNames',1);

resSig = gsePlotVIP;
pGenes = find(resSig{:,1}>2 & resSig{:,end}<0);

cg = clustergram(zscore(gse_data{x,:}'),'ColorMap','redbluecmap',...
                 'RowLabels',gse_data.Properties.VariableNames,...
                 'ColumnLabels',gse_data.Properties.RowNames(x))
plot(cg)
saveas(gcf,'Plots/Fig2_signatureHeatmap','epsc')

% Visualize PLSR Model Scores on Significant Phenotype Variables
%
figure % Infliximab Response
gscatter(XS(:,1),XS(:,2),Y_ifx,[],'.',80)
title('GSE16879 UC Infliximab Response')
xlabel([{pctvar(2,1) '%VarExp'}])
ylabel([{pctvar(2,2) '%VarExp'}])
saveas(gcf,'Plots/Fig2_trainPLS_scoresPlot','epsc')

%----------------------------------------------------------------------------------------------------------------------%
% Case 1: Supervised Learning, Predicting iHMP Binary Response and Construct Metabolomics/Microbiota Models
%----------------------------------------------------------------------------------------------------------------------%
Y_hmp = sign([ones(size(hmp_geneNorm,1),1) hmp_geneNorm]*beta); % Predicted Infliximab Response - Continuous
XS_hmp = hmp_geneNorm*stats.W;                                  % Projection of iHMP onto GSE-LV's

% Metaboloics Model
[XL_metab1, YL_metab1, XS_metab1, ~,~,pctvar_metab1, ~, stats_metab1]       = plsregress(hmp_metabNorm,Y_hmp,numLV,'cv',k);
metab1_Vips                                                                 = pls_vip(stats_metab1,XS_metab1,XL_metab1,YL_metab1,[1:numLV]);
metab1_VipTbl                                                               = table(metab1_Vips,'RowNames',hmp_metb.Properties.VariableNames(3:end));

figure
gscatter(XS_metab1(:,1),XS_metab1(:,2),Y_hmp,[],'.',80)
title('Case 1: iHMP Samples on Metabolomics LVs')
xlabel([{pctvar_metab1(2,1) '%VarExp'}])
ylabel([{pctvar_metab1(2,2) '%VarExp'}])
saveas(gcf,'Plots/Fig4_case1_metab1_scores','epsc')

% Microbiota Model
[XL_micro1, YL_micro1, XS_micro1, ~,~,pctvar_micro1, ~, stats_micro1]       = plsregress(hmp_microNorm,Y_hmp,numLV,'cv',k);
micro1_Vips                                                                 = pls_vip(stats_micro1,XS_micro1,XL_micro1,YL_micro1,[1:numLV]);
micro1_VipTbl                                                               = table(micro1_Vips,'RowNames',hmp_micr.Properties.VariableNames(3:end));

% Write significant features tables
writetable(metab1_VipTbl,'Plots/Case1_metab1_vipTable.txt','WriteRowNames',1)
m1 = find(metab1_VipTbl{:,1}>1);
writetable(metab1_VipTbl(m1,1),'Plots/Case1_metab1_vipTable_sig.txt','WriteRowNames',1)


writetable(micro1_VipTbl,'Plots/Case1_micro1_vipTable.txt','WriteRowNames',1)
i1 = find(micro1_VipTbl{:,1}>1);
writetable(micro1_VipTbl(i1,1),'Plots/Case1_micro1_vipTable_sig.txt','WriteRowNames',1)

figure
gscatter(XS_micro1(:,1),XS_micro1(:,2),Y_hmp,[],'.',80)
title('Case 1: iHMP Samples on Microbiota LVs')
xlabel([{pctvar_micro1(2,1) '%VarExp'}])
ylabel([{pctvar_micro1(2,2) '%VarExp'}])
saveas(gcf,'Plots/Fig3_case1_micro1_scores','epsc')

%----------------------------------------------------------------------------------------------------------------------%
% Case 2: Predicting 2D Continuous LV Score
%----------------------------------------------------------------------------------------------------------------------%

% Metaboloics Model
[XL_metab2, YL_metab2, XS_metab2, ~,~,pctvar_metab2, ~, stats_metab2]       = plsregress(hmp_metabNorm,XS_hmp,numLV,'cv',k);
metab2_Vips                                                                 = pls_vip(stats_metab2,XS_metab2,XL_metab2,YL_metab2,[1:numLV]);
metab2_VipTbl                                                               = table(metab2_Vips,'RowNames',hmp_metb.Properties.VariableNames(3:end));

figure
scatter(XS_metab2(:,1),XS_metab2(:,2),sz1,XS_hmp(:,1),'filled')
title('Case 2: iHMP Samples on Metabolomics LVs, XS1')
xlabel([{pctvar_metab2(2,1) '%VarExp'}])
ylabel([{pctvar_metab2(2,2) '%VarExp'}])
colorbar
saveas(gcf,'Plots/Fig4_case2_metab2_scores_XS1','epsc')

figure
scatter(XS_metab2(:,1),XS_metab2(:,2),sz1,XS_hmp(:,2),'filled')
title('Case 2: iHMP Samples on Metabolomics LVs, XS2')
xlabel([{pctvar_metab2(2,1) '%VarExp'}])
ylabel([{pctvar_metab2(2,2) '%VarExp'}])
colorbar
saveas(gcf,'Plots/Fig4_case2_metab2_scores_XS2','epsc')

% Microbiota Model
[XL_micro2, YL_micro2, XS_micro2, ~,~,pctvar_micro2, ~, stats_micro2]       = plsregress(hmp_microNorm,XS_hmp,numLV,'cv',k);
micro2_Vips                                                                 = pls_vip(stats_micro2,XS_micro2,XL_micro2,YL_micro2,[1:numLV]);
micro2_VipTbl                                                               = table(micro2_Vips,'RowNames',hmp_micr.Properties.VariableNames(3:end));

figure
scatter(XS_micro2(:,1),XS_micro2(:,2),sz1,XS_hmp(:,1),'filled')
title('Case 2: iHMP Samples on Microbiota LVs XS1')
xlabel([{pctvar_micro2(2,1) '%VarExp'}])
ylabel([{pctvar_micro2(2,2) '%VarExp'}])
colorbar
saveas(gcf,'Plots/Fig3_case2_micro2_scores_XS1','epsc')

figure
scatter(XS_micro2(:,1),XS_micro2(:,2),sz1,XS_hmp(:,2),'filled')
title('Case 2: iHMP Samples on Microbiota LVs XS2')
xlabel([{pctvar_micro2(2,1) '%VarExp'}])
ylabel([{pctvar_micro2(2,2) '%VarExp'}])
colorbar
saveas(gcf,'Plots/Fig3_case2_micro2_scores_XS2','epsc')

% Write significant features tables
writetable(metab2_VipTbl,'Plots/Case2_metab2_vipTable.txt','WriteRowNames',1)
m2 = find(metab2_VipTbl{:,1}>1);
writetable(metab2_VipTbl(m2,1),'Plots/Case2_metab2_vipTable_sig.txt','WriteRowNames',1)

writetable(micro2_VipTbl,'Plots/Case2_micro2_vipTable.txt','WriteRowNames',1)
i2 = find(micro2_VipTbl{:,1}>1);
writetable(micro2_VipTbl(i2,1),'Plots/Case2_micro2_vipTable_sig.txt','WriteRowNames',1)

%----------------------------------------------------------------------------------------------------------------------%
% Case 3: Predicting Infliximab Resistance Gene Score
%----------------------------------------------------------------------------------------------------------------------%
ifxScore = readtable('Case3_IFXgeneSore_ssGSEA.txt','Delimiter',',','ReadRowNames',1);
if isequal(ifxScore.Properties.VariableNames,hmp_gene.Properties.VariableNames)
    disp('Data Names Match for Case 3')
end

% Metaboloics Model
[XL_metab3, YL_metab3, XS_metab3, ~,~,pctvar_metab3, ~, stats_metab3]       = plsregress(hmp_metabNorm,-1*ifxScore{1,:}',numLV,'cv',k);
metab3_Vips                                                                 = pls_vip(stats_metab3,XS_metab3,XL_metab3,YL_metab3,[1:numLV]);
metab3_VipTbl                                                               = table(metab3_Vips,'RowNames',hmp_metb.Properties.VariableNames(3:end));

figure
scatter(XS_metab3(:,1),XS_metab3(:,2),sz1,ifxScore{1,:},'filled')
title('Case 3: iHMP Samples on Metabolomics LVs, IFX Score')
xlabel([{pctvar_metab3(2,1) '%VarExp'}])
ylabel([{pctvar_metab3(2,2) '%VarExp'}])
colorbar
saveas(gcf,'Plots/Fig4_case3_metab3_scores_ifxScore','epsc')

% Microbiota Model
[XL_micro3, YL_micro3, XS_micro3, ~,~,pctvar_micro3, ~, stats_micro3]       = plsregress(hmp_microNorm,ifxScore{1,:}',numLV,'cv',k);
micro3_Vips                                                                 = pls_vip(stats_micro3,XS_micro3,XL_micro3,YL_micro3,[1:numLV]);
micro3_VipTbl                                                               = table(micro3_Vips,'RowNames',hmp_micr.Properties.VariableNames(3:end));

figure
scatter(XS_micro3(:,1),XS_micro3(:,2),sz1,ifxScore{1,:},'filled')
title('Case 3: iHMP Samples on Microbiota LVs IFX Score')
xlabel([{pctvar_micro3(2,1) '%VarExp'}])
ylabel([{pctvar_micro3(2,2) '%VarExp'}])
colorbar
saveas(gcf,'Plots/Fig3_case3_micro3_scores_ifxScore','epsc')

% Write significant features tables
writetable(metab3_VipTbl,'Plots/Case3_metab3_vipTable.txt','WriteRowNames',1)
m3 = find(metab3_VipTbl{:,1}>1);
writetable(metab3_VipTbl(m3,1),'Plots/Case3_metab3_vipTable_sig.txt','WriteRowNames',1)

writetable(micro3_VipTbl,'Plots/Case3_micro3_vipTable.txt','WriteRowNames',1)
i3 = find(micro3_VipTbl{:,1}>1);
writetable(micro3_VipTbl(i3,1),'Plots/Case3_micro3_vipTable_sig.txt','WriteRowNames',1)

%----------------------------------------------------------------------------------------------------------------------%
% Compare Predicted Metabolomics and Microbiota Features Across Models
%----------------------------------------------------------------------------------------------------------------------%

% Case 1 Metabolomics x Case 2 Metabolomics
figure
scatter(metab1_VipTbl{:,1}, metab2_VipTbl{:,1},100,'k','filled')
xlabel('Metabolomics Case 1 VIP')
ylabel('Metabolomics Case 2 VIP')
title('VIP Compare')
yline(1)
xline(1)
saveas(gcf,'Plots/Fig4_case1_vs_case2_metabVIP','epsc')

% Case 1 Metabolomics x Case 3 Metabolomics
figure
scatter(metab1_VipTbl{:,1}, metab3_VipTbl{:,1},100,'k','filled')
xlabel('Metabolomics Case 1 VIP')
ylabel('Metabolomics Case 3 VIP')
title('VIP Compare')
yline(1)
xline(1)
saveas(gcf,'Plots/Fig4_case1_vs_case3_metabVIP','epsc')

% Case 2 Metabolomics x Case 3 Metabolomics
figure
scatter(metab2_VipTbl{:,1}, metab3_VipTbl{:,1},100,'k','filled')
xlabel('Metabolomics Case 2 VIP')
ylabel('Metabolomics Case 3 VIP')
title('VIP Compare')
yline(1)
xline(1)
saveas(gcf,'Plots/Fig4_case2_vs_case3_metabVIP','epsc')

% Case 1 Microbiota x Case 2 Microbiota
figure
scatter(micro1_VipTbl{:,1}, micro2_VipTbl{:,1},100,'k','filled')
xlabel('Microbiota Case 1 VIP')
ylabel('Microbiota Case 2 VIP')
title('VIP Compare')
yline(1)
xline(1)
saveas(gcf,'Plots/Fig3_case1_vs_case2_microVIP','epsc')

% Case 1 Microbiota x Case 3 Microbiota
figure
scatter(micro1_VipTbl{:,1}, micro3_VipTbl{:,1},100,'k','filled')
xlabel('Microbiota Case 1 VIP')
ylabel('Microbiota Case 3 VIP')
title('VIP Compare')
yline(1)
xline(1)
saveas(gcf,'Plots/Fig3_case1_vs_case3_microVIP','epsc')

% Case 2 Microbiota x Case 3 Microbiota
figure
scatter(micro2_VipTbl{:,1}, micro3_VipTbl{:,1},100,'k','filled')
xlabel('Microbiota Case 2 VIP')
ylabel('Microbiota Case 3 VIP')
title('VIP Compare')
yline(1)
xline(1)
saveas(gcf,'Plots/Fig3_case2_vs_case3_microVIP','epsc')

%----------------------------------------------------------------------------------------------------------------------%
% FUNCTIONS
%----------------------------------------------------------------------------------------------------------------------%

function [vipScores] = pls_vip(stats,XS,XL,YL,i)
% Calculate normalized PLS weights
W0 = bsxfun(@rdivide,stats.W(:,i),sqrt(sum(stats.W(:,i).^2,1)));
sumSq = sum(XS.^2,1).*sum(YL.^2,1);
% Calculate VIP scores for NCOMP components
vipScores = sqrt(size(XL,1) * sum(bsxfun(@times,sumSq,W0.^2),2) ./ sum(sumSq,2));
end
