function [] = Zonation_signature_genes(path_BaharHalpern_2017,genes_path, path_seq_data, output_path)

%%%%%%% choosing landmark genes used to reconstruct zonation profiles %%%%%%%%%%%%%
NUM_ZONES = 8; %number of zonation layers

load(path_BaharHalpern_2017); % zonation gene expression data from Bahar Halpern & Shenhav et al., 2017
nature_genes=all_genes;
fuzzy = Fuzzy.Mat;
qVal = KW.qval;
se = Fuzzy.SE_bootstrap;
seq_data_BaharHalpern_2017 = seq_data;
clear all_genes Fuzzy KW seq_data

% Take expressed, zonated genes with high dynamic range
indin = intersect(find(max(fuzzy,[],2) >= 0.0002 & qVal < 0.2), find( log2(max(fuzzy,[],2)) - log2(min(fuzzy,[],2)) >= log2(1.6) ) );
zon_genes=nature_genes(indin);


%%% finding genes with low inter-mouse/inter-timepoint variablility
load(genes_path);
all_genes=gene_names;
ind=find(ismember(lower(all_genes), lower(zon_genes))); %find genes detected in both Bahar Halpern & Shenhav et al., 2017 and present study.

mus={'ZT00A','ZT00B','ZT06A','ZT06B','ZT12A','ZT12B','ZT18A','ZT18B','ZT00C','ZT12C'};
mean_exp=NaN*ones(length(all_genes),10);
all_norm_data=[];
group=[];
for i=1:length(mus)
    load([path_seq_data mus{i} '.mat'])
    %%%%%%% normalise by non-mitochondrial and non-mup genes %%%%%%%%%%%%%
    non_mito_gene_ind=find(cellfun(@isempty,regexp(all_genes,'\<mt-'))); 
    non_mup_gene_ind=find(cellfun(@isempty,regexp(all_genes,'\<mup')));
    mat_norm =seq_data ./ repmat(sum(seq_data(intersect(non_mito_gene_ind,non_mup_gene_ind),:)),length(all_genes),1); %normalise data by summed expression per cell
    norm_genes = mat_norm ./ repmat(max(mat_norm,[],2),1,size(mat_norm,2));
    means=mean(mat_norm,2);
    mean_exp(:,i)=means;
    all_norm_data=[all_norm_data mat_norm];
    group=[group repmat(mus(i),1,size(mat_norm,2))];
end

%calculate the norm. mean expr over all mice and the noise between all mice
MoM=mean(mean_exp,2); %mean of means
norm_mean=mean_exp./min(mean_exp,[],2);
noise=(std(mean_exp,0,2)./MoM).^2;

%filter genes based on max expression in any mouse < 2.5 mean over all mice
ind3=intersect(find(max(norm_mean,[],2) < 2.5),ind);

indin=find(ismember(lower(nature_genes), lower(all_genes(ind3))));

Mat = fuzzy(indin,:); SE = se(indin,:);
com = (1:9) * Mat'; com = com ./ sum(Mat,2)';

% select the gene sets by max layer and center of mass 
[~,maxLayer] = max(Mat,[],2);
genes_cv = nature_genes(indin(com <= 4.5 & (maxLayer'==1)));
genes_pn = nature_genes(indin(com >= 5 & (maxLayer'>=9)));

display([num2str(length(genes_cv)) ' central landmark genes.'])
display([num2str(length(genes_pn)) ' portal landmark genes.'])


clearvars -except genes_cv genes_pn NUM_ZONES all_genes P Fuzzy output_path seq_data_BaharHalpern_2017 nature_genes

% --- load the hepatocytes from Bahar Halpern & Shenhav et al., 2017 and compute normalised zonation value (CLorig) ---
mat_norm = (seq_data_BaharHalpern_2017 ./ repmat(sum(seq_data_BaharHalpern_2017),length(nature_genes),1)); %normalise expr.matrix by sum of cell
norm_genes = mat_norm ./ repmat(max(mat_norm,[],2),1,size(mat_norm,2)); %norm. genes by max expr. of each gene
cv_vec = sum(norm_genes(ismember(lower(nature_genes),lower(genes_cv)),:)); %sum of central genes 
pn_vec = sum(norm_genes(ismember(lower(nature_genes),lower(genes_pn)),:)); %sum of portal genes

xx = pn_vec./(cv_vec+pn_vec); 
CLorig = (xx - min(xx)) ./ (max(xx) - min(xx));

% --- estimate the Gamma parameters for each layer (merge layers 9 and 8) ---
NUM_ZONES = 8;
Gamma_params = nan(NUM_ZONES,2);
P(P(:,10) == 9,10) = 8; % override layer 9
for i=1:NUM_ZONES
    cells_in_layer = (P(:,10) == i);
    Gamma_params(i,:) = gamfit(CLorig(cells_in_layer));
end

save(output_path, 'Gamma_params','NUM_ZONES', 'genes_cv', 'genes_pn')

end

