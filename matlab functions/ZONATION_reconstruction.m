function [] = ZONATION_reconstruction(fun_ZT,landmark_genes_path,filtered_seqdata_path,fun_output_path)

%load zonation landmark genes
load(landmark_genes_path)

% LOAD DATA:
nameFigure = fun_ZT;
load([filtered_seqdata_path nameFigure]);

%%%%%%% normalise by non-mitochondrial genes %%%%%%%%%%%%%
non_mito_gene_ind=find(cellfun(@isempty,regexp(gene_names,'\<mt-')));
non_mup_gene_ind=find(cellfun(@isempty,regexp(gene_names,'\<mup')));

mat_norm =seq_data ./ repmat(sum(seq_data(intersect(non_mito_gene_ind,non_mup_gene_ind),:)),length(gene_names),1); %normalise data by summed expression per cell
norm_genes = mat_norm ./ repmat(max(mat_norm,[],2),1,size(mat_norm,2)); %normalise normalised data by max expr. per gene
display('Normalising only by non-mitochondrial and non-mup genes!')

% sum of CV landmark genes 
cv_vec = sum(norm_genes(ismember(lower(gene_names), lower(genes_cv)),:));

% sum of PN landmark genes 
pn_vec = sum(norm_genes(ismember(lower(gene_names),lower(genes_pn)),:));

% CLorig
xx = pn_vec./(cv_vec+pn_vec);
CLorig = (xx - min(xx)) ./ (max(xx) - min(xx));


% =============== calculate the probabilities of each cell and each layer =====================
Pdf_times_prior_Mat = nan(length(CLorig), NUM_ZONES);

for i=1:NUM_ZONES
    pdf = gampdf(CLorig, Gamma_params(i,1),Gamma_params(i,2));
    pdf(pdf == 0) = min(pdf(pdf > 0)); % set 0 value to min pdf
    %     Pdf_times_prior_Mat(:,i) = pdf * P_Z(i); % take into consideration the prior
    Pdf_times_prior_Mat(:,i) = pdf;
end
Pmat = Pdf_times_prior_Mat ./ repmat(sum(Pdf_times_prior_Mat,2), 1, NUM_ZONES); % div by sum columns
Pmat = Pmat ./ repmat(sum(Pmat), size(Pmat,1), 1); % div by sum rows (weights)
MeanGeneExp = mat_norm * Pmat;


% === obtain standard error by bootstrapping ===
n = size(mat_norm,2);
reps = 500;
bootGenes = zeros(size(mat_norm,1), NUM_ZONES, reps); %  (genes x zones x iterations) matrix
for i=1:reps
    samples = randsample(1:n, n, 'true'); % sample n cells with replacement
    bootGenes(:,:,i) = mat_norm(:,samples) * Pmat(samples,:);
end
SE = std(bootGenes,[],3);

% Calc qvals by permutations:

vec = MeanGeneExp ./ repmat(mean(MeanGeneExp,2),1,size(MeanGeneExp,2));
stat_real = max(vec,[],2)-min(vec,[],2);
NUM_PERM=1000;
stat_rand=zeros(size(MeanGeneExp,1),NUM_PERM);
for i=1:NUM_PERM
    mat_normr=mat_norm(:,randperm(size(mat_norm,2)));
    MeanGeneExpr = mat_normr * Pmat;
    vec=MeanGeneExpr./repmat(mean(MeanGeneExpr,2),1,size(MeanGeneExpr,2));
    stat_rand(:,i)=max(vec,[],2)-min(vec,[],2);
end
Z=(stat_real-mean(stat_rand,2))./std(stat_rand,[],2);
pval_perm=1-normcdf(Z);

% Calc q-values
q_vals = mafdr(pval_perm,'BHFDR',true);

clearvars -except SE MeanGeneExp seq_data mat_norm gene_names maxP Pmat fun_ZT q_vals fun_output_path
save([fun_output_path fun_ZT])

end

