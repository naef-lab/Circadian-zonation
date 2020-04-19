%author: El Kholtei, Jakob
%2020


% this script was used for the spatial reconstruction of liver zonation profiles from
% scRNA-seq data in Colas, El Kholtei , Bahar Halpern et al., 2020

% UMI tables can be obtained from GEO with accession code GSE145197
% this script and expression data from Bahar Halpern & Shenhav et al., 2017
% can be obtained from https://c4science.ch/diffusion/10261/ .

% To use this script, download the UMI tables from the GEO and save the
% folder in the downloaded repo. Then unpack the GEO folder.

% Matlab R2018b was used to run the script.



%% 

dirname='./';
addpath([dirname 'matlab functions']);

%% load expression data and convert UMI tables into matlab data files

files=dir([dirname 'GSE145197_RAW/' '*.txt']);

mkdir([dirname 'exp_data/']);
for i=1:length(files)
    clearvars -except i files dirname
    A = importdata([files(i).folder '/' files(i).name], '	',1);
    gene_names = A.textdata(2:end,1);
    cell_barcodes = transpose(A.textdata(1,2:end));
    seq_data = A.data;
    temp = strsplit(files(i).name,'_');
    ZT = temp{end}(1:5);
    filename = strcat([dirname 'exp_data/' ZT '.mat']);
    save(filename, 'cell_barcodes', 'gene_names', 'seq_data', 'ZT')
    if i==1
        save(strcat([dirname 'all_genes.mat']))
    end
end

%% find zonation landmark genes

exp_data_path = [dirname 'exp_data/'];
genes_path = [dirname 'all_genes.mat'];
path_BaharHalpern_2017 = [dirname  'Datasets/BaharHalpern_Shenhav_2017/BaharHalpern_Shenhav_2017.mat']; % path to the file with zonation gene expression data from Bahar Halpern & Shenhav et al., 2017
output_path = [dirname 'landmark_genes.mat'];

Zonation_signature_genes(path_BaharHalpern_2017, genes_path, exp_data_path, output_path)


%% reconstruct zonation profiles for each mouse separately

ZT={'ZT00A','ZT00B','ZT06A','ZT06B','ZT12A','ZT12B','ZT18A','ZT18B','ZT00C','ZT12C'};
landmark_genes_path=[dirname 'landmark_genes.mat'];
mkdir([dirname 'reconstructed_profiles/']);
output_path=[dirname 'reconstructed_profiles/'];
for i=1:length(ZT)
    ZONATION_reconstruction(ZT{i},landmark_genes_path,exp_data_path,output_path)
end
