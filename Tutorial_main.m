%%%%%%%%%%%%%%%%%%%%%%%%%% TUTORIAL ON CELL ASSEMBLY DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD DATA
load('test_data.mat');                                  
nneu=size(spM,1);  % nneu is number of recorded units

BinSizes=[0.015 0.025 0.04 0.06 0.085 0.15 0.25 0.4 0.6 0.85 1.5];
MaxLags=[10 10 10 10 10 10 10 10 10 10 10];

% ASSEMBY DETECTION
[assembly]=Main_assemblies_detection(spM,MaxLags,BinSizes);

%% %%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%

% ASSEMBLY REORDERING
[As_across_bins,As_across_bins_index]=assemblies_across_bins(assembly,BinSizes);

display='raw';
% display='clustered';

% VISUALIZATION
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins, nneu, BinSizes, display);

%% %%%%%%%%%%%%%%%%%%%%%%%% PRUNING %%%%%%%%%%%%%%%%%%%%%%%%
clf
% PRUNING: criteria = 'biggest';
criteria = 'biggest';
[As_across_bins_pr,As_across_bins_index_pr]=pruning_across_bins(As_across_bins,As_across_bins_index,nneu,criteria);

display='raw';
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins_pr, nneu, BinSizes, display);

%%
clf

% PRUNING: criteria = 'distance';
criteria = 'distance';
% th=0.7;
th=0.3;

style='pvalue';
% style='occ';

[As_across_bins_pr,As_across_bins_index_pr]=pruning_across_bins(As_across_bins,As_across_bins_index,nneu,criteria,th,style);
display='raw';
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins_pr, nneu, BinSizes, display);


%% %%%%%%%%%%%%%%%%%%%%%%%% ASSEMBLY ACTVATION %%%%%%%%%%%%%%%%%%%%%%%%
clf

criteria = 'biggest';
[As_across_bins_pr,As_across_bins_index_pr]=pruning_across_bins(As_across_bins,As_across_bins_index,nneu,criteria);

lagChoice = 'beginning';
% lagChoice='duration';

act_count = 'full';
[assembly_activity]=Assembly_activity_function(As_across_bins_pr,assembly,spM,BinSizes,lagChoice,act_count);


for i=1:length(assembly_activity)
    subplot(5,1,i)
    plot(assembly_activity{i}(:,1),assembly_activity{i}(:,2));
    hold on
end













