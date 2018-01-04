function [As_across_bins_pr, As_across_bins_index_pr]=pruning_across_bins(As_across_bins,As_across_bins_index,nneu,criteria,th,style) 
% All pruning techniques compare assemblies with different temporal precisions (disregarding information on the bin size). 
% Two pruning criteria are implemented:
% criteria := 'biggest'  keeps assemblies with biggest dimension (I discard all assemblies whose elements are all present in a bigger
% assembly)
% criteria := 'distance' it groups assemblies according to their similarities in
% terms of constituting elements (all assemblies with a distance smaller
% than th are clustered). th can assume values inside the (0,1) interval (~0 -> no clustering; ~1 -> all in the same cluster). In each group only one assembly is kept according
% to the following criteria: 
% style := 'pvalue' - the assembly with lowest p-value is kept. 
% style := 'occ' - the assembly with most occurrences is kept.
%
%
%
%  © 2016 Russo, Durstewitz.
%  for information please contact eleonora.russo@zi-mannheim.de; daniel.durstewitz@zi-mannheim.de.
%
%  last update 11/01/2016


As_across_bins_pr=[];
As_across_bins_index_pr=[];

if nargin<4, fprintf('\n please, specify a pruning criterion\n'); return; end
if nargin>=4
    if ~strcmp(criteria,'distance') && ~strcmp(criteria,'biggest') 
        fprintf('\n ''biggest'', ''distance'' are the only acceptable criterion specifications\n'); return;
    end   
    if strcmp(criteria,'biggest')
        if nargin>4, fprintf('\n when criteria=''biggest'' is chosen no need for further specifications\n'); return; end               
    end    
    if strcmp(criteria,'distance')
        if nargin<6, fprintf('\n please, specify a style\n'); return; end        
        if nargin<5, fprintf('\n please, specify a cutoff (0 to 1, extremes excluded) and a style\n'); return; end 
        if ~strcmp(style,'pvalue') && ~strcmp(style,'occ') 
            fprintf('\n ''biggest'', ''distance'' are the only acceptable criterion specifications\n')
            return
        end        
        if nargin>5 && (th<=0 || th>=1), fprintf('\n cutoff can assume values inside (0,1) interval (~0 -> no clustering; ~1 -> all assemblies in the same cluster)\n'); return; end 
    end

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nA=length(As_across_bins);

    bin_vector=nan(length(As_across_bins));
    ATstructure=zeros(nneu,length(As_across_bins));
    for i=1:length(As_across_bins)
        ATstructure(As_across_bins{i}.elements,i)=1;
        bin_vector(i)=As_across_bins{i}.bin;
    end
    
%% pruning

switch criteria 
    case 'distance'
    % group assemblies on the base of the similarities between unit
    % composition. All configurations with a distance smaller than th will
    % be associated in the same group. Distances are computed in a N
    % dimensional space, with N being the number of all elements involved
    % in the two compared assemblies (i.e. without taking into account all 
    % the units not involved in those specific assemblies).
    


    D = nan(nA);
    for i=1:nA
        for j=1:nA
            a1=[ATstructure(:,i), ATstructure(:,j)];
            a1(sum(a1,2)==0,:)=[];
            D(i,j) = pdist(a1','cosine');
        end
        D(i,i)=0;
    end

    Z = squareform(D);
    Q2=linkage(Z,'average');
    [~,~,perm2]=dendrogram(Q2,0);

    C = cluster(Q2,'cutoff',th,'criterion','distance','depth',100);
    pat=cell(length(unique(C)),1);
    for i=1:length(unique(C)), 
        pat{i}=find(C==i); 
    end;
    L=cellfun(@length,pat);
    pat_multipli=pat(L>1);

    pat_to_remove=[];
    for i=1:length(pat_multipli) 
        pr=[];
        Nocc=[];
        aus=pat_multipli{i};
        for k=1:length(pat_multipli{i})
            A=aus(k);
            pr(k)=As_across_bins{A}.pr(end);
            Nocc(k)=As_across_bins{A}.Noccurrences(end);
        end
        if strcmp(style,'pvalue'), [~, b]=min(pr); end
        if strcmp(style,'occ'), [~, b]=max(Nocc); end        
        aus(b)=[];
        pat_to_remove=[pat_to_remove; aus];   
    end;
    As_across_bins_index_pr=As_across_bins_index;
    As_across_bins_pr=As_across_bins;
    As_across_bins_index_pr(pat_to_remove)=[];
    As_across_bins_pr(pat_to_remove)=[];

    
    
    case 'biggest'
    % irrespectively to the bin width at which the assembly is detected I
    % keep only assemblies with maximal dimension. Therefore if elements of
    % assembly A are ALL present in a bigger assembly B, I discard A.     
    
    
        [~, b]=sort(sum(ATstructure,1),'descend');
        As_across_bins_sorted=As_across_bins(b);
        As_across_bins_sorted_index=As_across_bins_index(b);
        bin_vector_sorted=bin_vector(b);
        
        to_remove=[];
        for i=1:nA
            test_elementsA=As_across_bins_sorted{i}.elements;
            for j=i+1:nA
                test_elementsB=As_across_bins_sorted{j}.elements;
                C = intersect(test_elementsA,test_elementsB);
                if length(C)==length(test_elementsA)          
                    if As_across_bins_sorted{j}.pr>As_across_bins_sorted{i}.pr   
                        to_remove=[to_remove, j];
                    else
                        to_remove=[to_remove, i];  
                    end
                elseif length(C)==length(test_elementsB)
                    to_remove=[to_remove, j];                
                end 
            end
        end       

        As_across_bins_sorted(unique(to_remove))=[];
        As_across_bins_sorted_index(unique(to_remove))=[];
        bin_vector_sorted(unique(to_remove))=[];
        
        [~, b]=sort(bin_vector_sorted);
        
        As_across_bins_index_pr=As_across_bins_sorted_index(b);
        As_across_bins_pr=As_across_bins_sorted(b);
  
    
    
end




















