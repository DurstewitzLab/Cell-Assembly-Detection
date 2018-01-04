function [assembly]=Main_assemblies_detection(spM,MaxLags,BinSizes,ref_lag,alph,No_th,O_th,bytelimit)
% this function returns cell assemblies detected in spM spike matrix binned 
% at a temporal resolution specified in 'BinSizes' vector and testing for all 
% lags between '-MaxLags(i)' and 'MaxLags(i)'
%
% USAGE: [assembly]=Main_assemblies_detection(spM, MaxLags, BinSizes, ref_lag, alph, Dc, No_th, O_th, bytelimit)
%
% ARGUMENTS:
% spM     := matrix with population spike trains; each row is the spike train (time stamps, not binned) relative to a unit. 
% BinSizes:= vector of bin sizes to be tested;
% MaxLags:= vector of maximal lags to be tested. For a binning dimension of BinSizes(i) the program will test all pairs configurations with a time shift between -MaxLags(i) and MaxLags(i);
% (optional) ref_lag      := reference lag. Default value 2
% (optional) alph      := alpha level. Default value 0.05
% (optional) No_th     := minimal number of occurrences required for an assembly (all assemblies, even if significant, with fewer occurrences than No_th are discarded). Default value 0.
% (optional) O_th      := maximal assembly order (the algorithm will return assemblies of composed by maximum O_th elements).
% (optional) bytelimit := maximal size (in bytes) allocated for all assembly structures detected with a bin dimension. When the size limit is reached the algorithm stops adding new units. 
%
% RETURNS:
% assembly - structure containing assembly information:
%     assembly.parameters       - parameters used to run Main_assemblies_detection
%     assembly.bin{i} contains information about assemblies detected with 
%                     'BinSizes(i)' bin size tested for all lags between 
%                     '-MaxLags(i)' and 'MaxLags(i)'
% 
%        assembly.bin{i}.bin_edges - bin edges (common to all assemblies in assembly.bin{i})
%        assembly.bin{i}.n{j} information about the j-th assembly detected with BinSizes(i) bin size 
%                   elements: vector of units taking part to the assembly (unit order correspond to the agglomeration order)
%                        lag: vector of time lags. '.lag(z)' is the activation delay between .elements(1) and .elements(z+1)
%                         pr: vector of pvalues. '.pr(z)' is the pvalue of the statistical test between performed adding .elements(z+1) to the structure .elements(1:z)
%                       Time: assembly activation time. It reports how many times the complete assembly activates in that bin. .Time always refers to the activation of the first listed assembly element (.elements(1)), that doesn't necessarily corresponds to the first unit firing.
%               Noccurrences: number of assembly occurrence. '.Noccurrences(z)' is the occurrence number of the structure composed by the units .elements(1:z+1) 
%
%
%
%  © 2016 Russo, Durstewitz.
%  for information please contact eleonora.russo@zi-mannheim.de; daniel.durstewitz@zi-mannheim.de.
%
%  last update 11/01/2016
%  last update 13/02/2016 TestPair.m update

if nargin<4 || isempty(ref_lag), ref_lag=2; end;  
if nargin<5 || isempty(alph), alph=0.05; end;  
if nargin<6 || isempty(No_th), No_th=0; end;      % no limitation on the number of assembly occurrences
if nargin<7 || isempty(O_th), O_th=Inf; end;     % no limitation on the assembly order (=number of elements in the assembly)
if nargin<8 || isempty(bytelimit), bytelimit=Inf; end;     % no limitation on assembly dimension

nneu=size(spM,1); % number of units
assemblybin=cell(1,length(BinSizes));  
Dc=100; %length (in # bins) of the segments in which the spike train is divided to compute #abba variance (parameter k).

parfor gg=1:length(BinSizes)


    int=BinSizes(gg);
    maxlag=MaxLags(gg);
    fprintf('%d - testing: bin size=%f sec; max tested lag=%d \n', gg, int, maxlag);    
    %%%%%%%%%%%%%%%%%%%%% Binning  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the spike train is binned with temporal resolution 'int'
    tb=min(spM(:)):int:max(spM(:));  
    binM=zeros(nneu,length(tb)-1,'uint8');
    
    for n=1:nneu
        [ binM(n,:),~] = histcounts(spM(n,:),tb);
    end 
    
    if size(binM,2)-MaxLags(gg)<100
        fprintf('Warning: testing bin size=%f. The time series is too short, consider taking a longer portion of spike train or diminish the bin size to be tested \n', int);
    else
        %%%%%%%%%%%%%%%%%%%%%%% Analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [assemblybin{gg}]=FindAssemblies_recursive(binM,maxlag,alph,gg,Dc,No_th,O_th,bytelimit,ref_lag);  % it returns assemblies detected at specific temporal resolution 'int' 
        if ~isempty(assemblybin{gg})
            assemblybin{gg}.bin_edges=tb;
        end
        fprintf('%d - testing done\n',gg); 
        fname = sprintf('assembly%d.mat', gg);
        parsave(fname,assemblybin{gg})
    end
end
assembly.bin=assemblybin;
assembly.parameters.alph=alph;
assembly.parameters.Dc=Dc;
assembly.parameters.No_th=No_th;
assembly.parameters.O_th=O_th;
assembly.parameters.bytelimit=bytelimit;

fprintf('\n');
end    

function parsave(fname,aus)
 save(fname,'aus')
end    
