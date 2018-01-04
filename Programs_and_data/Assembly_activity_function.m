function [assembly_activity]=Assembly_activity_function(As_across_bins, assembly, spM, BinSizes, lagChoice, act_count)
% This function returns assembly activity. Two option have to be chosen:
%
% ARGUMENTS:
% As_across_bins := structure containing assembly information output of ASSEMBLIES_ACROSS_BINS;
%     assembly :=  structure containing assembly information output of Main_assemblies_detection.m
%      spM := matrix with population spike trains; each row is the spike train (time stamps, not binned) relative to a unit. 
%      BinSizes := vector of bin sizes to be tested;
% 
%   (optional) lagChoice   -  The user can chose if to consider only the assembly activation (marked at 
% the occurrence of the first assembly unit spike) or the activity of the whole 
% assembly (from the spike of the first assembly unit to the spike of the last assembly unit)
% lagChoice='beginning' ---> only assembly activation 
% lagChoice='duration'  ---> activity stays up during the whole assembly duration 


%   (optional) act_count   - options on how to count activity
% The user can chose if to consider assemblies active when all assembly
% units are active or also if only part of them fire in a coordinated
% fashion:
% act_count='full' --> only full assembly occurrences are considered (when all assembly units are active)
% act_count='partial' --> partial activations are also taken into account
% act_count='combined' --> partial activations are weighted by the fraction of
% active assembly units.


% RETURNS:
% assembly_activity   - assembly activation profile across recording time 
% according to different  possible definitions. 
%
%
%
%  © 2016 Russo, Durstewitz.
%  for information please contact eleonora.russo@zi-mannheim.de; daniel.durstewitz@zi-mannheim.de.
%
%  last update 11/01/2016
%


if nargin<5 || isempty(lagChoice), lagChoice='duration'; end;  
if nargin<6 || isempty(act_count), act_count='full'; end;      


as=As_across_bins;
assembly_activity=cell(length(as),1);
 
for n=1:length(as)
   
   tb=assembly.bin{find(BinSizes==as{n}.bin)}.bin_edges; 
  
   switch act_count
       case 'full'
            [activity]=DetectAssembly(spM ,as{n},tb);  %  only full assembly occurrences are considered
       case 'partial'
            [activity]=DetectAssembly_additive(spM ,as{n},tb);  % partial activations are also took into account
       case 'combined'
            [activity]=DetectAssembly_additive_weighted(spM ,as{n},tb);  % partial activations are also took into account but weighting differently full assembly contributions from partial assembly contributions.  
   end 
   ta=tb(1:end-1)+as{n}.bin/2; %tb are the binning edges, I want the center
   WOO=[ta',activity']; 
     
 
   
   switch lagChoice
       case 'beginning'
           assembly_activity{n}=WOO;
       case 'duration'

           alag=As_across_bins{n}.lag(end);
           activity_lag=activity;
           for j=1:alag
                activity_lag=[activity_lag; [zeros(1,j), activity(1:end-j)]];
           end
           VOO=[ta',sum(activity_lag,1)'];     
 
           assembly_activity{n}=VOO;
   end

end

end






function [activity]=DetectAssembly(spM, Assembly, tb)
% only occurrences of full assemblies are considered. The 'activity' vector
% returns the number of times the complete set of assembly units is firing
% in a bin. 

eelements=Assembly.elements;
alag=Assembly.lag;
Sbin_assembly=zeros(length(eelements),length(tb)-1);
for i=1:length(eelements)
    [Sbin_assembly(i,:),~] = histcounts(spM(eelements(i),:),tb);
end

Sbin_shifted=nan(size(Sbin_assembly));
for e1=1:length(Assembly.elements)
    if alag(e1)==0
        Sbin_shifted(e1,:)=Sbin_assembly(e1,:);
    elseif alag(e1)>0
        Sbin_shifted(e1,1:end-alag(e1))=Sbin_assembly(e1,1+alag(e1):end);  
    else
        Sbin_shifted(e1,1-alag(e1):end)=Sbin_assembly(e1,1:end+alag(e1));
    end
end
aus=sum(Sbin_shifted,1);
Sbin_shifted(:,isnan(aus))=0; % cut at the borders
activity=min(Sbin_shifted,[],1);   
end


function [activationtot]=DetectAssembly_additive(spM, Assembly, tb)

eelements=Assembly.elements;
alag=Assembly.lag;
Sbin_assembly=zeros(length(eelements),length(tb)-1);
for i=1:length(eelements)
    [Sbin_assembly(i,:),~] = histcounts(spM(eelements(i),:),tb);
end

maxrate=max(Sbin_assembly(:)); 
Sbin_p=cell(1,maxrate);
for i=1:maxrate
    Sbin_p{i}=zeros(size(Sbin_assembly));
    Sbin_p{i}(find(Sbin_assembly>=i))=1;  
end

activity=zeros(maxrate,size(Sbin_assembly,2));
Sbin_shifted=nan(size(Sbin_assembly));
for i=1:maxrate
    for e1=1:length(Assembly.elements)
        if alag(e1)==0
            Sbin_shifted(e1,:)=Sbin_p{i}(e1,:);
        elseif alag(e1)>0
            Sbin_shifted(e1,1:end-alag(e1))=Sbin_p{i}(e1,1+alag(e1):end);  
        else
            Sbin_shifted(e1,1-alag(e1):end)=Sbin_p{i}(e1,1:end+alag(e1));
        end
    end
    activity(i,:)=sum(Sbin_shifted,1);   
    activity(i,find(activity(i,:)==1))=0;   
end
    
activationtot=sum(activity,1);   
activationtot(isnan(activationtot))=0;
activationtot=activationtot/length(eelements);
end



function [activationtot]=DetectAssembly_additive_weighted(spM, Assembly, tb)

eelements=Assembly.elements;
alag=Assembly.lag;
Sbin_assembly=zeros(length(eelements),length(tb)-1);
for i=1:length(eelements)
    [Sbin_assembly(i,:),~] = histcounts(spM(eelements(i),:),tb);
end

maxrate=max(Sbin_assembly(:));  
Sbin_p=cell(1,maxrate);
for i=1:maxrate
    Sbin_p{i}=zeros(size(Sbin_assembly));
    Sbin_p{i}(find(Sbin_assembly>=i))=1;  
end


activity=zeros(maxrate,size(Sbin_assembly,2));
Sbin_shifted=nan(size(Sbin_assembly));
for i=1:maxrate
    for e1=1:length(Assembly.elements)
        if alag(e1)==0
            Sbin_shifted(e1,:)=Sbin_p{i}(e1,:);
        elseif alag(e1)>0
            Sbin_shifted(e1,1:end-alag(e1))=Sbin_p{i}(e1,1+alag(e1):end);  
        else
            Sbin_shifted(e1,1-alag(e1):end)=Sbin_p{i}(e1,1:end+alag(e1));
        end
    end
    activity(i,:)=sum(Sbin_shifted,1);   
    activity(i,:)= activity(i,:).^5;
    activity(i,find(activity(i,:)==1))=0;   
end
    
activationtot=sum(activity,1);   
activationtot(isnan(activationtot))=0;
activationtot=activationtot/length(eelements)^5;
end


















