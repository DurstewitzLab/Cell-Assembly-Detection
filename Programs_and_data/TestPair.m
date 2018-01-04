function [assemD]=TestPair(ensemble,spikeTrain2,n2,maxlag,Dc,reference_lag)
% this function tests if the two spike trains have repetitive patterns occurring more
% frequently than chance.

% ensemble := structure with the previously formed assembly and its spike train
% spikeTrain2 := spike train of the new unit to be tested for significance (candidate to be a new assembly member)
% n2 := new unit tested 
% maxlag := maximum lag to be tested
% Dc := size (in bins) of chunks in which I divide the spike trains to compute the variance (to reduce non stationarity effects on variance estimation)
% reference_lag := lag of reference; if zero or negative reference lag=-l 
%
%
%
%  © 2016 Russo, Durstewitz.
%  for information please contact eleonora.russo@zi-mannheim.de; daniel.durstewitz@zi-mannheim.de.
%
%  last update 11/01/2016
%  last update 13/02/2016 added correction for continuity approximation (ExpAB limit and Yates correction) 

    couple=[ensemble.Time-min(ensemble.Time(:)) spikeTrain2-min(spikeTrain2(:))]; % spike train pair I am going to test
    
    nu=2;
    ntp=size(couple,1);  %%%% trial length
   
    % I divide in parallel trials with 0/1 elements
    maxrate=max(max(couple));
    Zaa=cell(1,maxrate);
    for i=1:maxrate
        Zaa{i}=zeros(size(couple), 'uint8');
        Zaa{i}(couple>=i)=1;
        ExpABi(i)=prod(sum(Zaa{i},1))/size(couple,1);
    end
    
    % % % decide which is the lag with most coincidences (l_:=best lag)
    ctAB=nan(1,maxlag+1);
    ctAB_=nan(1,maxlag+1);
    for l=0:maxlag 
        trAB=[couple(1:end-maxlag,1), couple(l+1:end-maxlag+l,2)];
        trBA=[couple(l+1:end-maxlag+l,1), couple(1:end-maxlag,2)];
        ctAB(l+1)=sum(min(trAB'));
        ctAB_(l+1)=sum(min(trBA'));  
    end   
    
    if reference_lag<=0
        aus=[ctAB; ctAB_];
        [a,b]=max(aus(:));
        [I,J] = ind2sub(size(aus),b); 
        l_=(I==1)*(J-1)-(I==2)*(J-1);  %% I select l_    
    else
        Hab_l=[ctAB_(end:-1:2), ctAB];
        [a,b]=max(Hab_l(:));
        lags=-maxlag:maxlag;
        l_=lags(b);
        Hab=a;
        if l_<0
           l_ref=l_+reference_lag;
           Hab_ref=Hab_l(find(lags==l_ref));
        else
           l_ref=l_-reference_lag;
           Hab_ref=Hab_l(find(lags==l_ref));
        end               
    end
    

 
%% HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH


ExpAB=sum(ExpABi);

if  a==0  || ExpAB<=5  || ExpAB>=(min(sum(couple))-5)     %case of no coincidences or limit for the F asinptotical distribution (too few coincidences)
    assemD.elements=[ensemble.elements n2];
    assemD.lag=[ensemble.lag, 99];
    assemD.pr=[ensemble.pr 1];  % setting pr=1 the tested pair will be discarded as assembly
    assemD.Time=[];
    assemD.Noccurrences=[ensemble.Noccurrences 0];
else

    % % % construct the activation time series for the couple
    len=size(couple,1);        %%%% trial length
    Time=uint8(zeros(len,1));  %%%% activation vector
    if reference_lag<=0 
        if l_==0
            for i=1:maxrate  
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);
                Time(1:len)=Time(1:len)+sA(1:end).*sB(1:end);
            end
            TPrMTot=[0 ctAB(1); ctAB_(3) 0]; % matrix with #AB and #BA
        elseif l_>0      
            for i=1:maxrate  
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);         
                Time(1:len-l_)=Time(1:len-l_)+sA(1:end-l_).*sB(l_+1:end);
            end
            TPrMTot=[0 ctAB(J); ctAB_(J) 0]; % matrix with #AB and #BA 
        else
            for i=1:maxrate  
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);        
                Time(-l_+1:end)=Time(-l_+1:end)+sA(-l_+1:end).*sB(1:end+l_);
            end
            TPrMTot=[0 ctAB(J); ctAB_(J) 0]; % matrix with #AB and #BA 
        end
    
    else 
    
        if l_==0
            for i=1:maxrate  
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);
                Time(1:len)=Time(1:len)+sA(1:end).*sB(1:end);
            end       
        elseif l_>0      
            for i=1:maxrate  
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);         
                Time(1:len-l_)=Time(1:len-l_)+sA(1:end-l_).*sB(l_+1:end);
            end

        else
            for i=1:maxrate  
                sA=Zaa{i}(:,1);
                sB=Zaa{i}(:,2);        
                Time(-l_+1:end)=Time(-l_+1:end)+sA(-l_+1:end).*sB(1:end+l_);
            end
        end
        TPrMTot=[0 Hab; Hab_ref 0]; % matrix with #AB and #BA 
    end
    
%% --------------------------------------------------------------------%
    % % % HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    % I cut the spike train in stationary segments
    %%%%%%%%%%%%%%%%%%%%% chunking  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nch=ceil((size(couple,1)-maxlag)/Dc);
    Dc=floor((size(couple,1)-maxlag)/nch); %% new chunk size, this is to have all chunks of rougly the same size 
    chunked=cell(nch,1);
    
    % in order to take into account the same time series part that I used
    % for MargPr I cut the time series according to l_:
    couple_cut=nan((size(couple,1)-maxlag),2);

    if l_==0
        couple_cut=couple(1:end-maxlag,:);  
    elseif l_>0      
        couple_cut(:,1)=couple(1:end-maxlag,1);
        couple_cut(:,2)=couple(l_+1:end-maxlag+l_,2);
    else
        couple_cut(:,1)=couple(1-l_:end-maxlag-l_,1);
        couple_cut(:,2)=couple(1:end-maxlag,2);
    end
  

    for iii=1:nch-1
        chunked{iii}=couple_cut((1+Dc*(iii-1)):Dc*iii,:);
    end
    chunked{nch}=couple_cut((1+Dc*(nch-1)):end,:); % last chunk can be of slightly different size
   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    MargPr_t=cell(nch,maxrate);   %% number of spikes in each parallel process in each unit
    maxrate_t=nan(nch,1);
    ch_nn=nan(nch,1);
    for iii=1:nch  
        couple_t=chunked{iii};        
        maxrate_t(iii)=max(max(couple_t));
        ch_nn(iii)=size(chunked{iii},1); 
        Zaa_t=cell(1,maxrate_t(iii));
        for i=1:maxrate_t(iii)
            Zaa_t{i}=zeros(size(couple_t), 'uint8');
            Zaa_t{i}(couple_t>=i)=1;
        end

        for i=1:maxrate_t(iii)  
            sA=Zaa_t{i}(:,1);
            sB=Zaa_t{i}(:,2); 
            MargPr_t{iii,i}=[sum(sA); sum(sB)];   
        end                
    end
    
    %%%%%%%%%%%%%%%%%%%%% chunks  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n=ntp-maxlag;
    
    Mx0=cell(nch,maxrate);
    covABAB=cell(1,nch);
    covABBA=cell(1,nch);
    varT=cell(1,nch);
    covX=cell(1,nch);
    varX=cell(1,nch);
    varXtot=zeros(2);
    for iii=1:nch  
        maxrate_t=max(chunked{iii}(:));
        ch_n=ch_nn(iii);    
        % % % HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
        % % % evaluation of  #AB 
        for i=1:maxrate_t 
            if  ~isempty(MargPr_t{iii,i})
            Mx0{iii,i}= MargPr_t{iii,i} *ones(1,2); 
            end
        end  

        % % variance & covariance
        varT{iii}=zeros(nu);
        covABAB{iii}=cell(maxrate_t,maxrate_t);
        for i=1:maxrate_t 
               Mx0{iii,i}= MargPr_t{iii,i} *ones(1,2);              
                covABAB{iii}{i,i}=(Mx0{iii,i}.*Mx0{iii,i}'./ch_n).*(ch_n-Mx0{iii,i}).*(ch_n-Mx0{iii,i}')./(ch_n*(ch_n-1));  
                varT{iii}=varT{iii}+covABAB{iii}{i,i};
            for j=(i+1):maxrate_t 
               covABAB{iii}{i,j}=2*(Mx0{iii,j}.*Mx0{iii,j}'./ch_n).*(ch_n-Mx0{iii,i}).*(ch_n-Mx0{iii,i}')./(ch_n*(ch_n-1));
               varT{iii}=varT{iii}+covABAB{iii}{i,j};
            end
        end
        
        %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
        % % evaluation of  X=#AB-#BA 
        covX{iii}=zeros(nu);
        covABBA{iii}=cell(maxrate_t,maxrate_t);

        for i=1:maxrate_t 
            covABBA{iii}{i,i}=(Mx0{iii,i}.*Mx0{iii,i}'./ch_n).*(ch_n-Mx0{iii,i}).*(ch_n-Mx0{iii,i}')./(ch_n*(ch_n-1)^2);
            covX{iii}=covX{iii}+covABBA{iii}{i,i};       
            for j=(i+1):maxrate_t 
                  covABBA{iii}{i,j}=2*(Mx0{iii,j}.*Mx0{iii,j}'./ch_n).*(ch_n-Mx0{iii,i}).*(ch_n-Mx0{iii,i}')./(ch_n*(ch_n-1)^2);
                  covX{iii}=covX{iii}+covABBA{iii}{i,j};
            end
        end
        
        varX{iii}=varT{iii}+varT{iii}'-covX{iii}-covX{iii}';
        varXtot=varXtot+varX{iii};
    end
    %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    X=TPrMTot-TPrMTot';
    if abs(X(1,2))>0
       X=abs(TPrMTot-TPrMTot')-0.5; %Yates correction
    end
    
    if varXtot(1,2)==0  % if variance is zero
        prF=1;  
    else
        F=X.^2./varXtot;
        prF=fcdf(F(1,2),1,n,'upper'); 
%         prF=fcdf(F(1,2),1,2*double(maxrate)*n-1,'upper'); 
    end


    
%%
    %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    %HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    %All information about the assembly and test are returned 
    
    assemD.elements=[ensemble.elements n2];
    assemD.lag=[ensemble.lag, l_];
    assemD.pr=[ensemble.pr prF];
    assemD.Time=Time;
    assemD.Noccurrences=[ensemble.Noccurrences sum(Time)];
        

end
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
end






