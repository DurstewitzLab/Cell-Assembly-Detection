function [AAspM,AAasspikes]=assembly_rasterplot(As_across_bins,assembly_activity,spM,show)
%  © 2016 Russo, Durstewitz.
%  for information please contact eleonora.russo@zi-mannheim.de; daniel.durstewitz@zi-mannheim.de.
%
%  last update 11/01/2016
%
    AAspM=cell(size(As_across_bins,2),1);
    AAasspikes=cell(size(As_across_bins,2),1);

    for ii=1:size(As_across_bins,2)

        Assembly_activity=assembly_activity{ii};
        A=As_across_bins{ii};

        activation=Assembly_activity(:,2);
        tb=Assembly_activity(:,1);

        Aelements=A.elements;
        Alags=A.lag;
        AspM=spM(Aelements,:);
        ATime=A.Time;

        ASassembly=nan(size(AspM));
        for i=1:size(AspM,1)  
            aus=[zeros(Alags(i),1); ATime(1:end-Alags(i))];
            activ_bins=tb(find(aus));
            unit_as_spikes=[];
            for j=1:size(activ_bins,1)       
                unit_as_spikes=[unit_as_spikes, AspM(i,find(AspM(i,:)>=activ_bins(j,1)-A.bin/2 & AspM(i,:)<activ_bins(j,1)+A.bin/2 ))];
            end
            ASassembly(i,1:length(unit_as_spikes))=unit_as_spikes;    
        end

        
        AAspM{ii}=AspM;
        AAasspikes{ii}=ASassembly;

        if show
            clf
            subplot(2,1,1)
            for i=1:size(AspM,1)
                aus=AspM(i,:);
                aus(isnan(aus))=[];
                plot(aus,i*ones(size(aus)),'.b')
                hold on

                aus=ASassembly(i,:);
                aus(isnan(aus))=[];
                plot(aus,i*ones(size(aus)),'.r')    
            end

            subplot(2,1,2)
            plot(tb(1:end),activation)

            fprintf('press a key\n');
            pause
        end

    end



end
