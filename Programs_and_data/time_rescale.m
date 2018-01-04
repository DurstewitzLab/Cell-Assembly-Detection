function [activity_tbref]=time_rescale(activity, ta, tbref)
%  to resample assembly activation
%  © 2016 Russo, Durstewitz.
%  for information please contact eleonora.russo@zi-mannheim.de; daniel.durstewitz@zi-mannheim.de.
%
%  last update 11/01/2016

    original_step=ta(2)-ta(1);
    reference_step=tbref(2)-tbref(1);

    if original_step>=reference_step    
        activity_tbref = interp1(ta,activity,tbref);    

    elseif original_step<reference_step  
        
        [activity_tbref] = whist(ta, activity, tbref);
        activity_tbref=activity_tbref';
   
    end

end



function [histw] = whist(v, w, tbref)
% Inputs:
% v - values
% w - weights
% tbref - time on which I want to scale
% v=ta;
% w=activity;

    reference_step=tbref(2)-tbref(1);
    
    tbref_edges=tbref-reference_step/2;
    tbref_edges=[tbref_edges, tbref_edges(end)+reference_step];
     
    minv=min(tbref_edges);
    aus=(v-minv)/reference_step;
    v(aus<0)=[];
    w(aus<0)=[];
    
    subs=floor((v-minv)/reference_step)+1;
    histw = accumarray(subs,w);
    


    
    [~,minind]=min(abs(tbref-min(v)));
    if minind>1
        histw=[zeros(minind-1,1);histw];
    end
    [~,maxind]=min(abs(tbref-max(v)));
    if (maxind+minind-1)<length(tbref)
        histw=[histw; zeros(length(tbref)-(maxind+minind-1),1)];
    end       
    


end









