% Created by Celine Lichtensteiger
% Calculates the contribution from each layer material to:
% plotdata.Material(k).expzm
% plotdata.ThicknessBelow
% plotdata.c
%*********************************

function expzm(k)
global plotdata;
    if strcmp(plotdata.Material(k).cDistribution,'constant')
        if strcmp(plotdata.orientation,'(001)')
            plotdata.Material(k).Thickness=plotdata.Material(k).N*plotdata.Material(k).c;
            plotdata.Material(k).F=FLayer(k,plotdata.Material(k).c);
        elseif strcmp(plotdata.orientation,'(111)') 
            plotdata.Material(k).Thickness=plotdata.Material(k).N*plotdata.Material(k).d;
            plotdata.Material(k).F=FLayer(k,plotdata.Material(k).d);
        else warning('Error with plotdata.orientation in function expzm');
        end;
        for zm=1:plotdata.Material(k).N,  
            if strcmp(plotdata.orientation,'(001)')
                plotdata.ThicknessBelow=plotdata.ThicknessBelow+plotdata.Material(k).c;
                plotdata.c=[plotdata.c plotdata.Material(k).c];
            elseif strcmp(plotdata.orientation,'(111)')
                plotdata.ThicknessBelow=plotdata.ThicknessBelow+plotdata.Material(k).d;
                plotdata.c=[plotdata.c plotdata.Material(k).d];
            else warning('Error with plotdata.orientation in function expzm');
            end;
            ThicknessAbove=plotdata.TotalThickness-plotdata.ThicknessBelow;
            plotdata.Material(k).expzm=plotdata.Material(k).expzm+plotdata.Material(k).F.*exp(1i*plotdata.Q.*plotdata.ThicknessBelow).*exp(-(ThicknessAbove)./plotdata.mu);
        end;
    elseif strcmp(plotdata.Material(k).cDistribution,'exp')
        plotdata.Material(k).Thickness=0;
        plotdata.Material(k).cexp=plotdata.Material(k).Meancexp-plotdata.Material(k).aexp/plotdata.Material(k).N*(1-exp(plotdata.Material(k).N/plotdata.Material(k).bexp))/(1-exp(1/plotdata.Material(k).bexp));
        for zm=1:plotdata.Material(k).N,
            czm=plotdata.Material(k).aexp*exp((zm-1)/plotdata.Material(k).bexp)+plotdata.Material(k).cexp;
            plotdata.ThicknessBelow=plotdata.ThicknessBelow+czm;
            ThicknessAbove=plotdata.TotalThickness-plotdata.ThicknessBelow;
            plotdata.Material(k).F=FLayer(k,czm);
            plotdata.Material(k).expzm=plotdata.Material(k).expzm+plotdata.Material(k).F.*exp(1i*plotdata.Q.*plotdata.ThicknessBelow).*exp(-(ThicknessAbove)./plotdata.mu);
            plotdata.Material(k).Thickness=plotdata.Material(k).Thickness+czm;
            plotdata.c=[plotdata.c czm];
        end;
    else
        warning('Error with toggle function for cDistribution in expzm')
    end;