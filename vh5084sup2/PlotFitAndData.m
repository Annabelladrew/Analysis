% Created by Celine Lichtensteiger
% Plots the measured and calculated XRD intensities
% Plots the lines delimitating the region over which the RMS is calculated
% Plots the c-axis profile of the simulated heterostructure
%*********************************
global plotdata;

figure(plotdata.FitDisplay.fig);
plotdata.hAxes=subplot(1,1,1);
    
if strcmp(plotdata.FitDisplay.haxis,'L')
    plotdata.fit.x.L=2*plotdata.Substrate.c*sin(plotdata.fit.x.TTheta/2*pi/180)/plotdata.lambda;
    plotdata.hplot=semilogy(plotdata.fit.x.L,plotdata.fit.y*plotdata.ScalingIntensity,'b',plotdata.measure.x,plotdata.measure.y,'r');    
    axis([plotdata.LMin plotdata.LMax 10^-9 1]);
    xlabel('L'); ylabel('Intensity (arbitrary units)');

    if size(plotdata.measure.x,1)~=0,   
        
    %delimitation of RMS calculation
    if isempty(plotdata.RMS.LMin),
        plotdata.RMS.LMin=max(plotdata.LMin,min(plotdata.measure.x));
        set(plotdata.chooseRMSLMin,'String',num2str(plotdata.RMS.LMin));
    end;
    if isempty(plotdata.RMS.LMax),
        plotdata.RMS.LMax=min(plotdata.LMax,max(plotdata.measure.x));
        set(plotdata.chooseRMSLMax,'String',num2str(plotdata.RMS.LMax));
    end;
    
    figure(plotdata.FitDisplay.fig);

    plotdata.hVerticalLines = line([plotdata.RMS.LMin,plotdata.RMS.LMin],get(plotdata.hAxes,'Ylim'),'Color','green');
    plotdata.hVerticalLines = [plotdata.hVerticalLines line([plotdata.RMS.LMax,plotdata.RMS.LMax],get(plotdata.hAxes,'Ylim'),'Color','green')];    
    RMS
    
    end;
        

elseif strcmp(plotdata.FitDisplay.haxis,'2Theta')
    semilogy(plotdata.fit.x.TTheta,plotdata.fit.y*plotdata.ScalingIntensity,'b',plotdata.measure.x,plotdata.measure.y,'r');    
    axis([plotdata.TThetaMin plotdata.TThetaMax 10^-9 1]);
    xlabel('2Theta (°)'); ylabel('Intensity (arbitrary units)');
    
    if size(plotdata.measure.x,1)~=0,   
        
    %delimitation of RMS calculation
    if isempty(plotdata.RMS.TThetaMin),
        plotdata.RMS.TThetaMin=max(plotdata.TThetaMin,min(plotdata.measure.x));
        set(plotdata.chooseRMSTThetaMin,'String',num2str(plotdata.RMS.TThetaMin));
    end;
    if isempty(plotdata.RMS.TThetaMax),
        plotdata.RMS.TThetaMax=min(plotdata.TThetaMax,max(plotdata.measure.x));
        set(plotdata.chooseRMSTThetaMax,'String',num2str(plotdata.RMS.TThetaMax));
    end;
    
    figure(plotdata.FitDisplay.fig);

    plotdata.hVerticalLines = line([plotdata.RMS.TThetaMin,plotdata.RMS.TThetaMin],get(plotdata.hAxes,'Ylim'),'Color','green');
    plotdata.hVerticalLines = [plotdata.hVerticalLines line([plotdata.RMS.TThetaMax,plotdata.RMS.TThetaMax],get(plotdata.hAxes,'Ylim'),'Color','green')];    
    RMS
        end;

else warning('Error in PlotFitAndData')    

end;

figure(plotdata.cDisplay.fig); clf
color=['rx';'bx';'mx';'cx';'kx'];
LegendPlot=[];
if plotdata.Material(1).N~=0,
    p(1)=plot([1:plotdata.Material(1).N],plotdata.c(1:plotdata.Material(1).N),color(1,:),'DisplayName','Bottom Layer'); hold on;
    LegendPlot=[LegendPlot,p(1)];
    legend(LegendPlot);
end;
x=plotdata.Material(1).N;
for rep=1:plotdata.Repetition.N;
    for k=2:4,
        if plotdata.Material(k).N~=0,
            if rep==1,
                p(k)=plot([x+1:x+plotdata.Material(k).N],plotdata.c(x+1:x+plotdata.Material(k).N),color(k,:),'DisplayName',['Material ',num2str(k-1)]); hold on;
                LegendPlot=[LegendPlot,p(k)];
                legend(LegendPlot);
            else
                plot([x+1:x+plotdata.Material(k).N],plotdata.c(x+1:x+plotdata.Material(k).N),color(k,:)); hold on;
            end;
            x=x+plotdata.Material(k).N;
        end;
    end;
end;
if plotdata.Material(5).N~=0,
    p(5)=plot([x+1:x+plotdata.Material(5).N],plotdata.c(x+1:x+plotdata.Material(5).N),color(5,:),'DisplayName','Top Layer'); hold on;
    LegendPlot=[LegendPlot,p(5)];
    legend(LegendPlot);
end;
xlabel('Sample thickness (u.c.)'); ylabel('c (Angstroms)');
hold off