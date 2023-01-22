% Created by Celine Lichtensteiger
% Plots the measured and calculated XRD intensities
% Plots the lines delimitating the region over which the RMS is calculated
% Plots the c-axis profile of the simulated heterostructure
%*********************************
global fitdata plotdata;

figure(plotdata.FitDisplay.fig);
plotdata.hAxes=subplot(1,1,1);
    
if strcmp(fitdata.FitDisplay.haxis,'L')
    plotdata.fit.x.L=2*fitdata.Substrate.d*sin(plotdata.fit.x.TTheta/2*pi/180)/plotdata.lambda;
    plotdata.hplot=semilogy(plotdata.fit.x.L,plotdata.fit.y*fitdata.ScalingIntensity,'b',plotdata.measure.x,plotdata.measure.y,'r');    
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
        

elseif strcmp(fitdata.FitDisplay.haxis,'2Theta')
    semilogy(plotdata.fit.x.TTheta,plotdata.fit.y*fitdata.ScalingIntensity,'b',plotdata.measure.x,plotdata.measure.y,'r');    
    axis([plotdata.TThetaMin plotdata.TThetaMax 10^-9 1]);
    xlabel('2Theta ({\circ})'); ylabel('Intensity (arbitrary units)');
    
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

figure(plotdata.dDisplay.fig); %clf
color=['rx';'bx';'mx';'cx';'kx';'gx'];
LegendPlot=[];
if fitdata.Material(1).N~=0,
    p(1)=plot([1:fitdata.Material(1).N],plotdata.d(1:fitdata.Material(1).N),color(1,:),'DisplayName',char(fitdata.Material(1).Type)); hold on;
    LegendPlot=[LegendPlot,p(1)];
    legend(LegendPlot);
end;
x=fitdata.Material(1).N;
for rep=1:fitdata.Repetition.N;
    for k=2:5,
        if fitdata.Material(k).N~=0,
            if rep==1,
                p(k)=plot([x+1:x+fitdata.Material(k).N],plotdata.d(x+1:x+fitdata.Material(k).N),color(k,:),'DisplayName',char(fitdata.Material(k).Type)); hold on;
                LegendPlot=[LegendPlot,p(k)];
                
            else
                plot([x+1:x+fitdata.Material(k).N],plotdata.d(x+1:x+fitdata.Material(k).N),color(k,:));
                hold on;
            end;
            x=x+fitdata.Material(k).N;
            legend(LegendPlot);
        end;
    end;
end;
if fitdata.Material(6).N~=0,
    p(6)=plot([x+1:x+fitdata.Material(6).N],plotdata.d(x+1:x+fitdata.Material(6).N),color(6,:),'DisplayName',char(fitdata.Material(6).Type)); hold on;
    LegendPlot=[LegendPlot,p(6)];
    legend(LegendPlot);
end;
xlabel('Sample thickness (u.c.)'); ylabel('c (Angstroms)');
hold off