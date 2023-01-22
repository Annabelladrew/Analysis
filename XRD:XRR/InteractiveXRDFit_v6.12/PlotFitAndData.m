% Created by Celine Lichtensteiger
% Plots the measured and calculated XRD intensities
% Plots the c-axis profile of the simulated heterostructure
%*********************************
global plotdata;

figure(plotdata.FitDisplay.fig);
plotdata.hAxes=subplot(1,1,1);
    
if strcmp(plotdata.FitDisplay.haxis,'L')
    plotdata.fit.x.L=2*plotdata.Substrate.d*sin(plotdata.fit.x.TTheta/2*pi/180)/plotdata.lambda;
    plotdata.hplot=semilogy(plotdata.fit.x.L,plotdata.fit.y*plotdata.ScalingIntensity,'b',plotdata.measure.x,plotdata.measure.y,'r');    
    axis([plotdata.LMin plotdata.LMax 10^-9 1]);
    xlabel('L'); ylabel('Intensity (arbitrary units)');

    if size(plotdata.measure.x,1)~=0,       
    
    figure(plotdata.FitDisplay.fig);
    
    end;
        

elseif strcmp(plotdata.FitDisplay.haxis,'2Theta')
    semilogy(plotdata.fit.x.TTheta,plotdata.fit.y*plotdata.ScalingIntensity,'b',plotdata.measure.x,plotdata.measure.y,'r');    
    axis([plotdata.TThetaMin plotdata.TThetaMax 10^-9 1]);
    xlabel('2Theta ({\circ})'); ylabel('Intensity (arbitrary units)');
    
    if size(plotdata.measure.x,1)~=0,   
    
    figure(plotdata.FitDisplay.fig);

        end;

else warning('Error in PlotFitAndData')    

end;

figure(plotdata.dDisplay.fig); %clf
color=['rx';'bx';'mx';'cx';'kx';'gx'];
LegendPlot=[];
if plotdata.Material(1).N~=0,
    p(1)=plot([1:plotdata.Material(1).N],plotdata.d(1:plotdata.Material(1).N),color(1,:),'DisplayName',char(plotdata.Material(1).Type)); hold on;
    LegendPlot=[LegendPlot,p(1)];
    legend(LegendPlot);
end;
x=plotdata.Material(1).N;
for rep=1:plotdata.Repetition.N;
    for k=2:5,
        if plotdata.Material(k).N~=0,
            if rep==1,
                p(k)=plot([x+1:x+plotdata.Material(k).N],plotdata.d(x+1:x+plotdata.Material(k).N),color(k,:),'DisplayName',char(plotdata.Material(k).Type)); hold on;
                LegendPlot=[LegendPlot,p(k)];
                
            else
                plot([x+1:x+plotdata.Material(k).N],plotdata.d(x+1:x+plotdata.Material(k).N),color(k,:));
                hold on;
            end;
            x=x+plotdata.Material(k).N;
            legend(LegendPlot);
        end;
    end;
end;
if plotdata.Material(6).N~=0,
    p(6)=plot([x+1:x+plotdata.Material(6).N],plotdata.d(x+1:x+plotdata.Material(6).N),color(6,:),'DisplayName',char(plotdata.Material(6).Type)); hold on;
    LegendPlot=[LegendPlot,p(6)];
    legend(LegendPlot);
end;
xlabel('Sample thickness (u.c.)'); ylabel('c (Angstroms)');
hold off