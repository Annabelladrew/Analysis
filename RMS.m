% Created by Celine Lichtensteiger
% Calculates the RMS (Root Mean Square) of the pair wise differences of the fit and measurement
% and RMS(log), the RMS of the pair wise differences of the log of the fit and the log of the measurement
% The user can change the region over which the RMS and RMS(log) values are
% calculated (shown by vertical green lines).
% This is done only if the user loads his measured data.
%*********************************
global fitdata plotdata;

if strcmp(fitdata.FitDisplay.haxis,'L')
        plotdata.fit.compare.x=plotdata.fit.x.L(and(plotdata.fit.x.L<max(plotdata.measure.x),plotdata.fit.x.L>min(plotdata.measure.x)));
        plotdata.fit.compare.y=fitdata.ScalingIntensity*plotdata.fit.y(and(plotdata.fit.x.L<max(plotdata.measure.x),plotdata.fit.x.L>min(plotdata.measure.x)));
    
        plotdata.fit.compare.y=plotdata.fit.compare.y(and(plotdata.fit.compare.x<plotdata.RMS.LMax,plotdata.fit.compare.x>plotdata.RMS.LMin));
        plotdata.fit.compare.x=plotdata.fit.compare.x(and(plotdata.fit.compare.x<plotdata.RMS.LMax,plotdata.fit.compare.x>plotdata.RMS.LMin));
    
        plotdata.measure.compare.y=plotdata.measure.y(and(plotdata.measure.x<plotdata.RMS.LMax,plotdata.measure.x>plotdata.RMS.LMin));
        plotdata.measure.compare.x=plotdata.measure.x(and(plotdata.measure.x<plotdata.RMS.LMax,plotdata.measure.x>plotdata.RMS.LMin));
elseif strcmp(fitdata.FitDisplay.haxis,'2Theta')
        plotdata.fit.compare.x=plotdata.fit.x.TTheta(and(plotdata.fit.x.TTheta<max(plotdata.measure.x),plotdata.fit.x.TTheta>min(plotdata.measure.x)));
        plotdata.fit.compare.y=fitdata.ScalingIntensity*plotdata.fit.y(and(plotdata.fit.x.TTheta<max(plotdata.measure.x),plotdata.fit.x.TTheta>min(plotdata.measure.x)));
    
        plotdata.fit.compare.y=plotdata.fit.compare.y(and(plotdata.fit.compare.x<plotdata.RMS.TThetaMax,plotdata.fit.compare.x>plotdata.RMS.TThetaMin));
        plotdata.fit.compare.x=plotdata.fit.compare.x(and(plotdata.fit.compare.x<plotdata.RMS.TThetaMax,plotdata.fit.compare.x>plotdata.RMS.TThetaMin));
    
        plotdata.measure.compare.y=plotdata.measure.y(and(plotdata.measure.x<plotdata.RMS.TThetaMax,plotdata.measure.x>plotdata.RMS.TThetaMin));
        plotdata.measure.compare.x=plotdata.measure.x(and(plotdata.measure.x<plotdata.RMS.TThetaMax,plotdata.measure.x>plotdata.RMS.TThetaMin));
else warning('Error in fitdata.FitDisplay.haxis value (2TTheta or L)');
end; 
      
if and(size(plotdata.measure.compare.x,1)>1,size(plotdata.fit.compare.x,1)>1)

    plotdata.measure.compare.y=interp1(plotdata.measure.compare.x,plotdata.measure.compare.y,plotdata.fit.compare.x);
    plotdata.measure.compare.y(isnan(plotdata.fit.compare.y))=[];
    plotdata.fit.compare.y(isnan(plotdata.fit.compare.y))=[];
    plotdata.fit.compare.y(isnan(plotdata.measure.compare.y))=[];
    plotdata.measure.compare.y(isnan(plotdata.measure.compare.y))=[];

    plotdata.RMS.lin=sum((plotdata.measure.compare.y'-plotdata.fit.compare.y).^2)/size(plotdata.fit.compare.y,2);
    plotdata.RMS.log=sum((log(plotdata.measure.compare.y')-log(plotdata.fit.compare.y)).^2)/size(plotdata.fit.compare.y,2);

else
    plotdata.RMS.lin=0;
    plotdata.RMS.log=0;
end;

uicontrol(plotdata.MainPanel.fig,'Style','text', ...
'String',['RMS=',num2str(plotdata.RMS.lin)],...
'Position',[0,0,100,14]);

uicontrol(plotdata.MainPanel.fig,'Style','text', ...
'String',['RMS(log)=',num2str(plotdata.RMS.log)],...
'Position',[0,20,100,14]);