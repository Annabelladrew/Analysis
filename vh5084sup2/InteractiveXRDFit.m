% Created by Celine Lichtensteiger
% This is the heart of the program
% In initializes all the variables and creates the Graphical User Interface
%*********************************

function InteractiveXRDFit(entry)

global plotdata;

if nargin == 0
   entry = 'makeGraph';
end

switch entry

case 'makeGraph'
%%Default values
   
   % lambda
   plotdata.lambda=1.5406;
   % Substrate
   load('Substrates.mat');
   plotdata.Substrate.c=3.905;                              
   plotdata.Substrate.g=gsubstrateSrTiO3001;                   
   plotdata.Substrate.N001=2e4;                             
   plotdata.Substrate.N111=plotdata.Substrate.N001*sqrt(3); %using more layers when calculating the substrate in (111) to avoid finite size oscillations
   plotdata.Substrate.d=plotdata.Substrate.c/sqrt(3);
   plotdata.Substrate.Type='SrTiO3';                       
   % Bottom Layer
   plotdata.Material(1).aexp=0;
   plotdata.Material(1).bexp=1;
   plotdata.Material(1).c=4;
   plotdata.Material(1).cDistribution='constant';
   plotdata.Material(1).d=2.309;%4/sqrt(3);
   plotdata.Material(1).Meancexp=4;
   plotdata.Material(1).N=0;
   plotdata.Material(1).Polarization=0;
   plotdata.Material(1).Type='none';
   plotdata.Material(1).xPST=0;
   plotdata.Material(1).xPZT=0;
   plotdata.Material(1).cexp=plotdata.Material(1).Meancexp-plotdata.Material(1).aexp/plotdata.Material(1).N*(1-exp(plotdata.Material(1).N/plotdata.Material(1).bexp))/(1-exp(1/plotdata.Material(1).bexp));
   % Min-Max:
   plotdata.Material(1).aexpMax=10;      
   plotdata.Material(1).aexpMin=-10;     
   plotdata.Material(1).bexpMax=1000;    
   plotdata.Material(1).bexpMin=-1000;   
   plotdata.Material(1).cMax=14;         
   plotdata.Material(1).cMin=0;          
   plotdata.Material(1).MeancexpMax=14;  
   plotdata.Material(1).MeancexpMin=0;   
   plotdata.Material(1).dMax=14;         
   plotdata.Material(1).dMin=0;          
   plotdata.Material(1).NMax=1000;       
   plotdata.Material(1).NMin=0;         
   plotdata.Material(1).xPSTMax=1;
   plotdata.Material(1).xPSTMin=0;
   plotdata.Material(1).xPZTMax=1;
   plotdata.Material(1).xPZTMin=0;
   plotdata.Material(1).PolarizationMax=2;     
   plotdata.Material(1).PolarizationMin=-2;    
   % Material 1 in superlattice
   plotdata.Material(2)=plotdata.Material(1);
   % Material 2 in superlattice
   plotdata.Material(3)=plotdata.Material(1);
   % Material 3 in superlattice
   plotdata.Material(4)=plotdata.Material(1);
   % Top Layer
   plotdata.Material(5)=plotdata.Material(1);
   % cDisplay
   plotdata.cDisplay.Width=660;
   plotdata.cDisplay.Height=400;
   plotdata.cDisplay.centerfig=centerfig(plotdata.cDisplay.Width,plotdata.cDisplay.Height)+[400,-250,0,0];
   % RMS
   plotdata.RMS.hVerticalLines=[];
   plotdata.RMS.LMax=[];
   plotdata.RMS.LMin=[];
   plotdata.RMS.TThetaMax=[];
   plotdata.RMS.TThetaMin=[];
   % FitDisplay
   plotdata.FitDisplay.haxis='2Theta';
   plotdata.FitDisplay.Width=660;
   plotdata.FitDisplay.Height=400;
   plotdata.FitDisplay.centerfig=centerfig(plotdata.FitDisplay.Width,plotdata.FitDisplay.Height)-[400,0,0,0];
   plotdata.LMin=0;
   plotdata.LMax=2*plotdata.Substrate.c/plotdata.lambda; 
   plotdata.ScalingIntensity=1;
   plotdata.TThetaMin=0;
   plotdata.TThetaMax=180;
   % Main Panel
   plotdata.MainPanel.Width=792;    
   plotdata.MainPanel.Height=420;   
   plotdata.MainPanel.centerfig=centerfig(plotdata.MainPanel.Width,plotdata.MainPanel.Height)+[400,+250,0,0];
   plotdata.orientation='(001)';
   plotdata.Repetition.N=1;
   plotdata.Repetition.Max=20;
   plotdata.Repetition.Min=0;
   % Misc
   plotdata.fit.compare.x=[];
   plotdata.fit.compare.y=[];
   plotdata.fit.x.TTheta=[0:0.01:180]';
   plotdata.fit.x.L=2*plotdata.Substrate.c*sin(plotdata.fit.x.TTheta/2*pi/180)/plotdata.lambda;
   plotdata.fit.y=10^10+zeros(size([0:0.01:180]'));
   plotdata.measure.compare.y=[];
   plotdata.measure.x=[];
   plotdata.measure.y=[];
   plotdata.measure.filename='none';
   plotdata.mu=1.5e4; %Penetration depth
   plotdata.prog = mfilename;
   plotdata.Q=2*pi*(2*sin([0:0.01:180]'/2*pi/180)/1.5406)'; %Momentum transfert
   plotdata.QQ=plotdata.Q.^2;
   plotdata.TotalThickness=plotdata.Substrate.N001*plotdata.Substrate.c+100000;
   plotdata.c=[];
   %% Create Main Panel
   plotdata.MainPanel.fig = figure('Position',plotdata.MainPanel.centerfig,...
      'Resize','on',...
      'NumberTitle','off',...
      'Name','XRD fitting parameters @Celine Lichtensteiger',...
      'Interruptible','off',...
      'Menubar','none',...
      'Color',get(0,'DefaultUIControlBackgroundColor'));
   set(plotdata.MainPanel.fig,'toolbar','figure');
   set(plotdata.MainPanel.fig,'menubar','figure');
     %==Text====================================
    uicontrol(plotdata.MainPanel.fig,'style','text',...
       'string','Simulates XRD diffraction for heterostructures.',...
       'position',[(plotdata.MainPanel.Width-450)/2,plotdata.MainPanel.Height-30,450,25]);
   
   %==file name====================================
   plotdata.filenameWrite=uicontrol(plotdata.MainPanel.fig,'Style','text', ...
      'String','Load the data you want to fit by pressing "Load data":',...
      'Position',[(plotdata.MainPanel.Width-400)/2,plotdata.MainPanel.Height-45,400,25]); 
   
   uicontrol(plotdata.MainPanel.fig,'Style','pushbutton',...
      'Position',[plotdata.MainPanel.Width-130 plotdata.MainPanel.Height-45 58 25],...
      'String','Load data',...
      'Interruptible','off',...
      'BusyAction','cancel',...
      'Callback',[plotdata.prog,' LoadDataFilename']);
  
   %===The quit button===============================
   uicontrol(plotdata.MainPanel.fig,'Style','pushbutton',...
      'Position',[plotdata.MainPanel.Width-65 plotdata.MainPanel.Height-25 45 25],...
      'String','Quit',...
      'Interruptible','off',...
      'BusyAction','cancel',...
      'Callback',[plotdata.prog,' quit']);
  
  %===The export button===============================
   uicontrol(plotdata.MainPanel.fig,'Style','pushbutton',...
      'Position',[plotdata.MainPanel.Width-75 plotdata.MainPanel.Height-45 55 25],...
      'String','Export fit',...
      'Interruptible','off',...
      'BusyAction','cancel',...
      'Callback',[plotdata.prog,' export']);
  
    uicontrol(plotdata.MainPanel.fig,'style','text',...
       'string','Export the fit as a .csv file.',...
       'position',[(plotdata.MainPanel.Width-450)/2,plotdata.MainPanel.Height-60,450,25]);
  
%==Substrate chooser=========================
plotdata.Substrate.Panel=uipanel(plotdata.MainPanel.fig,'Title','Substrate:',...
    'BorderType','etchedout',...
    'position',[0.005,0.75,0.16,0.1]);   
plotdata.Substrate.choose=uicontrol(plotdata.Substrate.Panel,'Style','popupmenu', ...
      'String','none|DyScO3 (pc c=3.9403A)|GdScO3 (pc c=3.9636A)|KTaO3 (c=3.989A)|LaAlO3 (c=3.789A)|LaSrAlO4 (c=12.6377A)|LSAT (c=3.868A)|NdAlO3 (c=3.74A)|NdGaO3 (c=3.864A)|Si (c=5.4307A)|SrTiO3 (c=3.905A)|TbScO3 (c=7.917A)|YAlO3 (c=3.71A)',...
      'value',11,...
      'Units','normalized',...
      'Position',[0,0.3,1,0.5],...
      'visible','on',...
      'CallBack',[plotdata.prog,' chooseSubstrate'] ); 

  %==Superlattice====================================
  uipanel(plotdata.MainPanel.fig,'BorderType','etchedout','position',[0.3308,0.0001,0.4988,0.86]);
  
  uicontrol(plotdata.MainPanel.fig,'style','text',...
       'string','|-------------------------------------------------------------------------------------------------------------------------------|',...
       'position',[260,30,400,20]);
   
   %slider for Repetition
   uicontrol(plotdata.MainPanel.fig,'Style','text', ...
      'String','Number of repetitions of the superlattice (Material 1 / Material 2 / Material 3):',...
      'Position',[260,20,400,20]);
   plotdata.RepetitionWrite = uicontrol(plotdata.MainPanel.fig,'Style','edit',...
      'String',num2str(plotdata.Repetition.N),...
      'Position',[475,2,50,20],...
      'callback', [plotdata.prog,' writeRepetition']);
   plotdata.RepetitionSlide = uicontrol(plotdata.MainPanel.fig,...
      'Style','slider',...
      'Min' ,0,'Max',20, ...
      'Position',[370,0,100,20], ...
      'Value', plotdata.Repetition.N,...
      'SliderStep',[1/20 1/20], ...
      'CallBack', [plotdata.prog,' slideRepetition']);
   
  %==Material chooser and slider ========================= 
  for k=1:5,
      
   plotdata.Material(k).Panel=uipanel(plotdata.MainPanel.fig,... 
       'BorderType','etchedout',...
       'position',[0.17+(k-1)*0.165,0.15,0.16,0.7]);
   
      plotdata.Material(k).choose=uicontrol(plotdata.Material(k).Panel,'Style','popupmenu', ...
     'String','none|AlO2|BaO|BaTiO3|BiFeO3|CaCuO2|LaAlO3|La2CuO4|LaFeO3|LaMnO3|LaNiO3|La2NiMnO6|LaO|LSMO (La0.67Sr0.33MnO3)|MnO|MnO2|MnTiO3|NiO2|NdNiO3|Nd2NiMnO6|NdO|PbO|PbTiO3|(Pb_x,Sr_{1-x})TiO3|Pb(Zr_x,Ti_{1-x})O3|PrBa2Cu3O7|RuO2|SmNiO3|SrO|SrO2|SrRuO3|SrTiO3|SrVO3|TiO2|VO2|YBa2Cu3O7|ZrO2',...      
      'value',1,...
      'Units','normalized',...
      'visible','on',...
      'Position',[0,0.88,1,0.1],...
      'CallBack',[plotdata.prog,' chooseLayerType(',num2str(k),')'] );
  
  %==polarization chooser and slider =========================
  plotdata.Material(k).TBPolarization=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','Polarization:',...
      'BorderType','etchedout',...
      'Position',[0,0.74,1,0.14]);
  
  plotdata.Material(k).TBPolarizationNo = uicontrol('Style','radiobutton','String','No',...
    'pos',[0,0,100,30],'parent',plotdata.Material(k).TBPolarization);
  
  plotdata.Material(k).TBPolarizationYes = uicontrol('Style','radiobutton','String','Yes',...
    'pos',[40,0,100,30],'parent',plotdata.Material(k).TBPolarization);

  plotdata.Material(k).PolarizationWrite = uicontrol(plotdata.Material(k).TBPolarization,'Style','edit',...
      'String',num2str(plotdata.Material(k).Polarization),...
      'Units','normalized',...
      'Position',[0.7,0.5,0.3,0.6],...
      'callback', [plotdata.prog,' writePolarization(',num2str(k),')']);
  
  plotdata.Material(k).PolarizationSlide = uicontrol(plotdata.Material(k).TBPolarization,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).PolarizationMin,'Max',plotdata.Material(k).PolarizationMax, ...
      'Units','normalized',...
      'Position',[0.7,0,0.3,0.6], ...
      'Value', plotdata.Material(k).Polarization,...
      'SliderStep',[0.1/(plotdata.Material(k).PolarizationMax-plotdata.Material(k).PolarizationMin) 0.1/(plotdata.Material(k).PolarizationMax-plotdata.Material(k).PolarizationMin)], ...
      'CallBack', [plotdata.prog,' slidePolarization(',num2str(k),')']);  
  
  set(plotdata.Material(k).TBPolarization,'SelectionChangeFcn',@togglePolarization);
  
  %==PST concentration chooser and slider =========================
  plotdata.Material(k).xPSTpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.74,1,0.12]);

  plotdata.Material(k).xPSTWrite = uicontrol(plotdata.Material(k).xPSTpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xPST),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexPST(',num2str(k),')']);
  
  plotdata.Material(k).xPSTSlide = uicontrol(plotdata.Material(k).xPSTpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xPSTMin,'Max',plotdata.Material(k).xPSTMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xPST,...
      'SliderStep',[0.1/(plotdata.Material(k).xPSTMax-plotdata.Material(k).xPSTMin) 0.1/(plotdata.Material(k).xPSTMax-plotdata.Material(k).xPSTMin)], ...
      'CallBack', [plotdata.prog,' slidexPST(',num2str(k),')']);
  
  %==PZT concentration chooser and slider =========================
  plotdata.Material(k).xPZTpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.74,1,0.12]);

  plotdata.Material(k).xPZTWrite = uicontrol(plotdata.Material(k).xPZTpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xPZT),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexPZT(',num2str(k),')']);
  
  plotdata.Material(k).xPZTSlide = uicontrol(plotdata.Material(k).xPZTpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xPZTMin,'Max',plotdata.Material(k).xPZTMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xPZT,...
      'SliderStep',[0.1/(plotdata.Material(k).xPZTMax-plotdata.Material(k).xPZTMin) 0.1/(plotdata.Material(k).xPZTMax-plotdata.Material(k).xPZTMin)], ...
      'CallBack', [plotdata.prog,' slidexPZT(',num2str(k),')']);
 
   %slider for N
   plotdata.Material(k).chooseN=uipanel(plotdata.Material(k).Panel, ...
      'Title','N:',...
      'BorderType','etchedout',...
      'Position',[0,0.6,1,0.12]);
   plotdata.Material(k).NWrite = uicontrol(plotdata.Material(k).chooseN,'Style','edit',...
      'String',num2str(plotdata.Material(k).N),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writeN(',num2str(k),')']);
   plotdata.Material(k).NSlide = uicontrol(plotdata.Material(k).chooseN,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).NMin,'Max',plotdata.Material(k).NMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).N,...
      'SliderStep',[1/(plotdata.Material(k).NMax-plotdata.Material(k).NMin) 1/(plotdata.Material(k).NMax-plotdata.Material(k).NMin)], ...
      'CallBack', [plotdata.prog,' slideN(',num2str(k),')']);  

   %slider for d in 111 orientation
   plotdata.Material(k).choosed=uipanel(plotdata.Material(k).Panel, ...
      'Title','d:',...
      'visible','off',...
      'BorderType','etchedout',...
      'Position',[0,0.4,1,0.15]);
   plotdata.Material(k).dWrite = uicontrol(plotdata.Material(k).choosed,'Style','edit',...
      'String',num2str(plotdata.Material(k).d),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writed(',num2str(k),')']);
   plotdata.Material(k).dSlide = uicontrol(plotdata.Material(k).choosed,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).dMin,'Max',plotdata.Material(k).dMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).d,...
      'SliderStep',[0.001/(plotdata.Material(k).dMax-plotdata.Material(k).dMin) 0.001/(plotdata.Material(k).dMax-plotdata.Material(k).dMin)], ...
      'CallBack', [plotdata.prog,' slided(',num2str(k),')']);  

%== c: Toggle between constant and exponential distribution =========================
% Create the button group.
plotdata.Material(k).TBc = uibuttongroup(plotdata.Material(k).Panel,'visible','on',...
    'Title','Choose c distribution:',...
    'BorderType','etchedout',...
    'Position',[0 0 1 0.55]);
% Create two buttons in the button group.
plotdata.Material(k).TBcConst = uicontrol('Style','radiobutton','String','c=constant:',...
    'pos',[0,120,100,30],'parent',plotdata.Material(k).TBc);

plotdata.Material(k).cWrite = uicontrol(plotdata.Material(k).TBc,'Style','edit',...
      'String',num2str(plotdata.Material(k).c),...
      'Position',[78,124,40,20],...
      'callback', [plotdata.prog,' writec(',num2str(k),')']);
  
plotdata.Material(k).cSlide = uicontrol(plotdata.Material(k).TBc,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).cMin,'Max',plotdata.Material(k).cMax, ...
      'Position',[2,102,116,20], ...
      'Value', plotdata.Material(k).c,...
      'SliderStep',[0.001*1/(plotdata.Material(k).cMax-plotdata.Material(k).cMin) 0.001*1/(plotdata.Material(k).cMax-plotdata.Material(k).cMin)], ...
      'CallBack', [plotdata.prog,' slidec(',num2str(k),')']);

plotdata.Material(k).TBcExp = uicontrol('Style','radiobutton','String','<html>c<sub>z</sub>=Aexp<sup>z/B</sup>+C</html>',...
    'pos',[0,75,100,30],'parent',plotdata.Material(k).TBc);

uicontrol(plotdata.Material(k).TBc,'Style','text', ...
      'String','  A',...
      'Position',[0,63,16,14]);

plotdata.Material(k).aexpWrite = uicontrol(plotdata.Material(k).TBc,'Style','edit',...
      'String',num2str(plotdata.Material(k).aexp),...
      'Position',[21,60,35,20],...
      'callback', [plotdata.prog,' writeaexp(',num2str(k),')']);
  
plotdata.Material(k).aexpSlide = uicontrol(plotdata.Material(k).TBc,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).aexpMin,'Max',plotdata.Material(k).aexpMax, ...
      'Position',[60,57,58,20], ...
      'Value', plotdata.Material(k).aexp,...
      'SliderStep',[0.001*1/(plotdata.Material(k).aexpMax-plotdata.Material(k).aexpMin) 0.001*1/(plotdata.Material(k).aexpMax-plotdata.Material(k).aexpMin)], ...
      'CallBack', [plotdata.prog,' slideaexp(',num2str(k),')']);

uicontrol(plotdata.Material(k).TBc,'Style','text', ...
      'String','  B',...
      'Position',[0,40,16,14]);

plotdata.Material(k).bexpWrite = uicontrol(plotdata.Material(k).TBc,'Style','edit',...
      'String',num2str(plotdata.Material(k).bexp),...
      'Position',[21,37,35,20],...
      'callback', [plotdata.prog,' writebexp(',num2str(k),')']);
  
plotdata.Material(k).bexpSlide = uicontrol(plotdata.Material(k).TBc,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).bexpMin,'Max',plotdata.Material(k).bexpMax, ...
      'Position',[60,34,58,20], ...
      'Value', plotdata.Material(k).bexp,...
      'SliderStep',[1/(plotdata.Material(k).bexpMax-plotdata.Material(k).bexpMin) 1/(plotdata.Material(k).bexpMax-plotdata.Material(k).bexpMin)], ...
      'CallBack', [plotdata.prog,' slidebexp(',num2str(k),')']);

uicontrol(plotdata.Material(k).TBc,'Style','text',...
      'String','<c>',...
      'Position',[0,17,22,14]);

plotdata.Material(k).MeancexpWrite = uicontrol(plotdata.Material(k).TBc,'Style','edit',...
      'String',num2str(plotdata.Material(k).Meancexp),...
      'Position',[21,14,35,20],...
      'callback', [plotdata.prog,' writeMeancexp(',num2str(k),')']);
  
plotdata.Material(k).MeancexpSlide = uicontrol(plotdata.Material(k).TBc,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).MeancexpMin,'Max',plotdata.Material(k).MeancexpMax, ...
      'Position',[60,10,58,20], ...
      'Value', plotdata.Material(k).Meancexp,...
      'SliderStep',[0.001*1/(plotdata.Material(k).MeancexpMax-plotdata.Material(k).MeancexpMin) 0.001*1/(plotdata.Material(k).MeancexpMax-plotdata.Material(k).MeancexpMin)], ...
      'CallBack', [plotdata.prog,' slideMeancexp(',num2str(k),')']);
  
uicontrol(plotdata.Material(k).TBc,'Style','text', ...
      'String',' C=',...
      'Position',[0,0,22,14]);
  
  plotdata.Material(k).cexpWrite=uicontrol(plotdata.Material(k).TBc,'Style','text', ...
      'String',num2str(plotdata.Material(k).cexp),...
      'Position',[25,0,50,14]);
  
set(plotdata.Material(k).TBc,'SelectionChangeFcn',@togglec);

  end;
  
    % Create the orientation button group.
plotdata.TBorientation = uibuttongroup(plotdata.MainPanel.fig,'visible','on',...
    'Title','Sample orientation:',...
    'BorderType','etchedout',...
    'Position',[0.005 0.86 0.16 0.13]);
% Create two buttons in the button group.
plotdata.TB001orientation = uicontrol('Style','radiobutton','String','(001)',...
    'pos',[20,20,100,14],'parent',plotdata.TBorientation);

plotdata.TB111orientation = uicontrol('Style','radiobutton','String','(111)',...
    'pos',[20,3,100,14],'parent',plotdata.TBorientation);   
   
set(plotdata.TBorientation,'SelectionChangeFcn',@toggleorientation);

  set(plotdata.Material(1).Panel,'Title','Bottom layer:','Tag','k=1');
  set(plotdata.Material(2).Panel,'Title','Material 1:','Tag','k=2');
  set(plotdata.Material(3).Panel,'Title','Material 2:','Tag','k=3');
  set(plotdata.Material(4).Panel,'Title','Material 3:','Tag','k=4');
  set(plotdata.Material(5).Panel,'Title','Top layer:','Tag','k=5');
  
%% Create Fit display
   plotdata.FitDisplay.fig = figure('Position',plotdata.FitDisplay.centerfig,...
      'Resize','on',...
      'NumberTitle','off',...
      'Name','Fit display',...
      'Interruptible','off',...
      'Menubar','none',...
      'Color',get(0,'DefaultUIControlBackgroundColor'));
   set(plotdata.FitDisplay.fig,'toolbar','figure');
   set(plotdata.FitDisplay.fig,'menubar','figure');
   subplot('Position',[0.07,0.1,0.9,0.8]);
   
   %==X min=========================
   plotdata.chooseXMinText=uicontrol('parent',plotdata.FitDisplay.fig,'style','text',...
       'string',[plotdata.FitDisplay.haxis,' min:'],... 
       'position',[0,0,100,20]);
   plotdata.chooseXMin=uicontrol('parent',plotdata.FitDisplay.fig,'Style','edit', ...
      'String',num2str(plotdata.TThetaMin),... 
      'Position',[80,3,40,20],...
      'CallBack',[plotdata.prog,' chooseXMin'] );
  
  %==X max=========================
   plotdata.chooseXMaxText=uicontrol(plotdata.FitDisplay.fig,'style','text',...
       'string',[plotdata.FitDisplay.haxis,' max:'],... 
       'position',[520,0,100,20]);
   plotdata.chooseXMax=uicontrol(plotdata.FitDisplay.fig,'Style','edit', ...
      'String',num2str(plotdata.TThetaMax),... 
      'Position',[600,3,40,20],...
      'CallBack',[plotdata.prog,' chooseXMax'] );
  
  %==Full scale=========================
   uicontrol(plotdata.FitDisplay.fig,'style','pushbutton',...
       'string','Full scale',... 
       'position',[380,3,100,20],...
       'CallBack',[plotdata.prog,' FullScale'] );
   
%==Toggle between 2Theta and L display=========================
% Create the button group.
h = uibuttongroup('parent',plotdata.FitDisplay.fig,'visible','on','Position',[0.005,0.94,0.15/plotdata.FitDisplay.Width*660,0.06/plotdata.FitDisplay.Height*400]);
% Create two buttons in the button group.
u0 = uicontrol('Style','radiobutton','String','2Theta',...
    'pos',[0,-3,100,30],'parent',h);
u1 = uicontrol('Style','radiobutton','String','L',...
    'pos',[60,-3,100,30],'parent',h);
set(h,'SelectionChangeFcn',@toggle);
   
  %==ScalingIntensity=========================
   uicontrol(plotdata.FitDisplay.fig,'style','text',...
       'string','Scaling int.:',... 
       'position',[120,0,100,20]);
   plotdata.chooseScalingIntensity=uicontrol(plotdata.FitDisplay.fig,'Style','edit', ...
      'String',num2str(plotdata.ScalingIntensity),... 
      'Position',[200,3,40,20],...
      'CallBack',[plotdata.prog,' chooseScalingIntensity'] );

  %==RMS min=========================
   plotdata.chooseRMSTThetaMinText=uicontrol('parent',plotdata.FitDisplay.fig,'style','text',...
       'visible','off',...
       'string',' RMS TThetamin:',... 
       'position',[plotdata.FitDisplay.Width/2,plotdata.FitDisplay.Height-27,100,20]);
   plotdata.chooseRMSTThetaMin=uicontrol('parent',plotdata.FitDisplay.fig,'Style','edit', ...
       'visible','off',...
      'String',num2str(plotdata.RMS.TThetaMin),... 
      'Position',[plotdata.FitDisplay.Width/2+90,plotdata.FitDisplay.Height-24,40,20],...
      'CallBack',[plotdata.prog,' chooseRMSMinMax'] );
   plotdata.chooseRMSLMinText=uicontrol('parent',plotdata.FitDisplay.fig,'style','text',...
       'visible','off',...
       'string',' RMS Lmin:',... 
       'position',[plotdata.FitDisplay.Width/2,plotdata.FitDisplay.Height-27,100,20]);
   plotdata.chooseRMSLMin=uicontrol('parent',plotdata.FitDisplay.fig,'Style','edit', ...
       'visible','off',...
      'String',num2str(plotdata.RMS.LMin),... 
      'Position',[plotdata.FitDisplay.Width/2+90,plotdata.FitDisplay.Height-24,40,20],...
      'CallBack',[plotdata.prog,' chooseRMSMinMax'] );
  %==RMS max=========================
   plotdata.chooseRMSTThetaMaxText=uicontrol('parent',plotdata.FitDisplay.fig,'style','text',...
       'string',' RMS TThetamax:',... 
       'visible','off',...
       'position',[plotdata.FitDisplay.Width/2+130,plotdata.FitDisplay.Height-27,100,20]);
   plotdata.chooseRMSTThetaMax=uicontrol('parent',plotdata.FitDisplay.fig,'Style','edit', ...
      'String',num2str(plotdata.RMS.TThetaMax),... 
      'visible','off',...
      'Position',[plotdata.FitDisplay.Width/2+220,plotdata.FitDisplay.Height-24,40,20],...
      'CallBack',[plotdata.prog,' chooseRMSMinMax'] );
   plotdata.chooseRMSLMaxText=uicontrol('parent',plotdata.FitDisplay.fig,'style','text',...
       'string',' RMS Lmax:',... 
       'visible','off',...
       'position',[plotdata.FitDisplay.Width/2+130,plotdata.FitDisplay.Height-27,100,20]);
   plotdata.chooseRMSLMax=uicontrol('parent',plotdata.FitDisplay.fig,'Style','edit', ...
      'String',num2str(plotdata.RMS.LMax),... 
      'visible','off',...
      'Position',[plotdata.FitDisplay.Width/2+220,plotdata.FitDisplay.Height-24,40,20],...
      'CallBack',[plotdata.prog,' chooseRMSMinMax'] );
     
   I=plotdata.Substrate.g.*conj(plotdata.Substrate.g);
   plotdata.fit.y=I/max(I);
   semilogy(plotdata.fit.x.TTheta,plotdata.fit.y*plotdata.ScalingIntensity,'b');    
   axis([plotdata.TThetaMin plotdata.TThetaMax 10^-9 1]);
   xlabel('2Theta (°)'); ylabel('Intensity (arbitrary units)');
  %% Create c display
   plotdata.cDisplay.fig = figure('Position',plotdata.cDisplay.centerfig,...
      'Resize','on',...
      'NumberTitle','off',...
      'Name','c display',...
      'Interruptible','off',...
      'Menubar','none',...
      'Color',get(0,'DefaultUIControlBackgroundColor'));
   set(plotdata.cDisplay.fig,'toolbar','figure');
   set(plotdata.cDisplay.fig,'menubar','figure');
   
%%
case 'quit'
   quit_reply = questdlg('Really quit this nice fitting progam?');
   if strcmp(quit_reply,'Yes')
      close all;
   end
   close all;
   clear all;
%%   
case 'export' 
   if strcmp(plotdata.FitDisplay.haxis,'L')
      x=plotdata.fit.x.L; 
   elseif strcmp(plotdata.FitDisplay.haxis,'2Theta')
      x=plotdata.fit.x.TTheta; 
   else
    warning('Error in export')
   end
   y=plotdata.fit.y*plotdata.ScalingIntensity;
   if exist('InteractiveXRDFit.mat', 'file')
       load('InteractiveXRDFit.mat','startfolder');
       startfolder = strcat(uigetdir(startfolder),filesep);
   else
       startfolder = strcat(uigetdir(),filesep);
   end
   
   if ~strcmp(startfolder,'/')
    save('InteractiveXRDFit.mat','startfolder');
    name = startfolder;
    if ~strcmp(plotdata.measure.filename,'none')
      name=[name,'fit ',strtok(plotdata.measure.filename,'.'),' with '];
                %name=strcat(name,'fit ',strtok(plotdata.measure.filename,'.'),' with ');

    end   
        name=[name,plotdata.orientation];
%                name=strcat(name,plotdata.orientation);

    if ~strcmp(plotdata.Substrate.Type,'none')
        name=[name,plotdata.Substrate.Type,' substrate '];
%        name=strcat(name,plotdata.Substrate.Type,' substrate ');

    end
    if ~strcmp(plotdata.Material(1).Type,'none')
        name=[name,char(plotdata.Material(1).Type),'_',char(num2str(plotdata.Material(1).N)),'x'];
    %        name=strcat(name,plotdata.Material(1).Type,'_',num2str(plotdata.Material(1).N),'x');

        if strcmp(plotdata.orientation,'(001)')
            name=[name,char(num2str(plotdata.Material(1).c)),'A '];
            %            name=strcat(name,num2str(plotdata.Material(1).c),'A ');

        elseif strcmp(plotdata.orientation,'(111)')
            name=[name,char(num2str(plotdata.Material(1).d)),'A '];
 %            name=strcat(name,num2str(plotdata.Material(1).d),'A ');

        else warning('Error in plotdata.orientation in Export case');
        end;
    end
    if plotdata.Repetition.N~=1
        name=[name,char(num2str(plotdata.Repetition.N)),'x( '];
        %        name=strcat(name,num2str(plotdata.Repetition.N),'x( ');

    end   
    for k=2:4,
        if ~strcmp(plotdata.Material(k).Type,'none')
            name=[name,char(plotdata.Material(k).Type),'_',char(num2str(plotdata.Material(k).N)),'x'];
            %            name=strcat(name,plotdata.Material(k).Type,'_',num2str(plotdata.Material(k).N),'x');

            if strcmp(plotdata.orientation,'(001)')
                name=[name,char(num2str(plotdata.Material(k).c)),'A '];
                %                name=strcat(name,num2str(plotdata.Material(k).c),'A ');

            elseif strcmp(plotdata.orientation,'(111)')
                name=[name,char(num2str(plotdata.Material(k).d)),'A '];
                %                name=strcat(name,num2str(plotdata.Material(k).c),'A ');

            else warning('Error in plotdata.orientation in Export case');
            end;    
        end;
    end;
    if plotdata.Repetition.N~=1
        name=[name,' ) '];
    end    
    if ~strcmp(plotdata.Material(5).Type,'none')
        name=[name,char(plotdata.Material(5).Type),'_',char(num2str(plotdata.Material(5).N)),'x'];
        if strcmp(plotdata.orientation,'(001)')
            name=[name,char(num2str(plotdata.Material(5).c)),'A '];
        elseif strcmp(plotdata.orientation,'(111)')
            name=[name,char(num2str(plotdata.Material(5).d)),'A '];
        else warning('Error in plotdata.orientation in Export case');
        end;    
    end
    
    fitname=[name,'.csv'];

    export_reply = questdlg(['Your fit will be saved in ',fitname,'. Do you want to choose a different name for your file?']);
    if strcmp(export_reply,'Yes')
      prompt={'Please indicate a name for your file:'};
      dlg_title='File name';
      num_lines=1;
      def={fitname};
      answer=inputdlg(prompt,dlg_title,num_lines,def);
      fitname=answer{1};
    end
    
    csvwrite(fitname,[x y']);
    
   end
%%   
case {'slidePolarization(1)','slidePolarization(2)','slidePolarization(3)','slidePolarization(4)','slidePolarization(5)'}
   k=str2num(entry(19));
   plotdata.Material(k).Polarization=get(plotdata.Material(k).PolarizationSlide,'Value');
   set(plotdata.Material(k).PolarizationWrite,'String',plotdata.Material(k).Polarization);
   set(plotdata.Material(k).TBPolarization,'SelectedObject',plotdata.Material(k).TBPolarizationYes);
   ProgcNSimu
%%     
case {'writePolarization(1)','writePolarization(2)','writePolarization(3)','writePolarization(4)','writePolarization(5)'}
   k=str2num(entry(19));
   plotdata.Material(k).Polarization=str2num(get(plotdata.Material(k).PolarizationWrite,'string'));
   if plotdata.Material(k).Polarization>plotdata.Material(k).PolarizationMax,
       plotdata.Material(k).Polarization=plotdata.Material(k).PolarizationMax;
       set(plotdata.Material(k).PolarizationWrite,'String',plotdata.Material(k).Polarization);
   elseif plotdata.Material(k).Polarization<plotdata.Material(k).PolarizationMin,
       plotdata.Material(k).Polarization=plotdata.Material(k).PolarizationMin;
       set(plotdata.Material(k).PolarizationWrite,'String',plotdata.Material(k).Polarization);
   end;
   set(plotdata.Material(k).PolarizationSlide,'value',plotdata.Material(k).Polarization);
   set(plotdata.Material(k).TBPolarization,'SelectedObject',plotdata.Material(k).TBPolarizationYes);
   ProgcNSimu
%%   
case {'slidexPST(1)','slidexPST(2)','slidexPST(3)','slidexPST(4)','slidePST(5)'}
   k=str2num(entry(11));
   plotdata.Material(k).xPST=get(plotdata.Material(k).xPSTSlide,'Value');
   set(plotdata.Material(k).xPSTWrite,'String',plotdata.Material(k).xPST);
   ProgcNSimu
%%     
case {'writexPST(1)','writexPST(2)','writexPST(3)','writexPST(4)','writexPST(5)'}
   k=str2num(entry(11));
   plotdata.Material(k).xPST=str2num(get(plotdata.Material(k).xPSTWrite,'string'));
   if plotdata.Material(k).xPST>plotdata.Material(k).xPSTMax,
       plotdata.Material(k).xPST=plotdata.Material(k).xPSTMax;
       set(plotdata.Material(k).xPSTWrite,'String',plotdata.Material(k).xPST);
   elseif plotdata.Material(k).xPST<plotdata.Material(k).xPSTMin,
       plotdata.Material(k).xPST=plotdata.Material(k).xPSTMin;
       set(plotdata.Material(k).xPSTWrite,'String',plotdata.Material(k).xPST);
   end;
   set(plotdata.Material(k).xPSTSlide,'value',plotdata.Material(k).xPST);
   ProgcNSimu
%%   
case {'slidexPZT(1)','slidexPZT(2)','slidexPZT(3)','slidexPZT(4)','slidePZT(5)'}
   k=str2num(entry(11));
   plotdata.Material(k).xPZT=get(plotdata.Material(k).xPZTSlide,'Value');
   set(plotdata.Material(k).xPZTWrite,'String',plotdata.Material(k).xPZT);
   ProgcNSimu
%%     
case {'writexPZT(1)','writexPZT(2)','writexPZT(3)','writexPZT(4)','writexPZT(5)'}
   k=str2num(entry(11));
   plotdata.Material(k).xPZT=str2num(get(plotdata.Material(k).xPZTWrite,'string'));
   if plotdata.Material(k).xPZT>plotdata.Material(k).xPZTMax,
       plotdata.Material(k).xPZT=plotdata.Material(k).xPZTMax;
       set(plotdata.Material(k).xPZTWrite,'String',plotdata.Material(k).xPZT);
   elseif plotdata.Material(k).xPZT<plotdata.Material(k).xPZTMin,
       plotdata.Material(k).xPZT=plotdata.Material(k).xPZTMin;
       set(plotdata.Material(k).xPZTWrite,'String',plotdata.Material(k).xPZT);
   end;
   set(plotdata.Material(k).xPZTSlide,'value',plotdata.Material(k).xPZT);
   ProgcNSimu
%% 
case {'slideN(1)','slideN(2)','slideN(3)','slideN(4)','slideN(5)'}
   k=str2num(entry(8));
   plotdata.Material(k).N=get(plotdata.Material(k).NSlide,'Value');
   set(plotdata.Material(k).NWrite,'String',plotdata.Material(k).N);
   ProgcNSimu
%%
case {'writeN(1)','writeN(2)','writeN(3)','writeN(4)','writeN(5)'}
   k=str2num(entry(8));
   plotdata.Material(k).N=str2num(get(plotdata.Material(k).NWrite,'string'));
   if plotdata.Material(k).N>plotdata.Material(k).NMax,
       plotdata.Material(k).N=plotdata.Material(k).NMax;
       set(plotdata.Material(k).NWrite,'String',plotdata.Material(k).N);
   elseif plotdata.Material(k).N<plotdata.Material(k).NMin,
       plotdata.Material(k).N=plotdata.Material(k).NMin;
       set(plotdata.Material(k).NWrite,'String',plotdata.Material(k).N);
   end;
   set(plotdata.Material(k).NSlide,'value',plotdata.Material(k).N);
   ProgcNSimu
%%   
case {'slidec(1)','slidec(2)','slidec(3)','slidec(4)','slidec(5)'}
   k=str2num(entry(8)); 
   plotdata.Material(k).c=get(plotdata.Material(k).cSlide,'Value');
   set(plotdata.Material(k).cWrite,'String',plotdata.Material(k).c);
   set(plotdata.Material(k).TBc,'SelectedObject',plotdata.Material(k).TBcConst);
   plotdata.Material(k).cDistribution='constant';
   ProgcNSimu
%%  
case {'writec(1)','writec(2)','writec(3)','writec(4)','writec(5)'}
   k=str2num(entry(8)); 
   plotdata.Material(k).c=str2num(get(plotdata.Material(k).cWrite,'string'));
   if plotdata.Material(k).c>plotdata.Material(k).cMax,
       plotdata.Material(k).c=plotdata.Material(k).cMax;
       set(plotdata.Material(k).cWrite,'String',plotdata.Material(k).c);
   elseif plotdata.Material(k).c<plotdata.Material(k).cMin,
       plotdata.Material(k).c=plotdata.Material(k).cMin;
       set(plotdata.Material(k).cWrite,'String',plotdata.Material(k).c);
   end;
   set(plotdata.Material(k).cSlide,'value',plotdata.Material(k).c);
   set(plotdata.Material(k).TBc,'SelectedObject',plotdata.Material(k).TBcConst);
   plotdata.Material(k).cDistribution='constant';
   ProgcNSimu
%%
case {'slideaexp(1)','slideaexp(2)','slideaexp(3)','slideaexp(4)','slideaexp(5)'}
   k=str2num(entry(11)); 
   plotdata.Material(k).aexp=get(plotdata.Material(k).aexpSlide,'Value');
   set(plotdata.Material(k).aexpWrite,'String',plotdata.Material(k).aexp);
   set(plotdata.Material(k).TBc,'SelectedObject',plotdata.Material(k).TBcExp);
   plotdata.Material(k).cDistribution='exp';   
   ProgcNSimu
   k=str2num(entry(11)); 
   set(plotdata.Material(k).cexpWrite,'String',num2str(plotdata.Material(k).cexp));
%%   
case {'writeaexp(1)','writeaexp(2)','writeaexp(3)','writeaexp(4)','writeaexp(5)'}
   k=str2num(entry(11));
   plotdata.Material(k).aexp=str2num(get(plotdata.Material(k).aexpWrite,'string'));
   if plotdata.Material(k).aexp>plotdata.Material(k).aexpMax,
       plotdata.Material(k).aexp=plotdata.Material(k).aexpMax;
       set(plotdata.Material(k).aexpWrite,'String',plotdata.Material(k).aexp);
   elseif plotdata.Material(k).aexp<plotdata.Material(k).aexpMin,
       plotdata.Material(k).aexp=plotdata.Material(k).aexpMin;
       set(plotdata.Material(k).aexpWrite,'String',plotdata.Material(k).aexp);
   end;
   set(plotdata.Material(k).aexpSlide,'value',plotdata.Material(k).aexp);
   set(plotdata.Material(k).TBc,'SelectedObject',plotdata.Material(k).TBcExp);
   plotdata.Material(k).cDistribution='exp';
   ProgcNSimu
   k=str2num(entry(11));
   set(plotdata.Material(k).cexpWrite,'String',num2str(plotdata.Material(k).cexp));
%%   
case {'slidebexp(1)','slidebexp(2)','slidebexp(3)','slidebexp(4)','slidebexp(5)'}
   k=str2num(entry(11));
   plotdata.Material(k).bexp=get(plotdata.Material(k).bexpSlide,'Value');
   set(plotdata.Material(k).bexpWrite,'String',plotdata.Material(k).bexp);
   set(plotdata.Material(k).TBc,'SelectedObject',plotdata.Material(k).TBcExp);
   plotdata.Material(k).cDistribution='exp';
   ProgcNSimu
   k=str2num(entry(11));
   set(plotdata.Material(k).cexpWrite,'String',num2str(plotdata.Material(k).cexp));
%%   
case {'writebexp(1)','writebexp(2)','writebexp(3)','writebexp(4)','writebexp(5)'}
   k=str2num(entry(11));
   plotdata.Material(k).bexp=str2num(get(plotdata.Material(k).bexpWrite,'string'));
   if plotdata.Material(k).bexp>plotdata.Material(k).bexpMax,
       plotdata.Material(k).bexp=plotdata.Material(k).bexpMax;
       set(plotdata.Material(k).bexpWrite,'String',plotdata.Material(k).bexp);
   elseif plotdata.Material(k).bexp<plotdata.Material(k).bexpMin,
       plotdata.Material(k).bexp=plotdata.Material(k).bexpMin;
       set(plotdata.Material(k).bexpWrite,'String',plotdata.Material(k).bexp);
   end;
   set(plotdata.Material(k).bexpSlide,'value',plotdata.Material(k).bexp);
   set(plotdata.Material(k).TBc,'SelectedObject',plotdata.Material(k).TBcExp);
   plotdata.Material(k).cDistribution='exp';
   ProgcNSimu
   k=str2num(entry(11));
   set(plotdata.Material(k).cexpWrite,'String',num2str(plotdata.Material(k).cexp));
%%   
case {'slideMeancexp(1)','slideMeancexp(2)','slideMeancexp(3)','slideMeancexp(4)','slideMeancexp(5)'}
   k=str2num(entry(15));
   plotdata.Material(k).Meancexp=get(plotdata.Material(k).MeancexpSlide,'Value');
   set(plotdata.Material(k).MeancexpWrite,'String',plotdata.Material(k).Meancexp);
   set(plotdata.Material(k).TBc,'SelectedObject',plotdata.Material(k).TBcExp);
   plotdata.Material(k).cDistribution='exp';   
   ProgcNSimu
   k=str2num(entry(15));
   set(plotdata.Material(k).cexpWrite,'String',num2str(plotdata.Material(k).cexp));
%%   
case {'writeMeancexp(1)','writeMeancexp(2)','writeMeancexp(3)','writeMeancexp(4)','writeMeancexp(5)'}
   k=str2num(entry(15));
   plotdata.Material(k).Meancexp=str2num(get(plotdata.Material(k).MeancexpWrite,'string'));
   if plotdata.Material(k).Meancexp>plotdata.Material(k).MeancexpMax,
       plotdata.Material(k).Meancexp=plotdata.Material(k).MeancexpMax;
       set(plotdata.Material(k).MeancexpWrite,'String',plotdata.Material(k).Meancexp);
   elseif plotdata.Material(k).Meancexp<plotdata.Material(k).MeancexpMin,
       plotdata.Material(k).Meancexp=plotdata.Material(k).MeancexpMin;
       set(plotdata.Material(k).MeancexpWrite,'String',plotdata.Material(k).Meancexp);
   end;
   set(plotdata.Material(k).MeancexpSlide,'value',plotdata.Material(k).Meancexp);
   set(plotdata.Material(k).TBc,'SelectedObject',plotdata.Material(k).TBcExp);
   plotdata.Material(k).cDistribution='exp';
   ProgcNSimu
   k=str2num(entry(15));
   set(plotdata.Material(k).cexpWrite,'String',num2str(plotdata.Material(k).cexp));
%%
case {'slided(1)','slided(2)','slided(3)','slided(4)','slided(5)'}
   k=str2num(entry(8));
   plotdata.Material(k).d=get(plotdata.Material(k).dSlide,'Value');
   set(plotdata.Material(k).dWrite,'String',plotdata.Material(k).d);
   ProgcNSimu
%%   
case {'writed(1)','writed(2)','writed(3)','writed(4)','writed(5)'}
   k=str2num(entry(8));
   plotdata.Material(k).d=str2num(get(plotdata.Material(k).dWrite,'string'));
   if plotdata.Material(k).d>plotdata.Material(k).dMax,
       plotdata.Material(k).d=plotdata.Material(k).dMax;
       set(plotdata.Material(k).dWrite,'String',plotdata.Material(k).d);
   elseif plotdata.Material(k).d<plotdata.Material(k).dMin,
       plotdata.Material(k).d=plotdata.Material(k).dMin;
       set(plotdata.Material(k).dWrite,'String',plotdata.Material(k).d);
   end;
   set(plotdata.Material(k).dSlide,'value',plotdata.Material(k).d);
   ProgcNSimu
%%   
case 'slideRepetition'
   plotdata.Repetition.N=get(plotdata.RepetitionSlide,'Value');
   set(plotdata.RepetitionWrite,'String',plotdata.Repetition.N);
   ProgcNSimu
%%    
case 'writeRepetition'
   plotdata.Repetition.N=str2num(get(plotdata.RepetitionWrite,'string'));
   if plotdata.Repetition.N>plotdata.MaxRepetition,
       plotdata.Repetition.N=plotdata.MaxRepetition;
       set(plotdata.RepetitionWrite,'String',plotdata.Repetition.N);
   elseif plotdata.Repetition.N<plotdata.MinRepetition,
       plotdata.Repetition.N=plotdata.MinRepetition;
       set(plotdata.RepetitionWrite,'String',plotdata.Repetition.N);
   end;
   set(plotdata.RepetitionSlide,'value',plotdata.Repetition.N);
   ProgcNSimu
%%
case 'LoadDataFilename'   
   if exist('InteractiveXRDFitload.mat', 'file')
       load('InteractiveXRDFitload.mat','loadfolder'); 
       if loadfolder==0;
           loadfolder=[];
       end;
       [plotdata.measure.filename, loadfolder] = uigetfile({'*.csv';'*.*'},'Select the data you wand to fit (.csv)',loadfolder);
   else
       [plotdata.measure.filename, loadfolder] = uigetfile({'*.csv';'*.*'},'Select the data you wand to fit (.csv)');
   end
   if loadfolder~=0;
    save('InteractiveXRDFitload.mat','loadfolder');
    plotdata.prog=[loadfolder,plotdata.measure.filename];
    set(plotdata.filenameWrite,'String',plotdata.prog);
    %open and load data
    fid=fopen([plotdata.prog],'r');
    data=[];
    if (fid==-1)
        msgbox(['ERROR: file not found or could not be opened for read. check file ' [plotdata.prog]]); 
        return;
    end 
    % skip first line of comments 
    sLine=fgetl(fid); 
    % read the data 
    while 1 
        sLine=fgetl(fid); 
        if (~ischar(sLine)) 
            % end of file 
            break; 
        end
        % add line to table
        data=[data;str2num(sLine)];
    end 
    fclose(fid); 
    plotdata.measure.x=data(:,1);
    plotdata.measure.y=data(:,2);
    plotdata.measure.y=plotdata.measure.y/max(plotdata.measure.y);
    
    PlotFitAndData
    figure(plotdata.FitDisplay.fig);
    if strcmp(plotdata.FitDisplay.haxis,'L')
        helpdlg('Displays your measurement in red. The two vertical green lines delimitate the region over which the RMS (root mean square of the pair wise differences of the fit and measurement) is calculated. You can modify the limits by changing the limiting values on the Fit display panel.');
        set(plotdata.chooseRMSLMinText,'visible','on');
        set(plotdata.chooseRMSLMin,'visible','on');
        set(plotdata.chooseRMSLMaxText,'visible','on');
        set(plotdata.chooseRMSLMax,'visible','on');
        set(plotdata.chooseRMSTThetaMinText,'visible','off');
        set(plotdata.chooseRMSTThetaMin,'visible','off');
        set(plotdata.chooseRMSTThetaMaxText,'visible','off');
        set(plotdata.chooseRMSTThetaMax,'visible','off');
    elseif strcmp(plotdata.FitDisplay.haxis,'2Theta')
        helpdlg('Displays your measurement in red. The two vertical green lines delimitate the region over which the RMS (root mean square of the pair wise differences of the fit and measurement) is calculated. You can modify the limits by changing the limiting values on the Fit display panel.');
        set(plotdata.chooseRMSTThetaMinText,'visible','on');
        set(plotdata.chooseRMSTThetaMin,'visible','on');
        set(plotdata.chooseRMSTThetaMaxText,'visible','on');
        set(plotdata.chooseRMSTThetaMax,'visible','on');
        set(plotdata.chooseRMSLMinText,'visible','off');
        set(plotdata.chooseRMSLMin,'visible','off');
        set(plotdata.chooseRMSLMaxText,'visible','off');
        set(plotdata.chooseRMSLMax,'visible','off');
    else warning('Error in plotdata.FitDisplay.haxis')
    end;
   end;
          
case 'chooseSubstrate'
   str=get(plotdata.Substrate.choose,'String');
   entry=get(plotdata.Substrate.choose,'value');
   SubstrateType=strsplit(str(entry,:));
   plotdata.Substrate.Type=SubstrateType(1);
   Substrate
   ProgcNSimu
%%   
case {'chooseLayerType(1)','chooseLayerType(2)','chooseLayerType(3)','chooseLayerType(4)','chooseLayerType(5)'}
   k=str2num(entry(17));
   str=get(plotdata.Material(k).choose,'String');
   entry=get(plotdata.Material(k).choose,'Value');
   MaterialType=strsplit(str(entry,:));
   plotdata.Material(k).Type=MaterialType(1);
   plotdata.Material(k).Polarization=0;
   plotdata.Material(k).xPST=0;
   plotdata.Material(k).xPZT=0;
   set(plotdata.Material(k).xPSTpanel,'Visible','off');
   set(plotdata.Material(k).xPZTpanel,'Visible','off');
   set(plotdata.Material(k).TBPolarization,'Visible','off');
   if strcmp(plotdata.Material(k).Type,'none'),
       plotdata.Material(k).N=0; plotdata.Material(k).c=0; plotdata.Material(k).d=0;
       set(plotdata.Material(k).choose,'value',1);
   elseif any(strcmp(plotdata.Material(k).Type,{'AlO2','BaO','CaCuO2','La2CuO4','LaO','MnO','MnO2','NiO2','NdO','PbO','PrBa2Cu3O7','RuO2','SrO','SrO2','TiO2','VO2','YBa2Cu3O7','ZrO2'})),
       if strcmp(plotdata.orientation,'(111)'),
           msgbox({'This program doesn''t allow you to calculate ' char(plotdata.Material(k).Type) '-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.'});
           plotdata.Material(k).Type='none';
       end;
   elseif strcmp(plotdata.Material(k).Type,'PbTiO3'),
       if strcmp(plotdata.orientation,'(001)'),
           set(plotdata.Material(k).TBPolarization,'Visible','on');
       end;
   elseif strcmp(plotdata.Material(k).Type,'(Pb_x,Sr_{1-x})TiO3'),
       set(plotdata.Material(k).xPSTpanel,'Visible','on');
   elseif strcmp(plotdata.Material(k).Type,'Pb(Zr_x,Ti_{1-x})O3'),
       set(plotdata.Material(k).xPZTpanel,'Visible','on');
   end;
   ProgcNSimu
%% 
case 'chooseXMin'   
   value=str2num(get(plotdata.chooseXMin,'string'));
   
   if strcmp(plotdata.FitDisplay.haxis,'L')
      if value<plotdata.LMax,
          plotdata.LMin=value;
          axis([plotdata.LMin plotdata.LMax 10^-9 1]);
          set(plotdata.chooseXMin,'String',plotdata.LMin);
      end;        
       
   elseif strcmp(plotdata.FitDisplay.haxis,'2Theta')
      if value<plotdata.TThetaMax,
          plotdata.TThetaMin=value;
          axis([plotdata.TThetaMin plotdata.TThetaMax 10^-9 1]);
          set(plotdata.chooseXMin,'String',plotdata.TThetaMin);
      end; 
   else warning('Error in chooseXMin')
   end;
%%      
case 'chooseXMax'   
   value=str2num(get(plotdata.chooseXMax,'string'));

   if strcmp(plotdata.FitDisplay.haxis,'L')
      if value>plotdata.LMin,
          plotdata.LMax=value;
          axis([plotdata.LMin plotdata.LMax 10^-9 1]);
          set(plotdata.chooseXMax,'String',plotdata.LMax);
      end;        
       
   elseif strcmp(plotdata.FitDisplay.haxis,'2Theta')
      if value>plotdata.TThetaMin,
          plotdata.TThetaMax=value;
          axis([plotdata.TThetaMin plotdata.TThetaMax 10^-9 1]);
          set(plotdata.chooseXMax,'String',plotdata.TThetaMax);
      end; 
   else warning('Error in chooseXMax')
   end;
%%
    case 'FullScale'
   if strcmp(plotdata.FitDisplay.haxis,'L')
      plotdata.LMin=0;
      set(plotdata.chooseXMin,'String',plotdata.LMin);
      plotdata.LMax=2*plotdata.Substrate.c/plotdata.lambda;
      set(plotdata.chooseXMax,'String',plotdata.LMax);
      axis([plotdata.LMin plotdata.LMax 10^-9 1]);
   elseif strcmp(plotdata.FitDisplay.haxis,'2Theta')
      plotdata.TThetaMin=0;
      set(plotdata.chooseXMin,'String',plotdata.TThetaMin);
      plotdata.TThetaMax=180;
      set(plotdata.chooseXMax,'String',plotdata.TThetaMax);
      axis([plotdata.TThetaMin plotdata.TThetaMax 10^-9 1]);
   else warning('Error in FullScale')
   end;
%%
case 'chooseScalingIntensity'
   plotdata.ScalingIntensity=str2num(get(plotdata.chooseScalingIntensity,'string'));
   PlotFitAndData
%%   
   case 'chooseRMSMinMax'   
   delete(plotdata.hVerticalLines);
   if strcmp(plotdata.FitDisplay.haxis,'L')
          plotdata.RMS.LMin=str2num(get(plotdata.chooseRMSLMin,'string'));
          plotdata.RMS.LMax=str2num(get(plotdata.chooseRMSLMax,'string'));
          plotdata.hVerticalLines = line([plotdata.RMS.LMin,plotdata.RMS.LMin],get(plotdata.hAxes,'Ylim'),'Color','green');
          plotdata.hVerticalLines = [plotdata.hVerticalLines line([plotdata.RMS.LMax,plotdata.RMS.LMax],get(plotdata.hAxes,'Ylim'),'Color','green')];   
   elseif strcmp(plotdata.FitDisplay.haxis,'2Theta')
          plotdata.RMS.TThetaMin=str2num(get(plotdata.chooseRMSTThetaMin,'string'));
          plotdata.RMS.TThetaMax=str2num(get(plotdata.chooseRMSTThetaMax,'string'));
          plotdata.hVerticalLines = line([plotdata.RMS.TThetaMin,plotdata.RMS.TThetaMin],get(plotdata.hAxes,'Ylim'),'Color','green');
          plotdata.hVerticalLines = [plotdata.hVerticalLines line([plotdata.RMS.TThetaMax,plotdata.RMS.TThetaMax],get(plotdata.hAxes,'Ylim'),'Color','green')];
   else warning('Error in plotdata.FitDisplay.haxis value');
      end; 
   RMS       
   
end
%===A utility to center the window on the screen============
function pos = centerfig(width,height)

% Find the screen size in pixels
screen_s = get(0,'ScreenSize');
pos = [screen_s(3)/2 - width/2, screen_s(4)/2 - height/2, width, height];

function toggleorientation(source,eventdata)
global plotdata;
if strcmp(get(get(source,'SelectedObject'),'String'),'(001)')
    plotdata.orientation='(001)';
    for k=1:5,    
    set(plotdata.Material(k).TBc,'visible','on');
    set(plotdata.Material(k).choosed,'visible','off');
    if strcmp(plotdata.Material(k).Type,'PbTiO3'),
        set(plotdata.Material(k).TBPolarization,'Visible','on');
    end;
    end;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'(111)')
    plotdata.orientation='(111)';
    for k=1:5,
        set(plotdata.Material(k).TBc,'visible','off');
        set(plotdata.Material(k).choosed,'visible','on');
        plotdata.Material(k).cDistribution='constant';
        set(plotdata.Material(k).TBPolarization,'Visible','off');
    end;
else
    warning('Error in Toggle function (orientation)')
end
Substrate
ProgcNSimu

function toggle(source,eventdata)
global plotdata;
if strcmp(get(get(source,'SelectedObject'),'String'),'L')
    plotdata.FitDisplay.haxis='L';
    plotdata.LMin=2*plotdata.Substrate.c*sin(plotdata.TThetaMin/2*pi/180)/plotdata.lambda;
    plotdata.LMax=2*plotdata.Substrate.c*sin(plotdata.TThetaMax/2*pi/180)/plotdata.lambda;
    set(plotdata.chooseXMin,'String',plotdata.LMin);
    set(plotdata.chooseXMax,'String',plotdata.LMax);
    set(plotdata.chooseRMSTThetaMin,'visible','off');
    set(plotdata.chooseRMSTThetaMinText,'visible','off');
    set(plotdata.chooseRMSTThetaMax,'visible','off');
    set(plotdata.chooseRMSTThetaMaxText,'visible','off');
    set(plotdata.chooseRMSLMin,'visible','on','String',num2str(plotdata.RMS.LMin));
    set(plotdata.chooseRMSLMinText,'visible','on');
    set(plotdata.chooseRMSLMax,'visible','on','String',num2str(plotdata.RMS.LMax));
    set(plotdata.chooseRMSLMaxText,'visible','on');
    
elseif strcmp(get(get(source,'SelectedObject'),'String'),'2Theta')
    plotdata.FitDisplay.haxis='2Theta';   
    set(plotdata.chooseXMin,'String',plotdata.TThetaMin);
    set(plotdata.chooseXMax,'String',plotdata.TThetaMax);
    set(plotdata.chooseRMSTThetaMin,'visible','on','String',num2str(plotdata.RMS.TThetaMin));
    set(plotdata.chooseRMSTThetaMinText,'visible','on');
    set(plotdata.chooseRMSTThetaMax,'visible','on','String',num2str(plotdata.RMS.TThetaMax));
    set(plotdata.chooseRMSTThetaMaxText,'visible','on');
    set(plotdata.chooseRMSLMin,'visible','off');
    set(plotdata.chooseRMSLMinText,'visible','off');
    set(plotdata.chooseRMSLMax,'visible','off');
    set(plotdata.chooseRMSLMaxText,'visible','off');
else
    warning('Error in Toggle function')
end
set(plotdata.chooseXMinText,'String',[plotdata.FitDisplay.haxis,' min:']);
set(plotdata.chooseXMaxText,'String',[plotdata.FitDisplay.haxis,' max:']);
PlotFitAndData
%%
function togglePolarization(source,eventdata)
global plotdata;
k=str2num(source.Parent.Tag(3));
if strcmp(get(get(source,'SelectedObject'),'String'),'Yes')
elseif strcmp(get(get(source,'SelectedObject'),'String'),'No')
    plotdata.Material(k).Polarization=0;
    set(plotdata.Material(k).PolarizationWrite,'String',plotdata.Material(k).Polarization);
    set(plotdata.Material(k).PolarizationSlide,'value',plotdata.Material(k).Polarization);    
else
    warning('Error in toggle function for Polarization')
end
ProgcNSimu;
%%
function togglec(source,eventdata)
global plotdata;
k=str2num(source.Parent.Tag(3));
if strcmp(get(get(source,'SelectedObject'),'String'),'c=constant:')
    plotdata.Material(k).cDistribution='constant';
elseif strcmp(get(get(source,'SelectedObject'),'String'),'<html>c<sub>z</sub>=Aexp<sup>z/B</sup>+C</html>')
    plotdata.Material(k).cDistribution='exp';
else
    warning('Error in toggle function for cDistribution')
end
ProgcNSimu;