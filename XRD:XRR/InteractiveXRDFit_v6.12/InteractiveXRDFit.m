% Created by Celine Lichtensteiger
% Journal of Applied Crystallography - Computer programs - 51(6)p.1745 (2018)
% https://doi.org/10.1107/S1600576718012840
% This is the heart of the program
% It initializes all the variables and creates the Graphical User Interface
%*********************************

function InteractiveXRDFit(entry)

global plotdata;

if nargin == 0
   entry = 'makeGraph';
end

switch entry

case 'makeGraph'
%%Default values
% Spatial units are Angströms
   
   % lambda
   plotdata.lambda=1.540562;%;1.54059;
   % Substrate
   load('Substrates.mat');
   SubstrateTable=struct2table(SubstrateData);
   plotdata.Substrate.Type='none';
   condition=find(strcmp(SubstrateTable.Name,plotdata.Substrate.Type));
   plotdata.Substrate.Orientation='';
   plotdata.Substrate.d=0;
   plotdata.Substrate.g=0;
   plotdata.Substrate.N=7000;
   plotdata.Substrate.NMax=100000;       
   plotdata.Substrate.NMin=0;
   plotdata.TotalThickness=100000;%Units: A
   % Bottom Layer
   plotdata.Material(1).Aexp=0;
   plotdata.Material(1).Bexp=1;
   plotdata.Material(1).Cexp=4;
   plotdata.Material(1).d=4;
   plotdata.Material(1).dBorder=2;
   plotdata.Material(1).dCenter=4;
   plotdata.Material(1).dDistribution='constant';
   %plotdata.Material(1).d=2.309;%4/sqrt(3);
   plotdata.Material(1).m=1; %modulus of the Jacobi Elliptic function
   plotdata.Material(1).Meandexp=4;
   plotdata.Material(1).N=0;
   plotdata.Material(1).Orientation='none';
   plotdata.Material(1).Polarization=0;
   plotdata.Material(1).Type='none';
   plotdata.Material(1).xBST=0;
   plotdata.Material(1).xCSA=3;
   plotdata.Material(1).xHZO=0;
   plotdata.Material(1).xLNO=1;
   plotdata.Material(1).yLNO=1;
   plotdata.Material(1).xLSCO=0;
   plotdata.Material(1).xNCM=0;
   plotdata.Material(1).xNLN=0;
   plotdata.Material(1).xNNO=1;
   plotdata.Material(1).yNNO=1;
   plotdata.Material(1).xNSM=0;
   plotdata.Material(1).xPST=0;
   plotdata.Material(1).xPZT=0;
   plotdata.Material(1).xYTmIG=0;
   plotdata.Material(1).xZMO=1;
   plotdata.Material(1).V=4^2; %unit cell volume in [A]
   % Min-Max:
   plotdata.Material(1).AexpMax=10;      
   plotdata.Material(1).AexpMin=-10;     
   plotdata.Material(1).BexpMax=1000;    
   plotdata.Material(1).BexpMin=-1000;   
   plotdata.Material(1).MeandexpMax=14;  
   plotdata.Material(1).MeandexpMin=0;   
   plotdata.Material(1).dMax=60;         
   plotdata.Material(1).dMin=0;          
   plotdata.Material(1).NMax=1000;       
   plotdata.Material(1).NMin=0;   
   plotdata.Material(1).xBSTMax=1;
   plotdata.Material(1).xBSTMin=0;
   plotdata.Material(1).xCSAMax=3;
   plotdata.Material(1).xCSAMin=0;
   plotdata.Material(1).xHZOMax=1;
   plotdata.Material(1).xHZOMin=0;
   plotdata.Material(1).xLNOMax=1;
   plotdata.Material(1).xLNOMin=0;
   plotdata.Material(1).yLNOMax=1;
   plotdata.Material(1).yLNOMin=0;
   plotdata.Material(1).xLSCOMax=1;
   plotdata.Material(1).xLSCOMin=0;
   plotdata.Material(1).xNCMMax=1;
   plotdata.Material(1).xNCMMin=0;
   plotdata.Material(1).xNLNMax=1;
   plotdata.Material(1).xNLNMin=0;
   plotdata.Material(1).xNNOMax=1;
   plotdata.Material(1).xNNOMin=0;
   plotdata.Material(1).yNNOMax=1;
   plotdata.Material(1).yNNOMin=0;
   plotdata.Material(1).xNSMMax=1;
   plotdata.Material(1).xNSMMin=0;
   plotdata.Material(1).xPSTMax=1;
   plotdata.Material(1).xPSTMin=0;
   plotdata.Material(1).xPZTMax=1;
   plotdata.Material(1).xPZTMin=0;
   plotdata.Material(1).xYTmIGMax=3;
   plotdata.Material(1).xYTmIGMin=0;
   plotdata.Material(1).xZMOMax=1;
   plotdata.Material(1).xZMOMin=0.7;
   plotdata.Material(1).PolarizationMax=2;     
   plotdata.Material(1).PolarizationMin=-2;    
   % Material 1, 2, 3 and top layer in superlattice
   for k=2:6,
       plotdata.Material(k)=plotdata.Material(1);
       plotdata.Material(k)=plotdata.Material(1);
   end;
   % dDisplay
   plotdata.dDisplay.Width=660;
   plotdata.dDisplay.Height=400;
   plotdata.dDisplay.centerfig=centerfig(plotdata.dDisplay.Width,plotdata.dDisplay.Height)+[400,-250,0,0];
   % FitDisplay
   plotdata.FitDisplay.haxis='2Theta';
   plotdata.FitDisplay.Width=660;
   plotdata.FitDisplay.Height=400;
   plotdata.FitDisplay.centerfig=centerfig(plotdata.FitDisplay.Width,plotdata.FitDisplay.Height)-[400,0,0,0];
   plotdata.LMin=0;
   plotdata.LMax=2*plotdata.Substrate.d/plotdata.lambda; 
   plotdata.ScalingIntensity=1;
   plotdata.TThetaMin=0;
   plotdata.TThetaMax=180;
   % Main Panel
   plotdata.MainPanel.Width=792;    
   plotdata.MainPanel.Height=420;   
   plotdata.MainPanel.centerfig=centerfig(plotdata.MainPanel.Width,plotdata.MainPanel.Height)+[400,+250,0,0];
   plotdata.Repetition.N=1;
   plotdata.Repetition.Max=50;
   plotdata.Repetition.Min=0;
   % Misc
   plotdata.epsilon=0.001;
   plotdata.epsilonMax=10000;       
   plotdata.epsilonMin=-10000;
   plotdata.fit.compare.x=[];
   plotdata.fit.compare.y=[];
   plotdata.fit.x.TTheta=[0:0.01:180]';
   plotdata.fit.x.L=2*plotdata.Substrate.d*sin(plotdata.fit.x.TTheta/2*pi/180)/plotdata.lambda;
   plotdata.fit.y=10^10+zeros(size([0:0.01:180]'));
   plotdata.mMax=1;
   plotdata.mMin=0;
   plotdata.measure.compare.y=[];
   plotdata.measure.x=[];
   plotdata.measure.y=[];
   plotdata.measure.filename='none';
   plotdata.prog = mfilename;
   plotdata.Q=2*pi*(2*sin([0:0.01:180]'/2*pi/180)/plotdata.lambda)'; %Momentum transfert
   plotdata.QQ=plotdata.Q.^2;
   plotdata.Substrate.g=0.*plotdata.Q;
   plotdata.d=[];
   plotdata.LayerNumber=0; %increment taking into account the layer number

   %% Create Main Panel
   plotdata.MainPanel.fig = figure('Position',plotdata.MainPanel.centerfig,...
      'Resize','on',...
      'NumberTitle','off',...
      'Name','XRD fitting parameters @Celine Lichtensteiger, J.Appl.Cryst.(2018)51(6)p.1745',...
      'Interruptible','off',...
      'Menubar','none',...
      'Color',get(0,'DefaultUIControlBackgroundColor'));
   set(plotdata.MainPanel.fig,'toolbar','figure');
   set(plotdata.MainPanel.fig,'menubar','figure');
     %==Text====================================
    %uicontrol(plotdata.MainPanel.fig,'style','text',...
    %   'string','Simulates XRD diffraction for heterostructures.',...
    %   'position',[(plotdata.MainPanel.Width-450)/2,plotdata.MainPanel.Height-30,450,25]);
   
   %==file name====================================
   plotdata.filenameWrite=uicontrol(plotdata.MainPanel.fig,'Style','text', ...
      'String','Load the data you want to fit by pressing "Load data":',...
      'Units','normalized',...
      'Position',[0.5 0.95 0.5 0.05]);
     % 'Position',[(plotdata.MainPanel.Width-400)/2,plotdata.MainPanel.Height-45,400,25]); 
   
   uicontrol(plotdata.MainPanel.fig,'Style','pushbutton',...
      'Units','normalized',...
      'Position',[0.84 0.9 0.08 0.05],...
      'String','Load data',...
      'Interruptible','off',...
      'BusyAction','cancel',...
      'Callback',[plotdata.prog,' LoadDataFilename']);
  
   %===The quit button===============================
   uicontrol(plotdata.MainPanel.fig,'Style','pushbutton',...
       'Unit','normalized',...
      'Position',[0.92 0.95 0.08 0.05],...
      'String','Quit',...
      'Interruptible','off',...
      'BusyAction','cancel',...
      'Callback',[plotdata.prog,' quit']);
  
  %===The export button===============================
   uicontrol(plotdata.MainPanel.fig,'Style','pushbutton',...
      'Units','normalized',...
      'Position',[0.92 0.9 0.08 0.05],...
      'String','Export fit',...
      'Interruptible','off',...
      'BusyAction','cancel',...
      'Callback',[plotdata.prog,' export']);
  
    %uicontrol(plotdata.MainPanel.fig,'style','text',...
    %   'string','Export the fit as a .csv file.',...
    %   'position',[(plotdata.MainPanel.Width-450)/2,plotdata.MainPanel.Height-60,450,25]);

%==Substrate chooser=========================
plotdata.Substrate.Panel=uipanel(plotdata.MainPanel.fig,'Title','Substrate:',...
    'BorderType','etchedout',...
    'position',[0.005,0.01,0.16,0.84]); %'position',[0.005,0.35,0.16,0.5]);  %'position',[0.17+(k-1)*0.135,0.01,0.13,0.84]);
plotdata.Substrate.choose=uicontrol(plotdata.Substrate.Panel,'Style','popupmenu', ...
      'String','none|Al2O3|DyScO3|beta-Ga2O3|Gd3Ga5O12|GdScO3|KTaO3|LaAlO3|LaGaO3|LSAT|MgO|Nb:SrTiO30.5%wt|NdAlO3|NdGaO3|NdScO3|PMN-PT|Si|SmScO3|SrLaAlO4|SrLaGaO4|SrTiO3|TbScO3|TiO2|YAlO3|YSZ|ZnO',...
      'value',1,...
      'Units','normalized',...
      'Position',[0,0.47,1,0.5],...
      'visible','on',...
      'CallBack',[plotdata.prog,' chooseSubstrate'] ); 
  
uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string',['Structure:'],...
    'Units','normalized','pos',[0,0.85,0.5,0.05]);
plotdata.Substrate.StructureWrite=uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string','Cubic',...
    'Units','normalized','pos',[0.5,0.85,0.5,0.05]);
uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string','a =',...
    'Units','normalized','pos',[0,0.8,0.35,0.05]);
plotdata.Substrate.aWrite=uicontrol(plotdata.Substrate.Panel,'style','edit',...
    'string',num2str(SubstrateData(condition).a),...
    'Units','normalized','pos',[0.35,0.8,0.35,0.05],...
    'callback', [plotdata.prog,' writeSubstrate_a']);
uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string',['A'],...
    'Units','normalized','pos',[0.7,0.8,0.2,0.05]); 
uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string','b =',...
    'Units','normalized','pos',[0,0.75,0.35,0.05]);
plotdata.Substrate.bWrite=uicontrol(plotdata.Substrate.Panel,'style','edit',...
    'string',num2str(SubstrateData(condition).b),...
    'Units','normalized','pos',[0.35,0.75,0.35,0.05],...
    'callback', [plotdata.prog,' writeSubstrate_b']);
uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string',['A'],...
    'Units','normalized','pos',[0.7,0.75,0.2,0.05]); 
uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string','c =',...
    'Units','normalized','pos',[0,0.7,0.35,0.05]);
plotdata.Substrate.cWrite=uicontrol(plotdata.Substrate.Panel,'style','edit',...
    'string',num2str(SubstrateData(condition).c),...
    'Units','normalized','pos',[0.35,0.7,0.35,0.05],...
    'callback', [plotdata.prog,' writeSubstrate_c']);
uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string',['A'],...
    'Units','normalized','pos',[0.7,0.7,0.2,0.05]); 
plotdata.Substrate.alphaWrite=uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string',['alpha=',num2str(SubstrateData(condition).alpha),'°'],...
    'Units','normalized','pos',[0,0.65,0.9,0.05]);
plotdata.Substrate.betaWrite=uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string',['beta=',num2str(SubstrateData(condition).beta),'°'],...
    'Units','normalized','pos',[0,0.6,0.9,0.05]);
plotdata.Substrate.gammaWrite=uicontrol(plotdata.Substrate.Panel,'style','text',...
    'string',['gamma=',num2str(SubstrateData(condition).gamma),'°'],...
    'Units','normalized','pos',[0,0.55,0.9,0.05]); 
%change N Substrate
uicontrol(plotdata.Substrate.Panel,'style','text',...
      'string','N of layers:',...
      'Units','normalized',...
      'Position',[0,0.5,0.5,0.05]);
plotdata.Substrate.NWrite = uicontrol(plotdata.Substrate.Panel,'Style','edit',...
      'String',num2str(plotdata.Substrate.N),...
      'Units','normalized',...
      'Position',[0.5,0.5,0.3,0.05],...
      'callback', [plotdata.prog,' writeNSubstrate']);
  
%%Cubic
plotdata.Substrate.TBorientation.Cubic=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','Choose orientation:',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.5]); 
% Create two buttons in the button group.
plotdata.Substrate.TBc001 = uicontrol('Style','radiobutton','String','(001)c',...
    'Units','normalized','pos',[0,0.8,1,0.15],'parent',plotdata.Substrate.TBorientation.Cubic);
%uicontrol(plotdata.Substrate.TBorientation.cubic,'style','text','string',['d=',num2str(round(plotdata.Substrate.d,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.TBc111 = uicontrol('Style','radiobutton','String','(111)c',...
    'Units','normalized','pos',[0,0.5,1,0.15],'parent',plotdata.Substrate.TBorientation.Cubic);
%uicontrol(plotdata.Substrate.TBorientation,'style','text','string',['d=',num2str(round(plotdata.Substrate.d,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.TBorientation.Cubic,'SelectionChangeFcn',@toggleSubstrateOrientation);

%%Hexagonal
plotdata.Substrate.TBorientation.Hexagonal=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','Hexagonal structure :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.5]);  
% Create three buttons in the button group.
plotdata.Substrate.TBh0001 = uicontrol('Style','radiobutton','String','(0001)h',...
    'Units','normalized','pos',[0,0.8,1,0.15],'parent',plotdata.Substrate.TBorientation.Hexagonal);
%uicontrol(plotdata.Substrate.TBorientation.hexagonal,'style','text','string',['d=',num2str(round(plotdata.Substrate.d,4)),'A'],'Units','normalized','pos',[0.2,0.7,0.8,0.15]);
plotdata.Substrate.TBh11bar20 = uicontrol('Style','radiobutton','String','(11-20)h',...
    'Units','normalized','pos',[0,0.5,1,0.15],'parent',plotdata.Substrate.TBorientation.Hexagonal);
%uicontrol(plotdata.Substrate.TBorientation,'style','text','string',['d=',num2str(round(plotdata.Substrate.d,4)),'A'],'Units','normalized','pos',[0.2,0.35,0.8,0.15]);
plotdata.Substrate.TBh10bar10 = uicontrol('Style','radiobutton','String','(10-10)h',...
    'Units','normalized','pos',[0,0.2,1,0.15],'parent',plotdata.Substrate.TBorientation.Hexagonal);
%uicontrol(plotdata.Substrate.ZnO.TBorientation,'style','text','string',['d=',num2str(round(plotdata.Substrate.d,4)),'A'],'Units','normalized','pos',[0.2,0,0.8,0.15]);
set(plotdata.Substrate.TBorientation.Hexagonal,'SelectionChangeFcn',@toggleSubstrateOrientation);

%%Monoclinic
plotdata.Substrate.TBorientation.Monoclinic=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','Choose orientation:',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.5]); 
% Create one button in the button group.
plotdata.Substrate.TBm010 = uicontrol('Style','radiobutton','String','(010)m',...
    'Units','normalized','pos',[0,0.8,1,0.15],'parent',plotdata.Substrate.TBorientation.Monoclinic);
 
%%Orthorhombic
plotdata.Substrate.TBorientation.Orthorhombic=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','Orthorhomic structure :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.5]);
uicontrol(plotdata.Substrate.TBorientation.Orthorhombic,'style','text',...
       'unit','normalized','position',[0,0.8,1,0.2]);
       %'string',['Orthorhombic a=',num2str(round(DyScO3o110.a,4)),'A b=',num2str(round(DyScO3o110.b,4)),'A c=',num2str(round(DyScO3o110.c,4)),'A'],...% Meley, APL Materials, 6, 046102 (2018) 1/d=sqrt(h^2/a^2+k^2/b^2+l^2/v^2)',...
% Create three buttons in the button group.
plotdata.Substrate.TBo001 = uicontrol('Style','radiobutton','String','(001)o',...
    'Units','normalized','pos',[0,0.8,1,0.15],'parent',plotdata.Substrate.TBorientation.Orthorhombic);
plotdata.Substrate.TBo110 = uicontrol('Style','radiobutton','String','(110)o',...
    'Units','normalized','pos',[0,0.5,1,0.15],'parent',plotdata.Substrate.TBorientation.Orthorhombic);
%uicontrol(plotdata.Substrate.TBorientation.Orthorhombic,'style','text','string',['equivalent to (001)pc, d=',num2str(round(plotdata.Substrate.d,4)),'A'],'Units','normalized','pos',[0,0.35,0.8,0.2]);
plotdata.Substrate.TBo101 = uicontrol('Style','radiobutton','String','(101)o',...
    'Units','normalized','pos',[0,0.2,1,0.15],'parent',plotdata.Substrate.TBorientation.Orthorhombic);
%uicontrol(plotdata.Substrate.TBorientation.Othorhombic,'style','text','string',['equivalent to (111)pc, d=',num2str(round(plotdata.Substrate.d,4)),'A'],'Units','normalized','pos',[0,0,0.8,0.2]);
set(plotdata.Substrate.TBorientation.Orthorhombic,'SelectionChangeFcn',@toggleSubstrateOrientation);

%%Pseudo-cubic
plotdata.Substrate.TBorientation.PseudoCubic=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','Pseudo-cubic structure :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.5]);  
% Create two buttons in the button group.
plotdata.Substrate.TBpc001 = uicontrol('Style','radiobutton','String','(001)pc',...
    'Units','normalized','pos',[0,0.8,1,0.15],'parent',plotdata.Substrate.TBorientation.PseudoCubic);
%uicontrol(plotdata.Substrate.TBorientation.PseudoCubic,'style','text','string',['d=',num2str(round(plotdata.Substrate.d,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.TBpc111 = uicontrol('Style','radiobutton','String','(111)pc',...
    'Units','normalized','pos',[0,0.5,1,0.15],'parent',plotdata.Substrate.TBorientation.PseudoCubic);
%uicontrol(plotdata.Substrate.TBorientation.PseudoCubic,'style','text','string',['d=',num2str(round(plotdata.Substrate.d,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.TBorientation.PseudoCubic,'SelectionChangeFcn',@toggleSubstrateOrientation);

%%Tetragonal
plotdata.Substrate.TBorientation.Tetragonal=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','Tetragonal structure :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.5]);  
% Create one button in the button group.
plotdata.Substrate.TBt001 = uicontrol('Style','radiobutton','String','(001)t',...
    'Units','normalized','pos',[0,0.8,1,0.15],'parent',plotdata.Substrate.TBorientation.Tetragonal);
%uicontrol(plotdata.Substrate.TBorientation.Tetragonal,'style','text','string',['d=',num2str(round(plotdata.Substrate.d,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
   
%change epsilon
   plotdata.chooseEpsilon=uipanel(plotdata.MainPanel.fig, ...
      'Title','Absorption coefficient:',...
      'BorderType','etchedout',...
      'Position',[0.01,0.9,0.15,0.1]);
   plotdata.epsilonWrite = uicontrol(plotdata.chooseEpsilon,'Style','edit',...
      'String',num2str(plotdata.epsilon),...
      'Units','normalized',...
      'Position',[0.01,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writeEpsilon']);
  
  %==Superlattice====================================
  %uipanel(plotdata.MainPanel.fig,'BorderType','etchedout','position',[0.3308,0.0001,0.4988,0.86]);
  plotdata.Superlattice.Panel=uipanel(plotdata.MainPanel.fig,'BorderType','etchedout','position',[0.3,0.0001,0.545,0.95],'Title','Superlattice:');
%  plotdata.Superlattice.Panel=uipanel(plotdata.MainPanel.fig,'BorderType','etchedout','position',[0.3308,0.0001,0.4988,0.95],'Title','Superlattice:');

  uicontrol(plotdata.Superlattice.Panel,'style','text',...
       'string','|--------------------------------------------------------------------------------------------------------|',...
       'Units','normalized',...
       'position',[0,0.9,1,0.05]);
   
   %slider for Repetition
   uicontrol(plotdata.Superlattice.Panel,'Style','text', ...
      'String','Number of repetitions:',...
      'Units','normalized',...
      'Position',[0.1,0.94,0.49,0.05]);
   plotdata.RepetitionWrite = uicontrol(plotdata.Superlattice.Panel,'Style','edit',...
      'String',num2str(plotdata.Repetition.N),...
      'Units','normalized',...
      'Position',[0.5,0.945,0.1,0.05],...
      'callback', [plotdata.prog,' writeRepetition']);
   plotdata.RepetitionSlide = uicontrol(plotdata.Superlattice.Panel,...
      'Style','slider',...
      'Min' ,plotdata.Repetition.Min,'Max',plotdata.Repetition.Max, ...
      'Units','normalized',...
      'Position',[0.6,0.94,0.3,0.05], ...
      'Value', plotdata.Repetition.N,...
      'SliderStep',[1/plotdata.Repetition.Max 1/plotdata.Repetition.Max], ...
      'CallBack', [plotdata.prog,' slideRepetition']);
   
  %==Material chooser and slider ========================= 
  for k=1:6,
      
   plotdata.Material(k).Panel=uipanel(plotdata.MainPanel.fig,... 
       'BorderType','etchedout',...
       'position',[0.17+(k-1)*0.135,0.01,0.13,0.84]);
    %'position',[0.17+(k-1)*0.165,0.01,0.16,0.84]);
   plotdata.Material(k).choose=uicontrol(plotdata.Material(k).Panel,'Style','popupmenu', ...
      'String','none|AlO2|BaBiO3|BaO|BaSnO3|(Ba_x,Sr_{1-x})TiO3|BaTiO3|BiFeO3|Bi2Te3|CaCuO2|Ca2RuO4|(Ca_{3-x},Sr_x)Al2O6|CaTiO3|CaVO3|CoO|mu-Fe2O3|beta-Ga2O3|GeTe|(Hf_{1-x},Zr_x)O2|In2O3|LaAlO3|LaCoO3|LaCrO3|La2CuO4|LaFeO3|LaMnO3|LaNiO3|La_xNi_yO3|La2NiMnO6|La2NiO4|LaO|(La_x,Sr_{1-x})CoO3|LaTiO3|LaVO3|LSMO (La0.67Sr0.33MnO3)|Mg3N2|MgO|MnO|MnO2|MnTiO3|MoS2(P.-3.m.1)|MoS2(P.63/m.m.c)|(Nd_{1-x},Ca_x)MnO3|NdCaMn2O6|(Nd_x,La_{1-x})NiO3|NdNiO2|NdNiO3|Nd_xNi_yO3|Nd2NiMnO6|NdO|(Nd_{1-x},Sr_x)MnO3|NiO2|PbO|PbNiO3|(Pb_x,Sr_{1-x})TiO3|PbTiO3|Pb(Zr_x,Ti_{1-x})O3|PrBa2Cu3O7|PrNiO2|PrNiO3|PrVO3|RuO2|SmNiO3|Sr3Al2O6|SrCoO2.5|SrCoO3|SrCrO3|SrCuO2|SrIrO3|SrMoO3|SrO|SrO2|SrRuO3|SrTiO3|SrVO3|TiO2|Tm3Fe5O12|VO2|WS2|YBa2Cu3O7|YBiO3|Y3Fe5O12|(Y_xTm_{3-x})Fe5O12|YNiO3|(Zn_x,Mg_{1-x})O|Zn3N2|ZnO',...      
      'value',1,...
      'Units','normalized',...
      'visible','on',...
      'Position',[0,0.88,1,0.1],...
      'CallBack',[plotdata.prog,' chooseLayerType(',num2str(k),')'] );
  
  %==orientation chooser =========================
  plotdata.Material(k).OrientationPanel=uipanel(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','Orientation:',...
      'BorderType','etchedout',...
      'Position',[0,0.8,1,0.1]);
  %Case hexagonal: 
  plotdata.Material(k).choosehOrientation=uicontrol('Parent',plotdata.Material(k).OrientationPanel,...
      'Style','popupmenu',...
      'String','h(0001)|h(11-20)|h(10-10)',... 
      'value',1,...
      'visible','off',...
      'Units','normalized',...
      'Position',[0,0,1,1],...
      'CallBack',[plotdata.prog,' choosehOrientation(',num2str(k),')'] );
  %Case monoclinic: 
  plotdata.Material(k).choosemOrientation=uicontrol('Parent',plotdata.Material(k).OrientationPanel,...
      'Style','popupmenu',...
      'String','m(010)',... 
      'value',1,...
      'visible','off',...
      'Units','normalized',...
      'Position',[0,0,1,1],...
      'CallBack',[plotdata.prog,' choosemOrientation(',num2str(k),')'] );
  %Case pseudo-cubic: 
  plotdata.Material(k).choosepcOrientation=uicontrol('Parent',plotdata.Material(k).OrientationPanel,...
      'Style','popupmenu',...
      'String','pc(001)|pc(111)',... 
      'value',1,...
      'visible','off',...
      'Units','normalized',...
      'Position',[0,0,1,1],...
      'CallBack',[plotdata.prog,' choosepcOrientation(',num2str(k),')'] );


  %==polarization chooser and slider =========================
  plotdata.Material(k).TBPolarization=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','Polarization:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);
  
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
  
  %==BST concentration chooser and slider =========================
  plotdata.Material(k).xBSTpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xBSTWrite = uicontrol(plotdata.Material(k).xBSTpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xBST),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexBST(',num2str(k),')']);
  
  plotdata.Material(k).xBSTSlide = uicontrol(plotdata.Material(k).xBSTpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xBSTMin,'Max',plotdata.Material(k).xBSTMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xBST,...
      'SliderStep',[0.1/(plotdata.Material(k).xBSTMax-plotdata.Material(k).xBSTMin) 0.1/(plotdata.Material(k).xBSTMax-plotdata.Material(k).xBSTMin)], ...
      'CallBack', [plotdata.prog,' slidexBST(',num2str(k),')']);

    %==CSA concentration chooser and slider =========================
  plotdata.Material(k).xCSApanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xCSAWrite = uicontrol(plotdata.Material(k).xCSApanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xCSA),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexCSA(',num2str(k),')']);
  
  plotdata.Material(k).xCSASlide = uicontrol(plotdata.Material(k).xCSApanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xCSAMin,'Max',plotdata.Material(k).xCSAMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xCSA,...
      'SliderStep',[0.1/(plotdata.Material(k).xCSAMax-plotdata.Material(k).xCSAMin) 0.1/(plotdata.Material(k).xCSAMax-plotdata.Material(k).xCSAMin)], ...
      'CallBack', [plotdata.prog,' slidexCSA(',num2str(k),')']);
  
   %==HZO concentration chooser and slider =========================
  plotdata.Material(k).xHZOpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xHZOWrite = uicontrol(plotdata.Material(k).xHZOpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xHZO),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexHZO(',num2str(k),')']);
  
  plotdata.Material(k).xHZOSlide = uicontrol(plotdata.Material(k).xHZOpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xHZOMin,'Max',plotdata.Material(k).xHZOMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xHZO,...
      'SliderStep',[0.01/(plotdata.Material(k).xHZOMax-plotdata.Material(k).xHZOMin) 0.01/(plotdata.Material(k).xHZOMax-plotdata.Material(k).xHZOMin)], ...
      'CallBack', [plotdata.prog,' slidexHZO(',num2str(k),')']);
    
  %==LSCO concentration chooser and slider =========================
  plotdata.Material(k).xLSCOpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xLSCOWrite = uicontrol(plotdata.Material(k).xLSCOpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xLSCO),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexLSCO(',num2str(k),')']);
  
  plotdata.Material(k).xLSCOSlide = uicontrol(plotdata.Material(k).xLSCOpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xLSCOMin,'Max',plotdata.Material(k).xLSCOMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xLSCO,...
      'SliderStep',[0.1/(plotdata.Material(k).xLSCOMax-plotdata.Material(k).xLSCOMin) 0.1/(plotdata.Material(k).xLSCOMax-plotdata.Material(k).xLSCOMin)], ...
      'CallBack', [plotdata.prog,' slidexLSCO(',num2str(k),')']);
  
  %==NCM concentration chooser and slider =========================
  plotdata.Material(k).xNCMpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xNCMWrite = uicontrol(plotdata.Material(k).xNCMpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xNCM),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexNCM(',num2str(k),')']);
  
  plotdata.Material(k).xNCMSlide = uicontrol(plotdata.Material(k).xNCMpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xNCMMin,'Max',plotdata.Material(k).xNCMMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xNCM,...
      'SliderStep',[0.01/(plotdata.Material(k).xNCMMax-plotdata.Material(k).xNCMMin) 0.01/(plotdata.Material(k).xNCMMax-plotdata.Material(k).xNCMMin)], ...
      'CallBack', [plotdata.prog,' slidexNCM(',num2str(k),')']);
  
  %==NLN concentration chooser and slider =========================
  plotdata.Material(k).xNLNpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xNLNWrite = uicontrol(plotdata.Material(k).xNLNpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xNLN),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexNLN(',num2str(k),')']);
  
  plotdata.Material(k).xNLNSlide = uicontrol(plotdata.Material(k).xNLNpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xNLNMin,'Max',plotdata.Material(k).xNLNMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xNLN,...
      'SliderStep',[0.01/(plotdata.Material(k).xNLNMax-plotdata.Material(k).xNLNMin) 0.01/(plotdata.Material(k).xNLNMax-plotdata.Material(k).xNLNMin)], ...
      'CallBack', [plotdata.prog,' slidexNLN(',num2str(k),')']);
  
    %==NSM concentration chooser and slider =========================
  plotdata.Material(k).xNSMpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xNSMWrite = uicontrol(plotdata.Material(k).xNSMpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xNSM),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexNSM(',num2str(k),')']);
  
  plotdata.Material(k).xNSMSlide = uicontrol(plotdata.Material(k).xNSMpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xNSMMin,'Max',plotdata.Material(k).xNSMMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xNSM,...
      'SliderStep',[0.01/(plotdata.Material(k).xNSMMax-plotdata.Material(k).xNSMMin) 0.01/(plotdata.Material(k).xNSMMax-plotdata.Material(k).xNSMMin)], ...
      'CallBack', [plotdata.prog,' slidexNSM(',num2str(k),')']);
  
  %==LNO off-stoichiometry chooser and slider =========================
  plotdata.Material(k).xLNOpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,0.5,0.1]);

  plotdata.Material(k).xLNOWrite = uicontrol(plotdata.Material(k).xLNOpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xLNO),...
      'Units','normalized',...
      'Position',[0,0.1,0.5,0.8],...
      'callback', [plotdata.prog,' writexLNO(',num2str(k),')']);
  
  plotdata.Material(k).xLNOSlide = uicontrol(plotdata.Material(k).xLNOpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xLNOMin,'Max',plotdata.Material(k).xLNOMax, ...
      'Units','normalized',...
      'Position',[0.5,0.3,0.5,0.5], ...
      'Value', plotdata.Material(k).xLNO,...
      'SliderStep',[0.01/(plotdata.Material(k).xLNOMax-plotdata.Material(k).xLNOMin) 0.01/(plotdata.Material(k).xLNOMax-plotdata.Material(k).xLNOMin)], ...
      'CallBack', [plotdata.prog,' slidexLNO(',num2str(k),')']);

   plotdata.Material(k).yLNOpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','y:',...
      'BorderType','etchedout',...
      'Position',[0.5,0.7,0.5,0.1]);
  
  plotdata.Material(k).yLNOWrite = uicontrol(plotdata.Material(k).yLNOpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).yLNO),...
      'Units','normalized',...
      'Position',[0,0.1,0.5,0.8],...
      'callback', [plotdata.prog,' writeyLNO(',num2str(k),')']);
  
  plotdata.Material(k).yLNOSlide = uicontrol(plotdata.Material(k).yLNOpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).yLNOMin,'Max',plotdata.Material(k).yLNOMax, ...
      'Units','normalized',...
      'Position',[0.5,0.3,0.5,0.5], ...
      'Value', plotdata.Material(k).yLNO,...
      'SliderStep',[0.01/(plotdata.Material(k).yLNOMax-plotdata.Material(k).yLNOMin) 0.01/(plotdata.Material(k).yLNOMax-plotdata.Material(k).yLNOMin)], ...
      'CallBack', [plotdata.prog,' slideyLNO(',num2str(k),')']);
  
    %==NNO off-stoichiometry chooser and slider =========================
  plotdata.Material(k).xNNOpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,0.5,0.1]);

  plotdata.Material(k).xNNOWrite = uicontrol(plotdata.Material(k).xNNOpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xNNO),...
      'Units','normalized',...
      'Position',[0,0.1,0.5,0.8],...
      'callback', [plotdata.prog,' writexNNO(',num2str(k),')']);
  
  plotdata.Material(k).xNNOSlide = uicontrol(plotdata.Material(k).xNNOpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xNNOMin,'Max',plotdata.Material(k).xNNOMax, ...
      'Units','normalized',...
      'Position',[0.5,0.3,0.5,0.5], ...
      'Value', plotdata.Material(k).xNNO,...
      'SliderStep',[0.01/(plotdata.Material(k).xNNOMax-plotdata.Material(k).xNNOMin) 0.01/(plotdata.Material(k).xNNOMax-plotdata.Material(k).xNNOMin)], ...
      'CallBack', [plotdata.prog,' slidexNNO(',num2str(k),')']);

   plotdata.Material(k).yNNOpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','y:',...
      'BorderType','etchedout',...
      'Position',[0.5,0.7,0.5,0.1]);
  
  plotdata.Material(k).yNNOWrite = uicontrol(plotdata.Material(k).yNNOpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).yNNO),...
      'Units','normalized',...
      'Position',[0,0.1,0.5,0.8],...
      'callback', [plotdata.prog,' writeyNNO(',num2str(k),')']);
  
  plotdata.Material(k).yNNOSlide = uicontrol(plotdata.Material(k).yNNOpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).yNNOMin,'Max',plotdata.Material(k).yNNOMax, ...
      'Units','normalized',...
      'Position',[0.5,0.3,0.5,0.5], ...
      'Value', plotdata.Material(k).yNNO,...
      'SliderStep',[0.01/(plotdata.Material(k).yNNOMax-plotdata.Material(k).yNNOMin) 0.01/(plotdata.Material(k).yNNOMax-plotdata.Material(k).yNNOMin)], ...
      'CallBack', [plotdata.prog,' slideyNNO(',num2str(k),')']);
  
  %==PST concentration chooser and slider =========================
  plotdata.Material(k).xPSTpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

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
      'Position',[0,0.7,1,0.1]);

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
 
  %==YTmIG concentration chooser and slider =========================
  plotdata.Material(k).xYTmIGpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xYTmIGWrite = uicontrol(plotdata.Material(k).xYTmIGpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xYTmIG),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexYTmIG(',num2str(k),')']);
  
  plotdata.Material(k).xYTmIGSlide = uicontrol(plotdata.Material(k).xYTmIGpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xYTmIGMin,'Max',plotdata.Material(k).xYTmIGMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xYTmIG,...
      'SliderStep',[0.1/(plotdata.Material(k).xYTmIGMax-plotdata.Material(k).xYTmIGMin) 0.1/(plotdata.Material(k).xYTmIGMax-plotdata.Material(k).xYTmIGMin)], ...
      'CallBack', [plotdata.prog,' slidexYTmIG(',num2str(k),')']);
  
    %==ZnMgO concentration chooser and slider =========================
  plotdata.Material(k).xZMOpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xZMOWrite = uicontrol(plotdata.Material(k).xZMOpanel,'Style','edit',...
      'String',num2str(plotdata.Material(k).xZMO),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexZMO(',num2str(k),')']);
  
  plotdata.Material(k).xZMOSlide = uicontrol(plotdata.Material(k).xZMOpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xZMOMin,'Max',plotdata.Material(k).xZMOMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', plotdata.Material(k).xZMO,...
      'SliderStep',[0.01/(plotdata.Material(k).xZMOMax-plotdata.Material(k).xZMOMin) 0.01/(plotdata.Material(k).xZMOMax-plotdata.Material(k).xZMOMin)], ...
      'CallBack', [plotdata.prog,' slidexZMO(',num2str(k),')']);

   %slider for N
   plotdata.Material(k).chooseN=uipanel(plotdata.Material(k).Panel, ...
      'Title','N:',...
      'BorderType','etchedout',...
      'Position',[0,0.6,1,0.1]);
   plotdata.Material(k).NWrite = uicontrol(plotdata.Material(k).chooseN,'Style','edit',...
      'String',num2str(plotdata.Material(k).N),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writeN(',num2str(k),')']);
   plotdata.Material(k).NSlide = uicontrol(plotdata.Material(k).chooseN,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).NMin,'Max',plotdata.Material(k).NMax, ...
      'Units','normalized',...
      'Position',[0.35,0.35,0.65,0.5], ...
      'Value', plotdata.Material(k).N,...
      'SliderStep',[1/(plotdata.Material(k).NMax-plotdata.Material(k).NMin) 1/(plotdata.Material(k).NMax-plotdata.Material(k).NMin)], ...
      'CallBack', [plotdata.prog,' slideN(',num2str(k),')']);   

%== d: choose between constant, exponential and sinusoidal distribution =========================

plotdata.Material(k).chooseddistrib = uipanel(plotdata.Material(k).Panel,'visible','on',...
    'Title','Choose d distribution:',...
    'BorderType','etchedout',...
    'Position',[0 0 1 0.55]);

plotdata.Material(k).ddistrib = uicontrol(plotdata.Material(k).chooseddistrib,'Style','popupmenu',...
    'String','d=constant|<html>c<sub>z</sub>=Aexp<sup>z/B</sup>+C</html>',...%    'String','d=constant|<html>d<sub>z</sub>=Aexp<sup>z/B</sup>+C</html>|Jacobi elliptic function|upload|',...
    'value',1,...
    'Units','normalized',...
    'visible','on',...
    'pos',[0,0.475,1,0.5],...
    'CallBack',[plotdata.prog,' chooseddistrib(',num2str(k),')'] );
%% d=constant
plotdata.Material(k).dconst=uipanel(plotdata.Material(k).chooseddistrib,...
    'visible','on','Position',[0 0 1 0.825]);
  
uicontrol(plotdata.Material(k).dconst,...
    'Style','text','String','  d:','Units','normalized','Position',[0 0.9 1 0.1]);

plotdata.Material(k).dWrite = uicontrol(plotdata.Material(k).dconst,'Style','edit',...
      'String',num2str(plotdata.Material(k).d),...
      'Units','normalized',...
      'Position',[0.6,0.8,0.4,0.2],...
      'callback', [plotdata.prog,' writed(',num2str(k),')']);
  
plotdata.Material(k).dSlide = uicontrol(plotdata.Material(k).dconst,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).dMin,'Max',plotdata.Material(k).dMax, ...
      'Units','normalized',...
      'Position',[0,0.6,1,0.2], ...
      'Value', plotdata.Material(k).d,...
      'SliderStep',[0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin) 0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin)], ...
      'CallBack', [plotdata.prog,' slided(',num2str(k),')']);
%% d=exp
plotdata.Material(k).dexpdistrib=uipanel(plotdata.Material(k).chooseddistrib,...
      'visible','off',...
      'BorderType','etchedout',...
      'Position',[0,0,1,0.825]);
%A  
plotdata.Material(k).chooseAexp=uipanel(plotdata.Material(k).dexpdistrib,...
    'Title','A:',...
    'BorderType','etchedout',...
    'Position',[0 0.7 1 0.3]);

plotdata.Material(k).AexpWrite = uicontrol(plotdata.Material(k).chooseAexp,'Style','edit',...
      'String',num2str(plotdata.Material(k).Aexp),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writeAexp(',num2str(k),')']);
  
plotdata.Material(k).AexpSlide = uicontrol(plotdata.Material(k).chooseAexp,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).AexpMin,'Max',plotdata.Material(k).AexpMax, ...
      'Units','normalized',...
      'Position',[0.35,0.35,0.65,0.5], ...
      'Value', plotdata.Material(k).Aexp,...
      'SliderStep',[0.001*1/(plotdata.Material(k).AexpMax-plotdata.Material(k).AexpMin) 0.001*1/(plotdata.Material(k).AexpMax-plotdata.Material(k).AexpMin)], ...
      'CallBack', [plotdata.prog,' slideAexp(',num2str(k),')']);
%B
plotdata.Material(k).chooseBexp=uipanel(plotdata.Material(k).dexpdistrib,...
    'Title','B:',...
    'BorderType','etchedout',...
    'Position',[0 0.4 1 0.3]);

plotdata.Material(k).BexpWrite = uicontrol(plotdata.Material(k).chooseBexp,'Style','edit',...
      'String',num2str(plotdata.Material(k).Bexp),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writeBexp(',num2str(k),')']);
  
plotdata.Material(k).BexpSlide = uicontrol(plotdata.Material(k).chooseBexp,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).BexpMin,'Max',plotdata.Material(k).BexpMax,...
      'Units','normalized',...
      'Position',[0.35,0.35,0.65,0.5], ...
      'Value', plotdata.Material(k).Bexp,...
      'SliderStep',[1/(plotdata.Material(k).BexpMax-plotdata.Material(k).BexpMin) 1/(plotdata.Material(k).BexpMax-plotdata.Material(k).BexpMin)], ...
      'CallBack', [plotdata.prog,' slideBexp(',num2str(k),')']);
% average d
plotdata.Material(k).chooseAveraged=uipanel(plotdata.Material(k).dexpdistrib,...
    'Title','Average d:',...
    'BorderType','etchedout',...
    'Position',[0 0.1 1 0.3]);

plotdata.Material(k).MeandexpWrite = uicontrol(plotdata.Material(k).chooseAveraged,'Style','edit',...
      'String',num2str(plotdata.Material(k).Meandexp),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writeMeandexp(',num2str(k),')']);
  
plotdata.Material(k).MeandexpSlide = uicontrol(plotdata.Material(k).chooseAveraged,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).MeandexpMin,'Max',plotdata.Material(k).MeandexpMax, ...
      'Units','normalized',...
      'Position',[0.35,0.35,0.65,0.5], ...
      'Value', plotdata.Material(k).Meandexp,...
      'SliderStep',[0.001*1/(plotdata.Material(k).MeandexpMax-plotdata.Material(k).MeandexpMin) 0.001*1/(plotdata.Material(k).MeandexpMax-plotdata.Material(k).MeandexpMin)], ...
      'CallBack', [plotdata.prog,' slideMeandexp(',num2str(k),')']);
  
uicontrol(plotdata.Material(k).dexpdistrib,'Style','text', ...
      'String','   C:',...
      'Position',[0,0,22,12]);
  
plotdata.Material(k).CexpWrite=uicontrol(plotdata.Material(k).dexpdistrib,'Style','text', ...
      'String',num2str(plotdata.Material(k).Cexp),...
      'Position',[25,0,50,12]);
%% Jacobi elliptic function
plotdata.Material(k).dJacobidistrib=uipanel(plotdata.Material(k).chooseddistrib,...
    'visible','off','Position',[0 0 1 0.825]);
  
uicontrol(plotdata.Material(k).dJacobidistrib,...
    'Style','text','String','d @center:','Units','normalized','Position',[0 0.8 0.7 0.15]);%[0 0.8 1 0.15]);

plotdata.Material(k).dCenterWrite = uicontrol(plotdata.Material(k).dJacobidistrib,'Style','edit',...
      'String',num2str(plotdata.Material(k).dCenter),...
      'Units','normalized',...
      'Position',[0.6,0.8,0.4,0.2],...
      'callback', [plotdata.prog,' writedCenter(',num2str(k),')']);
  
plotdata.Material(k).dCenterSlide = uicontrol(plotdata.Material(k).dJacobidistrib,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).dMin,'Max',plotdata.Material(k).dMax, ...
      'Units','normalized',...
      'Position',[0,0.6,1,0.2], ...
      'Value', plotdata.Material(k).dCenter,...
      'SliderStep',[0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin) 0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin)], ...
      'CallBack', [plotdata.prog,' slidedCenter(',num2str(k),')']);
  
uicontrol(plotdata.Material(k).dJacobidistrib,...
    'Style','text','String','d @border:','Units','normalized','Position',[0 0.45 0.7 0.15]);%[0 0.8 1 0.15]);

plotdata.Material(k).dBorderWrite = uicontrol(plotdata.Material(k).dJacobidistrib,'Style','edit',...
      'String',num2str(plotdata.Material(k).dBorder),...
      'Units','normalized',...
      'Position',[0.6,0.45,0.4,0.2],...
      'callback', [plotdata.prog,' writedBorder(',num2str(k),')']);
  
plotdata.Material(k).dBorderSlide = uicontrol(plotdata.Material(k).dJacobidistrib,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).dMin,'Max',plotdata.Material(k).dMax, ...
      'Units','normalized',...
      'Position',[0,0.25,1,0.2], ...
      'Value', plotdata.Material(k).dBorder,...
      'SliderStep',[0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin) 0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin)], ...
      'CallBack', [plotdata.prog,' slidedBorder(',num2str(k),')']);
  
uicontrol(plotdata.Material(k).dJacobidistrib,...
    'Style','text','String','modulus m:','Units','normalized','Position',[0 0.1 0.7 0.15]);

plotdata.Material(k).mWrite = uicontrol(plotdata.Material(k).dJacobidistrib,'Style','edit',...
      'String',num2str(plotdata.Material(k).m),...
      'Units','normalized',...
      'Position',[0.6,0.1,0.4,0.2],...
      'callback', [plotdata.prog,' writem(',num2str(k),')']);
  
plotdata.Material(k).mSlide = uicontrol(plotdata.Material(k).dJacobidistrib,...
      'Style','slider',...
      'Min' ,plotdata.mMin,'Max',plotdata.mMax, ...
      'Units','normalized',...
      'Position',[0,-0.1,1,0.2], ...
      'Value', plotdata.Material(k).m,...
      'SliderStep',[0.001*1/(plotdata.mMax-plotdata.mMin) 0.001*1/(plotdata.mMax-plotdata.mMin)], ...
      'CallBack', [plotdata.prog,' slidem(',num2str(k),')']);
  
  %% upload d-spacing in Angströms
plotdata.Material(k).dupload=uipanel(plotdata.Material(k).chooseddistrib,...
    'visible','off','Position',[0 0 1 0.825]);
  
uicontrol(plotdata.Material(k).dupload,...
    'Style','text','String','Load file containing d spacing in Angströms, number of rows = N, first row = bottom layer:','Units','normalized',...
    'Position',[0 0.5 1 0.5]);%[0 0.8 1 0.15]);

   uicontrol(plotdata.Material(k).dupload,...
      'Style','pushbutton',...
      'Units','normalized',...
      'Position',[0.25 0.25 0.5 0.2],...
      'String','Load d ',...
      'Interruptible','off',...
      'BusyAction','cancel',...
      'Callback',[plotdata.prog,' Loadd(',num2str(k),')']);

  end;

  set(plotdata.Material(1).Panel,'Title','Bottom layer:','Tag','k=1');
  set(plotdata.Material(2).Panel,'Title','Material 1:','Tag','k=2');
  set(plotdata.Material(3).Panel,'Title','Material 2:','Tag','k=3');
  set(plotdata.Material(4).Panel,'Title','Material 3:','Tag','k=4');
  set(plotdata.Material(5).Panel,'Title','Material 4:','Tag','k=5');
  set(plotdata.Material(6).Panel,'Title','Top layer:','Tag','k=6');
  
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
   subplot('Position',[0.07,0.2,0.9,0.78]);
   %subplot('Position',[0.07,0.1,0.9,0.8]);
   
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
%h = uibuttongroup('parent',plotdata.FitDisplay.fig,'visible','on','Position',[0.005,0.94,0.15/plotdata.FitDisplay.Width*660,0.06/plotdata.FitDisplay.Height*400]);
h = uibuttongroup('parent',plotdata.FitDisplay.fig,'visible','on','Units','pixel','Position',[10,30,100,25]);
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
     
   I=plotdata.Substrate.g.*conj(plotdata.Substrate.g);
   plotdata.fit.y=I/max(I);
   semilogy(plotdata.fit.x.TTheta,plotdata.fit.y*plotdata.ScalingIntensity,'b');    
   axis([plotdata.TThetaMin plotdata.TThetaMax 10^-9 1]);
   xlabel('2Theta ({\circ})'); ylabel('Intensity (arbitrary units)');
  %% Create d display
   plotdata.dDisplay.fig = figure('Position',plotdata.dDisplay.centerfig,...
      'Resize','on',...
      'NumberTitle','off',...
      'Name','d display',...
      'Interruptible','off',...
      'Menubar','none',...
      'Color',get(0,'DefaultUIControlBackgroundColor'));
   set(plotdata.dDisplay.fig,'toolbar','figure');
   set(plotdata.dDisplay.fig,'menubar','figure');
   
%%
case 'quit'
   quit_reply = questdlg('Really quit this nice fitting progam?');
   if strcmp(quit_reply,'Yes')
      close all;
      clear all;
   end
   %close all;
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
    export_reply = questdlg(strcat('Your fit (fit.csv) and paremeters (parameters.csv and parameters.mat) will be saved in ',startfolder,'. WARNING: If this folder already contains files with the same name, they will be overwritten.'));     
    csvwrite([startfolder,'fit.csv'],[x y']);
    save([startfolder,'parameters'],'plotdata');
    %save important parameters
    parameters=[];
    %Ref to measurement file
    if ~strcmp(plotdata.measure.filename,'none')
      parameters=[parameters;{'Measurement file name: '},plotdata.measure.filename];
    end
    %Substrate
    parameters=[parameters;{'Substrate type: '},{char(plotdata.Substrate.Type)};{'Substrate d[A]: '},{plotdata.Substrate.d};{'Substrate N: '},{plotdata.Substrate.N}];
    %Bottom layer
    if ~strcmp(plotdata.Material(1).Type,'none')
        parameters=[parameters;{'Bottom layer type: '},{char(plotdata.Material(1).Type)}];
        parameters=[parameters;{'Bottom layer orientation: '},{plotdata.Material(1).Orientation}];
        if strcmp(plotdata.Material(1).Type,'(Ba_x,Sr_{1-x})TiO3'),
            parameters=[parameters;{'Bottom layer Ba atomic concentration x: '},{plotdata.Material(1).xBST}];
        elseif strcmp(plotdata.Material(1).Type,'(Ca_{3-x},Sr_x)Al2O6'),
            parameters=[parameters;{'Bottom layer Sr atomic concentration x: '},{plotdata.Material(1).xCSA}];
        elseif strcmp(plotdata.Material(1).Type,'(Hf_{1-x},Zr_x)O2'),
            parameters=[parameters;{'Bottom layer Zr atomic concentration x: '},{plotdata.Material(1).xHZO}];
        elseif strcmp(plotdata.Material(1).Type,'La_xNi_yO3'),
            parameters=[parameters;{'Bottom layer La off-stochiometry x: '},{plotdata.Material(1).xLNO}];
            parameters=[parameters;{'Bottom layer Ni off-stochiometry y: '},{plotdata.Material(1).yLNO}];
        elseif strcmp(plotdata.Material(1).Type,'(La_x,Sr_{1-x})CoO3'),
            parameters=[parameters;{'Bottom layer La atomic concentration x: '},{plotdata.Material(1).xLSCO}];
        elseif strcmp(plotdata.Material(1).Type,'(Nd_{1-x},Ca_x)MnO3'),
            parameters=[parameters;{'Bottom layer Ca atomic concentration x: '},{plotdata.Material(1).xNCM}];
        elseif strcmp(plotdata.Material(1).Type,'(Nd_x,La_{1-x})NiO3'),
            parameters=[parameters;{'Bottom layer Nd atomic concentration x: '},{plotdata.Material(1).xNLN}];
        elseif strcmp(plotdata.Material(1).Type,'Nd_xNi_yO3'),
            parameters=[parameters;{'Bottom layer Nd off-stochiometry x: '},{plotdata.Material(1).xNNO}];
            parameters=[parameters;{'Bottom layer Ni off-stochiometry y: '},{plotdata.Material(1).yNNO}];
        elseif strcmp(plotdata.Material(1).Type,'(Nd_{1-x},Sr_x)MnO3'),
            parameters=[parameters;{'Bottom layer Sr atomic concentration x: '},{plotdata.Material(1).xNSM}];
        elseif strcmp(plotdata.Material(1).Type,'(Pb_x,Sr_{1-x})TiO3'),
            parameters=[parameters;{'Bottom layer Pb atomic concentration x: '},{plotdata.Material(1).xPST}];
        elseif strcmp(plotdata.Material(1).Type,'Pb(Zr_x,Ti_{1-x})O3'),
            parameters=[parameters;{'Bottom layer Zr atomic concentration x: '},{plotdata.Material(1).xPZT}];
        elseif strcmp(plotdata.Material(1).Type,'(Y_xTm_{3-x})Fe5O12'),
            parameters=[parameters;{'Bottom layer Y atomic concentration x: '},{plotdata.Material(1).xYTmIG}];
        elseif strcmp(plotdata.Material(1).Type,'(Zn_x,Mg_{1-x})O'),
            parameters=[parameters;{'Bottom layer Zn atomic concentration x: '},{plotdata.Material(1).xZMO}];
        end;
        if strcmp(plotdata.Material(1).dDistribution,'constant')
                parameters=[parameters;{'Bottom layer d[A]: '},{round(plotdata.Material(1).d,4)}];
        elseif strcmp(plotdata.Material(1).dDistribution,'exp')
            parameters=[parameters;{'Bottom layer d(z)[A]: '},{[num2str(plotdata.Material(1).Aexp),' x exp(z/',num2str(plotdata.Material(1).Bexp),')+',num2str(plotdata.Material(1).Cexp)]};...
                {'Bottom layer <d>(average)[A]: '},plotdata.Material(1).Meandexp];
        else warning('Error in plotdata.Material(1).dDistribution in Export case');
        end;
        parameters=[parameters;{'Bottom layer N: '},{plotdata.Material(1).N}];
    end
    %Number of repetitions (superlattice) 
    if plotdata.Repetition.N~=1,
        parameters=[parameters;{'Number of superlattice repetitions [Material 1+2+3]: '},{plotdata.Repetition.N}];
    end;
    %Material 1+2+3+4
    for k=2:5,
        if ~strcmp(plotdata.Material(k).Type,'none')
            parameters=[parameters;{['Material ',num2str(k-1),' type : ']},{char(plotdata.Material(k).Type)}];
            parameters=[parameters;{['Material ',num2str(k-1),' orientation : ']},{char(plotdata.Material(k).Orientation)}];
            if strcmp(plotdata.Material(k).Type,'(Ba_x,Sr_{1-x})TiO3'),
                parameters=[parameters;{['Material ',num2str(k-1),' Ba atomic concentration x : ']},{plotdata.Material(k).xBST}];
            elseif strcmp(plotdata.Material(k).Type,'(Ca_{3-x},Sr_x)Al2O6'),
                parameters=[parameters;{['Material ',num2str(k-1),' Sr atomic concentration x : ']},{plotdata.Material(k).xCSA}];
            elseif strcmp(plotdata.Material(k).Type,'(Hf_{1-x},Zr_x)O2'),
                parameters=[parameters;{['Material ',num2str(k-1),' Zr atomic concentration x: ']},{plotdata.Material(k).xHZO}];
            elseif strcmp(plotdata.Material(k).Type,'La_xNi_yO3'),
                parameters=[parameters;{'Material ',num2str(k-1),' La off-stochiometry x: '},{plotdata.Material(k).xLNO}];
                parameters=[parameters;{'Material ',num2str(k-1),' Ni off-stochiometry y: '},{plotdata.Material(k).yLNO}];
            elseif strcmp(plotdata.Material(k).Type,'(La_x,Sr_{1-x})CoO3'),
                parameters=[parameters;{['Material ',num2str(k-1),' La atomic concentration x : ']},{plotdata.Material(k).xLSCO}];
            elseif strcmp(plotdata.Material(k).Type,'(Nd_{1-x},Ca_x)MnO3'),
                parameters=[parameters;{['Material ',num2str(k-1),' Ca atomic concentration x: ']},{plotdata.Material(k).xNCM}];
            elseif strcmp(plotdata.Material(k).Type,'(Nd_x,La_{1-x})NiO3'),
                parameters=[parameters;{['Material ',num2str(k-1),' Nd atomic concentration x: ']},{plotdata.Material(k).xNLN}];
            elseif strcmp(plotdata.Material(k).Type,'Nd_xNi_yO3'),
                parameters=[parameters;{'Material ',num2str(k-1),' Nd off-stochiometry x: '},{plotdata.Material(k).xNNO}];
                parameters=[parameters;{'Material ',num2str(k-1),' Ni off-stochiometry y: '},{plotdata.Material(k).yNNO}];
            elseif strcmp(plotdata.Material(k).Type,'(Nd_{1-x},Sr_x)MnO3'),
                parameters=[parameters;{['Material ',num2str(k-1),' Sr atomic concentration x: ']},{plotdata.Material(k).xNSM}];
            elseif strcmp(plotdata.Material(k).Type,'(Pb_x,Sr_{1-x})TiO3'),
                parameters=[parameters;{['Material ',num2str(k-1),' Pb atomic concentration x : ']},{plotdata.Material(k).xPST}];
            elseif strcmp(plotdata.Material(k).Type,'Pb(Zr_x,Ti_{1-x})O3'),
                parameters=[parameters;{['Material ',num2str(k-1),' Zr atomic concentration x : ']},{plotdata.Material(k).xPZT}];
            elseif strcmp(plotdata.Material(k).Type,'(Y_xTm_{3-x})Fe5O12'),
                parameters=[parameters;{['Material ',num2str(k-1),' Y atomic concentration x: ']},{plotdata.Material(k).xYTmIG}];
            elseif strcmp(plotdata.Material(k).Type,'(Zn_x,Mg_{1-x})O'),
                parameters=[parameters;{['Material ',num2str(k-1),' Zn atomic concentration x: ']},{plotdata.Material(k).xZMO}];
            end;
            if strcmp(plotdata.Material(k).dDistribution,'constant')
                    parameters=[parameters;{['Material ',num2str(k-1),' d[A]: ']},{round(plotdata.Material(k).d,4)}];
        elseif strcmp(plotdata.Material(k).dDistribution,'exp')
            parameters=[parameters;{['Material ',num2str(k-1),' d(z)[A]: ']},{[num2str(plotdata.Material(k).Aexp),' x exp(z/',num2str(plotdata.Material(k).Bexp),')+',num2str(plotdata.Material(k).Cexp)]};...
                {['Material ',num2str(k-1),' <d>(average)[A]: ']},plotdata.Material(k).Meandexp];
        else warning('Error in plotdata.Material(k).dDistribution in Export case');
        end;
        parameters=[parameters;{['Material ',num2str(k-1),' N: ']},{plotdata.Material(k).N}];
        end
    end;
    %Top layer
    if ~strcmp(plotdata.Material(6).Type,'none')
        parameters=[parameters;{'Top layer type: '},{char(plotdata.Material(6).Type)}];
        parameters=[parameters;{'Top layer orientation: '},{char(plotdata.Material(6).Orientation)}];
        if strcmp(plotdata.Material(6).Type,'(Ba_x,Sr_{1-x})TiO3'),
            parameters=[parameters;{'Top layer Ba atomic concentration x: '},{plotdata.Material(6).xBST}];
        elseif strcmp(plotdata.Material(6).Type,'(Ca_{3-x},Sr_x)Al2O6'),
            parameters=[parameters;{'Top layer Sr atomic concentration x: '},{plotdata.Material(6).xCSA}];
        elseif strcmp(plotdata.Material(6).Type,'(Hf_{1-x},Zr_x)O2'),
            parameters=[parameters;{'Top layer Zr atomic concentration x: '},{plotdata.Material(6).xHZO}]; 
        elseif strcmp(plotdata.Material(6).Type,'La_xNi_yO3'),
            parameters=[parameters;{'Top layer La off-stochiometry x: '},{plotdata.Material(6).xLNO}];
            parameters=[parameters;{'Top layer Ni off-stochiometry y: '},{plotdata.Material(6).yLNO}];
                    elseif strcmp(plotdata.Material(6).Type,'(La_x,Sr_{1-x})CoO3'),
            parameters=[parameters;{'Top layer La atomic concentration x: '},{plotdata.Material(6).xLSCO}];
        elseif strcmp(plotdata.Material(6).Type,'(Nd_{1-x},Ca_x)MnO3'),
            parameters=[parameters;{'Top layer Ca atomic concentration x: '},{plotdata.Material(6).xNCM}];  
        elseif strcmp(plotdata.Material(6).Type,'(Nd_x,La_{1-x})NiO3'),
            parameters=[parameters;{'Top layer Nd atomic concentration x: '},{plotdata.Material(6).xNLN}];     
        elseif strcmp(plotdata.Material(6).Type,'Nd_xNi_yO3'),
            parameters=[parameters;{'Top layer Nd off-stochiometry x: '},{plotdata.Material(6).xNNO}];
            parameters=[parameters;{'Top layer Ni off-stochiometry y: '},{plotdata.Material(6).yNNO}];
        elseif strcmp(plotdata.Material(6).Type,'(Nd_{1-x},Sr_x)MnO3'),
            parameters=[parameters;{'Top layer Sr atomic concentration x: '},{plotdata.Material(6).xNSM}];
        elseif strcmp(plotdata.Material(6).Type,'(Pb_x,Sr_{1-x})TiO3'),
            parameters=[parameters;{'Top layer Pb atomic concentration x: '},{plotdata.Material(6).xPST}];
        elseif strcmp(plotdata.Material(6).Type,'Pb(Zr_x,Ti_{1-x})O3'),
            parameters=[parameters;{'Top layer Zr atomic concentration x: '},{plotdata.Material(6).xPZT}];
        elseif strcmp(plotdata.Material(6).Type,'(Y_xTm_{3-x})Fe5O12'),
            parameters=[parameters;{'Top layer Y atomic concentration x: '},{plotdata.Material(6).xYTmIG}]; 
        elseif strcmp(plotdata.Material(6).Type,'(Zn_x,Mg_{1-x})O'),
            parameters=[parameters;{'Top layer Zn atomic concentration x: '},{plotdata.Material(6).xZMO}];
        end;
        if strcmp(plotdata.Material(6).dDistribution,'constant')
                parameters=[parameters;{'Top layer d[A]: '},{round(plotdata.Material(6).d,4)}];
        elseif strcmp(plotdata.Material(6).dDistribution,'exp')
            parameters=[parameters;{'Top layer d(z)[A]: '},{[num2str(plotdata.Material(6).Aexp),' x exp(z/',num2str(plotdata.Material(6).Bexp),')+',num2str(plotdata.Material(6).Cexp)]};...
                {'Top layer <d>(average)[A]: '},plotdata.Material(6).Meandexp];
        else warning('Error in plotdata.Material(6).dDistribution in Export case');
        end;
        parameters=[parameters;{'Top layer N: '},{plotdata.Material(6).N}];
    end
    parameters=cell2table(parameters,'VariableNames',{'Parameters','Values'});
    writetable(parameters,[startfolder,'parameters.csv']);
   end
   %%
    case {'writeSubstrate_a'}
        a=str2num(get(plotdata.Substrate.aWrite,'string'));
        if strcmp(plotdata.Substrate.Type,'none')
            set(plotdata.Substrate.aWrite,'String','');
            set(plotdata.Substrate.bWrite,'String','');
            set(plotdata.Substrate.cWrite,'String','');
            set(plotdata.Substrate.alphaWrite,'String',['alpha=','','°']);
            set(plotdata.Substrate.betaWrite,'String',['beta=','','°']);
            set(plotdata.Substrate.gammaWrite,'String',['gamma=','','°']);
        else
        load('Substrates.mat')
        SubstrateTable=struct2table(SubstrateData);
        condition=find(strcmp(SubstrateTable.Name,plotdata.Substrate.Type));
        for k=1:size(condition),
            SubstrateData(condition(k)).a=a;
        end;
        if strcmp(SubstrateData(condition(1)).Structure,'Cubic'),
            for k=1:size(condition),
                SubstrateData(condition(k)).b=a;
                SubstrateData(condition(k)).c=a;
            end
            set(plotdata.Substrate.bWrite,'String',num2str(a));
            set(plotdata.Substrate.cWrite,'String',num2str(a));
        end,
        save('Substrates.mat','SubstrateData','-append');
        Substrate
        ProgcNSimu
        end,
        
    case {'writeSubstrate_b'}
        b=str2num(get(plotdata.Substrate.bWrite,'string'));        
        if strcmp(plotdata.Substrate.Type,'none')
            set(plotdata.Substrate.aWrite,'String','');
            set(plotdata.Substrate.bWrite,'String','');
            set(plotdata.Substrate.cWrite,'String','');
            set(plotdata.Substrate.alphaWrite,'String',['alpha=','','°']);
            set(plotdata.Substrate.betaWrite,'String',['beta=','','°']);
            set(plotdata.Substrate.gammaWrite,'String',['gamma=','','°']);
        else
        load('Substrates.mat')
        SubstrateTable=struct2table(SubstrateData);
        condition=find(strcmp(SubstrateTable.Name,plotdata.Substrate.Type));
        for k=1:size(condition),
            SubstrateData(condition(k)).b=b;
        end;
        if strcmp(SubstrateData(condition(1)).Structure,'Cubic'),
            for k=1:size(condition),
                SubstrateData(condition(k)).a=b;
                SubstrateData(condition(k)).c=b;
            end
            set(plotdata.Substrate.aWrite,'String',num2str(b));
            set(plotdata.Substrate.cWrite,'String',num2str(b));
        end,
        save('Substrates.mat','SubstrateData','-append');
        Substrate
        ProgcNSimu
        end
    
    case {'writeSubstrate_c'}
        c=str2num(get(plotdata.Substrate.cWrite,'string'));        
        if strcmp(plotdata.Substrate.Type,'none')
            set(plotdata.Substrate.aWrite,'String','');
            set(plotdata.Substrate.bWrite,'String','');
            set(plotdata.Substrate.cWrite,'String','');
            set(plotdata.Substrate.alphaWrite,'String',['alpha=','','°']);
            set(plotdata.Substrate.betaWrite,'String',['beta=','','°']);
            set(plotdata.Substrate.gammaWrite,'String',['gamma=','','°']);
        else
        load('Substrates.mat')
        SubstrateTable=struct2table(SubstrateData);
        condition=find(strcmp(SubstrateTable.Name,plotdata.Substrate.Type));
        for k=1:size(condition),
            SubstrateData(condition(k)).c=c;
        end;
        if strcmp(SubstrateData(condition(1)).Structure,'Cubic'),
            for k=1:size(condition),
                SubstrateData(condition(k)).a=c;
                SubstrateData(condition(k)).b=c;
            end
            set(plotdata.Substrate.aWrite,'String',num2str(c));
            set(plotdata.Substrate.bWrite,'String',num2str(c));
        end,
        save('Substrates.mat','SubstrateData','-append');
        Substrate
        ProgcNSimu
        end
        
    case {'writeNSubstrate'}
   plotdata.Substrate.N=str2num(get(plotdata.Substrate.NWrite,'string'));
   if plotdata.Substrate.N>plotdata.Substrate.NMax,
       plotdata.Substrate.N=plotdata.Substrate.NMax;
       set(plotdata.Substrate.NWrite,'String',plotdata.Substrate.N);
   elseif plotdata.Substrate.N<plotdata.Substrate.NMin,
       plotdata.Substrate.N=plotdata.Substrate.NMin;
       set(plotdata.Substrate.NWrite,'String',plotdata.Substrate.N);
   end
   Substrate
   ProgcNSimu
   
       case {'writeEpsilon'}
   plotdata.epsilon=str2num(get(plotdata.epsilonWrite,'string'));
   if plotdata.epsilon>plotdata.epsilonMax,
       plotdata.epsilon=plotdata.epsilonMax;
       set(plotdata.espilonWrite,'String',plotdata.epsilon);
   elseif plotdata.epsilon<plotdata.epsilonMin,
       plotdata.epsilon=plotdata.epsilonMin;
       set(plotdata.epsilonWrite,'String',plotdata.epsilon);
   end
   Substrate
   ProgcNSimu
        
case {'slidePolarization(1)','slidePolarization(2)','slidePolarization(3)','slidePolarization(4)','slidePolarization(5)','slidePolarization(6)'}
   k=str2num(entry(19));
   plotdata.Material(k).Polarization=get(plotdata.Material(k).PolarizationSlide,'Value');
   set(plotdata.Material(k).PolarizationWrite,'String',plotdata.Material(k).Polarization);
   set(plotdata.Material(k).TBPolarization,'SelectedObject',plotdata.Material(k).TBPolarizationYes);
   ProgcNSimu
%%     
case {'writePolarization(1)','writePolarization(2)','writePolarization(3)','writePolarization(4)','writePolarization(5)','writePolarization(6)'}
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
case {'slidexBST(1)','slidexBST(2)','slidexBST(3)','slidexBST(4)','slidexBST(5)','slidexBST(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xBST=get(plotdata.Material(k).xBSTSlide,'Value');
   set(plotdata.Material(k).xBSTWrite,'String',plotdata.Material(k).xBST);
   ProgcNSimu
%%     
case {'writexBST(1)','writexBST(2)','writexBST(3)','writexBST(4)','writexBST(5)','writexBST(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xBST=str2num(get(plotdata.Material(k).xBSTWrite,'string'));
   if plotdata.Material(k).xBST>plotdata.Material(k).xBSTMax,
       plotdata.Material(k).xBST=plotdata.Material(k).xBSTMax;
       set(plotdata.Material(k).xBSTWrite,'String',plotdata.Material(k).xBST);
   elseif plotdata.Material(k).xBST<plotdata.Material(k).xBSTMin,
       plotdata.Material(k).xBST=plotdata.Material(k).xBSTMin;
       set(plotdata.Material(k).xBSTWrite,'String',plotdata.Material(k).xBST);
   end;
   set(plotdata.Material(k).xBSTSlide,'value',plotdata.Material(k).xBST);
   ProgcNSimu
%%  
case {'slidexCSA(1)','slidexCSA(2)','slidexCSA(3)','slidexCSA(4)','slidexCSA(5)','slidexCSA(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xCSA=get(plotdata.Material(k).xCSASlide,'Value');
   set(plotdata.Material(k).xCSAWrite,'String',plotdata.Material(k).xCSA);
   ProgcNSimu
%%     
case {'writexCSA(1)','writexCSA(2)','writexCSA(3)','writexCSA(4)','writexCSA(5)','writexCSA(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xCSA=str2num(get(plotdata.Material(k).xCSAWrite,'string'));
   if plotdata.Material(k).xCSA>plotdata.Material(k).xCSAMax,
       plotdata.Material(k).xCSA=plotdata.Material(k).xCSAMax;
       set(plotdata.Material(k).xCSAWrite,'String',plotdata.Material(k).xCSA);
   elseif plotdata.Material(k).xCSA<plotdata.Material(k).xCSAMin,
       plotdata.Material(k).xCSA=plotdata.Material(k).xCSAMin;
       set(plotdata.Material(k).xCSAWrite,'String',plotdata.Material(k).xCSA);
   end;
   set(plotdata.Material(k).xCSASlide,'value',plotdata.Material(k).xCSA);
   ProgcNSimu
%%  
case {'slidexHZO(1)','slidexHZO(2)','slidexHZO(3)','slidexHZO(4)','slidexHZO(5)','slidexHZO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xHZO=get(plotdata.Material(k).xHZOSlide,'Value');
   set(plotdata.Material(k).xHZOWrite,'String',plotdata.Material(k).xHZO);
   ProgcNSimu
   %%  
case {'writexHZO(1)','writexHZO(2)','writexHZO(3)','writexHZO(4)','writexHZO(5)','writexHZO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xHZO=str2num(get(plotdata.Material(k).xHZOWrite,'string'));
   if plotdata.Material(k).xHZO>plotdata.Material(k).xHZOMax,
       plotdata.Material(k).xHZO=plotdata.Material(k).xHZOMax;
       set(plotdata.Material(k).xHZOWrite,'String',plotdata.Material(k).xHZO);
   elseif plotdata.Material(k).xHZO<plotdata.Material(k).xHZOMin,
       plotdata.Material(k).xHZO=plotdata.Material(k).xHZOMin;
       set(plotdata.Material(k).xHZOWrite,'String',plotdata.Material(k).xHZO);
   end;
   set(plotdata.Material(k).xHZOSlide,'value',plotdata.Material(k).xHZO);
   ProgcNSimu
%% 
case {'slidexLNO(1)','slidexLNO(2)','slidexLNO(3)','slidexLNO(4)','slidexLNO(5)','slidexLNO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xLNO=get(plotdata.Material(k).xLNOSlide,'Value');
   set(plotdata.Material(k).xLNOWrite,'String',plotdata.Material(k).xLNO);
   ProgcNSimu
%%     
case {'writexLNO(1)','writexLNO(2)','writexLNO(3)','writexLNO(4)','writexLNO(5)','writexLNO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xLNO=str2num(get(plotdata.Material(k).xLNOWrite,'string'));
   if plotdata.Material(k).xLNO>plotdata.Material(k).xLNOMax,
       plotdata.Material(k).xLNO=plotdata.Material(k).xLNOMax;
       set(plotdata.Material(k).xLNOWrite,'String',plotdata.Material(k).xLNO);
   elseif plotdata.Material(k).xLNO<plotdata.Material(k).xLNOMin,
       plotdata.Material(k).xLNO=plotdata.Material(k).xLNOMin;
       set(plotdata.Material(k).xLNOWrite,'String',plotdata.Material(k).xLNO);
   end;
   set(plotdata.Material(k).xLNOSlide,'value',plotdata.Material(k).xLNO);
   ProgcNSimu
%%  
case {'slideyLNO(1)','slideyLNO(2)','slideyLNO(3)','slideyLNO(4)','slideyLNO(5)','slideyLNO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).yLNO=get(plotdata.Material(k).yLNOSlide,'Value');
   set(plotdata.Material(k).yLNOWrite,'String',plotdata.Material(k).yLNO);
   ProgcNSimu
%%     
case {'writeyLNO(1)','writeyLNO(2)','writeyLNO(3)','writeyLNO(4)','writeyLNO(5)','writeyLNO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).yLNO=str2num(get(plotdata.Material(k).yLNOWrite,'string'));
   if plotdata.Material(k).yLNO>plotdata.Material(k).yLNOMax,
       plotdata.Material(k).yLNO=plotdata.Material(k).yLNOMax;
       set(plotdata.Material(k).yLNOWrite,'String',plotdata.Material(k).yLNO);
   elseif plotdata.Material(k).yLNO<plotdata.Material(k).yLNOMin,
       plotdata.Material(k).yLNO=plotdata.Material(k).yLNOMin;
       set(plotdata.Material(k).yLNOWrite,'String',plotdata.Material(k).yLNO);
   end;
   set(plotdata.Material(k).yLNOSlide,'value',plotdata.Material(k).yLNO);
   ProgcNSimu  
   %% 
case {'slidexLSCO(1)','slidexLSCO(2)','slidexLSCO(3)','slidexLSCO(4)','slidexLSCO(5)','slidexLSCO(6)'}
   k=str2num(entry(12));
   plotdata.Material(k).xLSCO=get(plotdata.Material(k).xLSCOSlide,'Value');
   set(plotdata.Material(k).xLSCOWrite,'String',plotdata.Material(k).xLSCO);
   ProgcNSimu    
%% 
case {'writexLSCO(1)','writexLSCO(2)','writexLSCO(3)','writexLSCO(4)','writexLSCO(5)','writexLSCO(5)'}
   k=str2num(entry(12));
   plotdata.Material(k).xLSCO=str2num(get(plotdata.Material(k).xLSCOWrite,'string'));
   if plotdata.Material(k).xLSCO>plotdata.Material(k).xLSCOMax,
       plotdata.Material(k).xLSCO=plotdata.Material(k).xLSCOMax;
       set(plotdata.Material(k).xLSCOWrite,'String',plotdata.Material(k).xLSCO);
   elseif plotdata.Material(k).xLSCO<plotdata.Material(k).xLSCOMin,
       plotdata.Material(k).xLSCO=plotdata.Material(k).xLSCOMin;
       set(plotdata.Material(k).xLSCOWrite,'String',plotdata.Material(k).xLSCO);
   end;
   set(plotdata.Material(k).xLSCOSlide,'value',plotdata.Material(k).xLSCO);
   ProgcNSimu
%%   
case {'slidexNCM(1)','slidexNCM(2)','slidexNCM(3)','slidexNCM(4)','slidexNCM(5)','slidexNCM(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xNCM=get(plotdata.Material(k).xNCMSlide,'Value');
   set(plotdata.Material(k).xNCMWrite,'String',plotdata.Material(k).xNCM);
   ProgcNSimu
   %%  
case {'writexNCM(1)','writexNCM(2)','writexNCM(3)','writexNCM(4)','writexNCM(5)','writexNCM(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xNCM=str2num(get(plotdata.Material(k).xNCMWrite,'string'));
   if plotdata.Material(k).xNCM>plotdata.Material(k).xNCMMax,
       plotdata.Material(k).xNCM=plotdata.Material(k).xNCMMax;
       set(plotdata.Material(k).xNCMWrite,'String',plotdata.Material(k).xNCM);
   elseif plotdata.Material(k).xNCM<plotdata.Material(k).xNCMMin,
       plotdata.Material(k).xNCM=plotdata.Material(k).xNCMMin;
       set(plotdata.Material(k).xNCMWrite,'String',plotdata.Material(k).xNCM);
   end;
   set(plotdata.Material(k).xNCMSlide,'value',plotdata.Material(k).xNCM);
   ProgcNSimu
%% 
case {'slidexNLN(1)','slidexNLN(2)','slidexNLN(3)','slidexNLN(4)','slidexNLN(5)','slidexNLN(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xNLN=get(plotdata.Material(k).xNLNSlide,'Value');
   set(plotdata.Material(k).xNLNWrite,'String',plotdata.Material(k).xNLN);
   ProgcNSimu
%%     
case {'writexNLN(1)','writexNLN(2)','writexNLN(3)','writexNLN(4)','writexNLN(5)','writexNLN(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xNLN=str2num(get(plotdata.Material(k).xNLNWrite,'string'));
   if plotdata.Material(k).xNLN>plotdata.Material(k).xNLNMax,
       plotdata.Material(k).xNLN=plotdata.Material(k).xNLNMax;
       set(plotdata.Material(k).xNLNWrite,'String',plotdata.Material(k).xNLN);
   elseif plotdata.Material(k).xNLN<plotdata.Material(k).xNLNMin,
       plotdata.Material(k).xNLN=plotdata.Material(k).xNLNMin;
       set(plotdata.Material(k).xNLNWrite,'String',plotdata.Material(k).xNLN);
   end;
   set(plotdata.Material(k).xNLNSlide,'value',plotdata.Material(k).xNLN);
   ProgcNSimu
%% 
case {'slidexNNO(1)','slidexNNO(2)','slidexNNO(3)','slidexNNO(4)','slidexNNO(5)','slidexNNO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xNNO=get(plotdata.Material(k).xNNOSlide,'Value');
   set(plotdata.Material(k).xNNOWrite,'String',plotdata.Material(k).xNNO);
   ProgcNSimu
%%     
case {'writexNNO(1)','writexNNO(2)','writexNNO(3)','writexNNO(4)','writexNNO(5)','writexNNO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xNNO=str2num(get(plotdata.Material(k).xNNOWrite,'string'));
   if plotdata.Material(k).xNNO>plotdata.Material(k).xNNOMax,
       plotdata.Material(k).xNNO=plotdata.Material(k).xNNOMax;
       set(plotdata.Material(k).xNNOWrite,'String',plotdata.Material(k).xNNO);
   elseif plotdata.Material(k).xNNO<plotdata.Material(k).xNNOMin,
       plotdata.Material(k).xNNO=plotdata.Material(k).xNNOMin;
       set(plotdata.Material(k).xNNOWrite,'String',plotdata.Material(k).xNNO);
   end;
   set(plotdata.Material(k).xNNOSlide,'value',plotdata.Material(k).xNNO);
   ProgcNSimu
%%  
case {'slideyNNO(1)','slideyNNO(2)','slideyNNO(3)','slideyNNO(4)','slideyNNO(5)','slideyNNO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).yNNO=get(plotdata.Material(k).yNNOSlide,'Value');
   set(plotdata.Material(k).yNNOWrite,'String',plotdata.Material(k).yNNO);
   ProgcNSimu
%%     
case {'writeyNNO(1)','writeyNNO(2)','writeyNNO(3)','writeyNNO(4)','writeyNNO(5)','writeyNNO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).yNNO=str2num(get(plotdata.Material(k).yNNOWrite,'string'));
   if plotdata.Material(k).yNNO>plotdata.Material(k).yNNOMax,
       plotdata.Material(k).yNNO=plotdata.Material(k).yNNOMax;
       set(plotdata.Material(k).yNNOWrite,'String',plotdata.Material(k).yNNO);
   elseif plotdata.Material(k).yNNO<plotdata.Material(k).yNNOMin,
       plotdata.Material(k).yNNO=plotdata.Material(k).yNNOMin;
       set(plotdata.Material(k).yNNOWrite,'String',plotdata.Material(k).yNNO);
   end;
   set(plotdata.Material(k).yNNOSlide,'value',plotdata.Material(k).yNNO);
   ProgcNSimu
%%  
case {'slidexNSM(1)','slidexNSM(2)','slidexNSM(3)','slidexNSM(4)','slidexNSM(5)','slidexNSM(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xNSM=get(plotdata.Material(k).xNSMSlide,'Value');
   set(plotdata.Material(k).xNSMWrite,'String',plotdata.Material(k).xNSM);
   ProgcNSimu
   %%  
case {'writexNSM(1)','writexNSM(2)','writexNSM(3)','writexNSM(4)','writexNSM(5)','writexNSM(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xNSM=str2num(get(plotdata.Material(k).xNSMWrite,'string'));
   if plotdata.Material(k).xNSM>plotdata.Material(k).xNSMMax,
       plotdata.Material(k).xNSM=plotdata.Material(k).xNSMMax;
       set(plotdata.Material(k).xNSMWrite,'String',plotdata.Material(k).xNSM);
   elseif plotdata.Material(k).xNSM<plotdata.Material(k).xNSMMin,
       plotdata.Material(k).xNSM=plotdata.Material(k).xNSMMin;
       set(plotdata.Material(k).xNSMWrite,'String',plotdata.Material(k).xNSM);
   end;
   set(plotdata.Material(k).xNSMSlide,'value',plotdata.Material(k).xNSM);
   ProgcNSimu
%% 
case {'slidexPST(1)','slidexPST(2)','slidexPST(3)','slidexPST(4)','slidexPST(5)','slidexPST(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xPST=get(plotdata.Material(k).xPSTSlide,'Value');
   set(plotdata.Material(k).xPSTWrite,'String',plotdata.Material(k).xPST);
   ProgcNSimu    
%% 
case {'writexPST(1)','writexPST(2)','writexPST(3)','writexPST(4)','writexPST(5)','writexPST(5)'}
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
case {'slidexPZT(1)','slidexPZT(2)','slidexPZT(3)','slidexPZT(4)','slidexPZT(5)','slidexPZT(5)'}
   k=str2num(entry(11));
   plotdata.Material(k).xPZT=get(plotdata.Material(k).xPZTSlide,'Value');
   set(plotdata.Material(k).xPZTWrite,'String',plotdata.Material(k).xPZT);
   ProgcNSimu
%%     
case {'writexPZT(1)','writexPZT(2)','writexPZT(3)','writexPZT(4)','writexPZT(5)','writexPZT(5)'}
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
case {'slidexYTmIG(1)','slidexYTmIG(2)','slidexYTmIG(3)','slidexYTmIG(4)','slidexYTmIG(5)','slidexYTmIG(5)'}
   k=str2num(entry(13));
   plotdata.Material(k).xYTmIG=get(plotdata.Material(k).xYTmIGSlide,'Value');
   set(plotdata.Material(k).xYTmIGWrite,'String',plotdata.Material(k).xYTmIG);
   ProgcNSimu
%%     
case {'writexYTmIG(1)','writexYTmIG(2)','writexYTmIG(3)','writexYTmIG(4)','writexYTmIG(5)','writexYTmIG(6)'}
   k=str2num(entry(13));
   plotdata.Material(k).xYTmIG=str2num(get(plotdata.Material(k).xYTmIGWrite,'string'));
   if plotdata.Material(k).xYTmIG>plotdata.Material(k).xYTmIGMax,
       plotdata.Material(k).xYTmIG=plotdata.Material(k).xYTmIGMax;
       set(plotdata.Material(k).xYTmIGWrite,'String',plotdata.Material(k).xYTmIG);
   elseif plotdata.Material(k).xYTmIG<plotdata.Material(k).xYTmIGMin,
       plotdata.Material(k).xYTmIG=plotdata.Material(k).xYTmIGMin;
       set(plotdata.Material(k).xYTmIGWrite,'String',plotdata.Material(k).xYTmIG);
   end;
   set(plotdata.Material(k).xYTmIGSlide,'value',plotdata.Material(k).xYTmIG);
   ProgcNSimu
%%     
case {'slidexZMO(1)','slidexZMO(2)','slidexZMO(3)','slidexZMO(4)','slidexZMO(5)','slidexZMO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xZMO=get(plotdata.Material(k).xZMOSlide,'Value');
   set(plotdata.Material(k).xZMOWrite,'String',plotdata.Material(k).xZMO);
   ProgcNSimu
%%     
case {'writexZMO(1)','writexZMO(2)','writexZMO(3)','writexZMO(4)','writexZMO(5)','writexZMO(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).xZMO=str2num(get(plotdata.Material(k).xZMOWrite,'string'));
   if plotdata.Material(k).xZMO>plotdata.Material(k).xZMOMax,
       plotdata.Material(k).xZMO=plotdata.Material(k).xZMOMax;
       set(plotdata.Material(k).xZMOWrite,'String',plotdata.Material(k).xZMO);
   elseif plotdata.Material(k).xZMO<plotdata.Material(k).xZMOMin,
       plotdata.Material(k).xZMO=plotdata.Material(k).xZMOMin;
       set(plotdata.Material(k).xZMOWrite,'String',plotdata.Material(k).xZMO);
   end;
   set(plotdata.Material(k).xZMOSlide,'value',plotdata.Material(k).xZMO);
   ProgcNSimu
%%
case {'slideN(1)','slideN(2)','slideN(3)','slideN(4)','slideN(5)','slideN(6)'}
   k=str2num(entry(8));
   plotdata.Material(k).N=get(plotdata.Material(k).NSlide,'Value');
   set(plotdata.Material(k).NWrite,'String',plotdata.Material(k).N);
   ProgcNSimu
%%
case {'writeN(1)','writeN(2)','writeN(3)','writeN(4)','writeN(5)','writeN(6)'}
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
case {'slided(1)','slided(2)','slided(3)','slided(4)','slided(5)','slided(6)'}
   k=str2num(entry(8)); 
   plotdata.Material(k).d=get(plotdata.Material(k).dSlide,'Value');
   set(plotdata.Material(k).dWrite,'String',plotdata.Material(k).d);
   ProgcNSimu
%%  
case {'writed(1)','writed(2)','writed(3)','writed(4)','writed(5)','writed(6)'}
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
case {'slidedCenter(1)','slidedCenter(2)','slidedCenter(3)','slidedCenter(4)','slidedCenter(5)','slidedCenter(6)'}
   k=str2num(entry(14)); 
   plotdata.Material(k).dCenter=get(plotdata.Material(k).dCenterSlide,'Value');
   set(plotdata.Material(k).dCenterWrite,'String',plotdata.Material(k).dCenter);
   ProgcNSimu
%%  
case {'writedCenter(1)','writedCenter(2)','writedCenter(3)','writedCenter(4)','writedCenter(5)','writedCenter(6)'}
   k=str2num(entry(14)); 
   plotdata.Material(k).dCenter=str2num(get(plotdata.Material(k).dCenterWrite,'string'));
   if plotdata.Material(k).dCenter>plotdata.Material(k).dMax,
       plotdata.Material(k).dCenter=plotdata.Material(k).dMax;
       set(plotdata.Material(k).dCenterWrite,'String',plotdata.Material(k).dCenter);
   elseif plotdata.Material(k).dCenter<plotdata.Material(k).dMin,
       plotdata.Material(k).dCenter=plotdata.Material(k).dMin;
       set(plotdata.Material(k).dCenterWrite,'String',plotdata.Material(k).dCenter);
   end;
   set(plotdata.Material(k).dCenterSlide,'value',plotdata.Material(k).dCenter);
   ProgcNSimu
%%
case {'slidedBorder(1)','slidedBorder(2)','slidedBorder(3)','slidedBorder(4)','slidedBorder(5)','slidedBorder(6)'}
   k=str2num(entry(14)); 
   plotdata.Material(k).dBorder=get(plotdata.Material(k).dBorderSlide,'Value');
   set(plotdata.Material(k).dBorderWrite,'String',plotdata.Material(k).dBorder);
   ProgcNSimu
%%  
case {'writedBorder(1)','writedBorder(2)','writedBorder(3)','writedBorder(4)','writedBorder(5)','writedBorder(6)'}
   k=str2num(entry(14)); 
   plotdata.Material(k).dBorder=str2num(get(plotdata.Material(k).dBorderWrite,'string'));
   if plotdata.Material(k).dBorder>plotdata.Material(k).dMax,
       plotdata.Material(k).dBorder=plotdata.Material(k).dMax;
       set(plotdata.Material(k).dBorderWrite,'String',plotdata.Material(k).dBorder);
   elseif plotdata.Material(k).dBorder<plotdata.Material(k).dMin,
       plotdata.Material(k).dBorder=plotdata.Material(k).dMin;
       set(plotdata.Material(k).dBorderWrite,'String',plotdata.Material(k).dBorder);
   end;
   set(plotdata.Material(k).dBorderSlide,'value',plotdata.Material(k).dBorder);
   ProgcNSimu
%%
case {'slideAexp(1)','slideAexp(2)','slideAexp(3)','slideAexp(4)','slideAexp(5)','slideAexp(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).Aexp=get(plotdata.Material(k).AexpSlide,'Value');
   set(plotdata.Material(k).AexpWrite,'String',plotdata.Material(k).Aexp);
   ProgcNSimu
   k=str2num(entry(11)); 
   set(plotdata.Material(k).CexpWrite,'String',plotdata.Material(k).Cexp);
%%   
case {'writeAexp(1)','writeAexp(2)','writeAexp(3)','writeAexp(4)','writeAexp(5)','writeAexp(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).Aexp=str2num(get(plotdata.Material(k).AexpWrite,'string'));
   if plotdata.Material(k).Aexp>plotdata.Material(k).AexpMax,
       plotdata.Material(k).Aexp=plotdata.Material(k).AexpMax;
       set(plotdata.Material(k).AexpWrite,'String',plotdata.Material(k).Aexp);
   elseif plotdata.Material(k).Aexp<plotdata.Material(k).AexpMin,
       plotdata.Material(k).Aexp=plotdata.Material(k).AexpMin;
       set(plotdata.Material(k).AexpWrite,'String',plotdata.Material(k).Aexp);
   end;
   set(plotdata.Material(k).AexpSlide,'value',plotdata.Material(k).Aexp);
   ProgcNSimu
   k=str2num(entry(11));
   set(plotdata.Material(k).CexpWrite,'String',num2str(plotdata.Material(k).Cexp));
%%   
case {'slideBexp(1)','slideBexp(2)','slideBexp(3)','slideBexp(4)','slideBexp(5)','slideBexp(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).Bexp=get(plotdata.Material(k).BexpSlide,'Value');
   set(plotdata.Material(k).BexpWrite,'String',plotdata.Material(k).Bexp);
   ProgcNSimu
   k=str2num(entry(11));
   set(plotdata.Material(k).CexpWrite,'String',num2str(plotdata.Material(k).Cexp));
%%   
case {'writeBexp(1)','writeBexp(2)','writeBexp(3)','writeBexp(4)','writeBexp(5)','writeBexp(6)'}
   k=str2num(entry(11));
   plotdata.Material(k).Bexp=str2num(get(plotdata.Material(k).BexpWrite,'string'));
   if plotdata.Material(k).Bexp>plotdata.Material(k).BexpMax,
       plotdata.Material(k).Bexp=plotdata.Material(k).BexpMax;
       set(plotdata.Material(k).BexpWrite,'String',plotdata.Material(k).Bexp);
   elseif plotdata.Material(k).Bexp<plotdata.Material(k).BexpMin,
       plotdata.Material(k).Bexp=plotdata.Material(k).BexpMin;
       set(plotdata.Material(k).BexpWrite,'String',plotdata.Material(k).Bexp);
   end;
   set(plotdata.Material(k).BexpSlide,'value',plotdata.Material(k).Bexp);
   ProgcNSimu
   k=str2num(entry(11));
   set(plotdata.Material(k).CexpWrite,'String',num2str(plotdata.Material(k).Cexp));
%%   
case {'slideMeandexp(1)','slideMeandexp(2)','slideMeandexp(3)','slideMeandexp(4)','slideMeandexp(5)','slideMeandexp(6)'}
   k=str2num(entry(15));
   plotdata.Material(k).Meandexp=get(plotdata.Material(k).MeandexpSlide,'Value');
   set(plotdata.Material(k).MeandexpWrite,'String',plotdata.Material(k).Meandexp);
   ProgcNSimu
   k=str2num(entry(15));
   set(plotdata.Material(k).CexpWrite,'String',num2str(plotdata.Material(k).Cexp));
%%   
case {'writeMeandexp(1)','writeMeandexp(2)','writeMeandexp(3)','writeMeandexp(4)','writeMeandexp(5)','writeMeandexp(6)'}
   k=str2num(entry(15));
   plotdata.Material(k).Meandexp=str2num(get(plotdata.Material(k).MeandexpWrite,'string'));
   if plotdata.Material(k).Meandexp>plotdata.Material(k).MeandexpMax,
       plotdata.Material(k).Meandexp=plotdata.Material(k).MeandexpMax;
       set(plotdata.Material(k).MeandexpWrite,'String',plotdata.Material(k).Meandexp);
   elseif plotdata.Material(k).Meandexp<plotdata.Material(k).MeandexpMin,
       plotdata.Material(k).Meandexp=plotdata.Material(k).MeandexpMin;
       set(plotdata.Material(k).MeandexpWrite,'String',plotdata.Material(k).Meandexp);
   end;
   set(plotdata.Material(k).MeandexpSlide,'value',plotdata.Material(k).Meandexp);
   ProgcNSimu
   k=str2num(entry(15));
   set(plotdata.Material(k).CexpWrite,'String',num2str(plotdata.Material(k).Cexp));
%%   
case {'slidem(1)','slidem(2)','slidem(3)','slidem(4)','slidem(5)','slidem(6)'}
   k=str2num(entry(8)); 
   plotdata.Material(k).m=get(plotdata.Material(k).mSlide,'Value');
   set(plotdata.Material(k).mWrite,'String',plotdata.Material(k).m);
   ProgcNSimu
%%  
case {'writem(1)','writem(2)','writem(3)','writem(4)','writem(5)','writem(6)'}
   k=str2num(entry(8)); 
   plotdata.Material(k).m=str2num(get(plotdata.Material(k).mWrite,'string'));
   if plotdata.Material(k).m>plotdata.mMax,
       plotdata.Material(k).m=plotdata.mMax;
       set(plotdata.Material(k).mWrite,'String',plotdata.Material(k).m);
   elseif plotdata.Material(k).m<plotdata.mMin,
       plotdata.Material(k).m=plotdata.mMin;
       set(plotdata.Material(k).mWrite,'String',plotdata.Material(k).m);
   end;
   set(plotdata.Material(k).mSlide,'value',plotdata.Material(k).m);
   ProgcNSimu
%%
case 'slideRepetition'
   plotdata.Repetition.N=get(plotdata.RepetitionSlide,'Value');
   set(plotdata.RepetitionWrite,'String',plotdata.Repetition.N);
   ProgcNSimu
%%    
case 'writeRepetition'
   plotdata.Repetition.N=str2num(get(plotdata.RepetitionWrite,'string'));
   if plotdata.Repetition.N>plotdata.Repetition.Max,
       plotdata.Repetition.N=plotdata.Repetition.Max;
       set(plotdata.RepetitionWrite,'String',plotdata.Repetition.N);
   elseif plotdata.Repetition.N<plotdata.Repetition.Min,
       plotdata.Repetition.N=plotdata.Repetition.Min;
       set(plotdata.RepetitionWrite,'String',plotdata.Repetition.N);
   end;
   set(plotdata.RepetitionSlide,'value',plotdata.Repetition.N);
   ProgcNSimu
%%   
case {'Loadd(1)','Loadd(2)','Loadd(3)','Loadd(4)','Loadd(5)','Loadd(6)'} 
   k=str2num(entry(7));
   if exist('InteractiveXRDFitd.mat', 'file')
       load('InteractiveXRDFitd.mat','loadfolderd'); 
       if loadfolderd==0;
           loadfolderd=[];
       end;
       [plotdata.measure.dfilename, loadfolderd] = uigetfile({'*.csv';'*.*'},'Select the d you wand to upload (.csv)',loadfolderd);
   else
       [plotdata.measure.dfilename, loadfolderd] = uigetfile({'*.csv';'*.*'},'Select the d you wand to upload (.csv)');
   end
   if loadfolderd~=0;
    save('InteractiveXRDFitd.mat','loadfolderd');
    %plotdata.prog=[loadfolderd,plotdata.measure.dfilename];
    %set(plotdata.filenameWrite,'String',plotdata.prog);
    %open and load data
    fid=fopen([loadfolderd,plotdata.measure.dfilename],'r'); %open file [] with permission r = open file for reading
    data=[];
    %idxAngle=[];
    %idxIntensity=[];
    if (fid==-1) % if fopen cannot open the file
        msgbox(['ERROR: file not found or could not be opened for read. check file ' [loadfolderd,plotdata.measure.dfilename]]); 
        return;
    end
    %% skip first line of comments
    while 1 
        sLine=fgetl(fid); %read file one line at a time, removing newline characters
        if (~ischar(sLine))
            % end of file 
            break; 
        end
        % add numerical lines to table
        data=[data;str2num(sLine)];
    end 
    fclose(fid); 
    plotdata.Material(k).dupload=data(:,1);    
    plotdata.Material(k).N=size(data,1);
    set(plotdata.Material(k).NWrite,'String',plotdata.Material(k).N);
    set(plotdata.Material(k).NSlide,'value',plotdata.Material(k).N);
    ProgcNSimu

   end;
        
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
    fid=fopen([plotdata.prog],'r'); %open file [plotdata.prog] with permission r = open file for reading
    data=[];
    idxAngle=[];
    idxIntensity=[];
    if (fid==-1) % if fopen cannot open the file
        msgbox(['ERROR: file not found or could not be opened for read. check file ' [plotdata.prog]]); 
        return;
    end
    %% skip first line of comments and find Angle and Intensity columns if specified
    while 1 
        sLine=fgetl(fid); %read file one line at a time, removing newline characters
        if (~ischar(sLine))
            % end of file 
            break; 
        end
        % find line for "Angle" and "Intensity", if specified
        tfAngle=regexp(sLine,'Angle');
        tfIntensity=regexp(sLine,'Intensity');
        if (~isempty(tfAngle)),
            header=strsplit(sLine,',');
            idxAngle=find(strcmp('Angle',header));
        end;
        if (~isempty(tfIntensity)),
            header=strsplit(sLine,',');
            idxIntensity=find(strcmp('Intensity',header));
            if isempty(idxIntensity),
                idxIntensity=find(strcmp(' Intensity',header));
            end;
        end;
        % add numerical lines to table
        data=[data;str2num(sLine)];
    end 
    fclose(fid); 
    if isempty(idxAngle),
        plotdata.measure.x=data(:,1); 
    else plotdata.measure.x=data(:,idxAngle);
    end;
    if isempty(idxIntensity),
        plotdata.measure.y=data(:,2);
    else plotdata.measure.y=data(:,idxIntensity);
    end;
    plotdata.measure.y=plotdata.measure.y/max(plotdata.measure.y);
    
    PlotFitAndData
    figure(plotdata.FitDisplay.fig);
   end;
          
case 'chooseSubstrate'
   load('Substrates.mat');
   str=get(plotdata.Substrate.choose,'String');
   entry=get(plotdata.Substrate.choose,'value');
   SubstrateType=strsplit(str(entry,:));
   plotdata.Substrate.Type=SubstrateType(1);
   SubstrateTable=struct2table(SubstrateData);
   condition=find(strcmp(SubstrateTable.Name,plotdata.Substrate.Type));
   plotdata.Substrate.Orientation=SubstrateData(condition(1)).Orientation;
   set([plotdata.Substrate.TBorientation.Cubic,...
       plotdata.Substrate.TBorientation.Hexagonal,...
       plotdata.Substrate.TBorientation.Monoclinic,...
       plotdata.Substrate.TBorientation.Orthorhombic,...
       plotdata.Substrate.TBorientation.PseudoCubic,...
       plotdata.Substrate.TBorientation.Tetragonal],...
       'visible','off');
   set(plotdata.Substrate.StructureWrite,'String',SubstrateData(condition(1)).Structure);
   set(plotdata.Substrate.aWrite,'String',num2str(SubstrateData(condition(1)).a));
   set(plotdata.Substrate.bWrite,'String',num2str(SubstrateData(condition(1)).b));
   set(plotdata.Substrate.cWrite,'String',num2str(SubstrateData(condition(1)).c));
   set(plotdata.Substrate.alphaWrite,'String',['alpha=',num2str(SubstrateData(condition(1)).alpha),'°']);
   set(plotdata.Substrate.betaWrite,'String',['beta=',num2str(SubstrateData(condition(1)).beta),'°']);
   set(plotdata.Substrate.gammaWrite,'String',['gamma=',num2str(SubstrateData(condition(1)).gamma),'°']);
   switch char(SubstrateData(condition(1)).Structure)
       case 'Cubic'
           set(plotdata.Substrate.TBorientation.Cubic,'visible','on');
           plotdata.Substrate.Orientation='(001)c';
           plotdata.Substrate.TBc001.Value=true;
       case 'Hexagonal'
           set(plotdata.Substrate.TBorientation.Hexagonal,'visible','on');
           plotdata.Substrate.Orientation='(0001)h';
           plotdata.Substrate.TBh0001.Value=true;
       case 'Monoclinic'
           set(plotdata.Substrate.TBorientation.Monoclinic,'visible','on');
           plotdata.Substrate.Orientation='(010)m';
           plotdata.Substrate.TBm010.Value=true;
       case 'Orthorhombic'
           set(plotdata.Substrate.TBorientation.Orthorhombic,'visible','on');
           plotdata.Substrate.Orientation='(001)o';
           plotdata.Substrate.TBo001.Value=true;
           if and(strcmp(char(plotdata.Substrate.Type),'NdGaO3'),strcmp(char(plotdata.Substrate.Orientation),'(001)o')),
            SubstrateWarning(plotdata.Substrate.Orientation,plotdata.Substrate.Type);
            plotdata.Substrate.Orientation='(110)o';
            plotdata.Substrate.TBo110.Value=true;
           end
       case 'PseudoCubic'
           set(plotdata.Substrate.TBorientation.PseudoCubic,'visible','on');
           plotdata.Substrate.Orientation='(001)pc';
           plotdata.Substrate.TBpc001.Value=true;
       case 'Tetragonal'
           set(plotdata.Substrate.TBorientation.Tetragonal,'visible','on');
           plotdata.Substrate.Orientation='(001)t';
           plotdata.Substrate.TBt001.Value=true;
   end;
   Substrate
   ProgcNSimu
%%   
case {'chooseLayerType(1)','chooseLayerType(2)','chooseLayerType(3)','chooseLayerType(4)','chooseLayerType(5)','chooseLayerType(6)'}
   k=str2num(entry(17));
   str=get(plotdata.Material(k).choose,'String');
   entry=get(plotdata.Material(k).choose,'Value');
   MaterialType=strsplit(str(entry,:));
   plotdata.Material(k).Type=MaterialType(1);
   plotdata.Material(k).Polarization=0;
   set(plotdata.Material(k).xBSTpanel,'Visible','off');
   set(plotdata.Material(k).xCSApanel,'Visible','off');
   set(plotdata.Material(k).xHZOpanel,'Visible','off');
   set(plotdata.Material(k).xLSCOpanel,'Visible','off');
   set(plotdata.Material(k).xNCMpanel,'Visible','off');
   set(plotdata.Material(k).xNLNpanel,'Visible','off');
   set(plotdata.Material(k).xNNOpanel,'Visible','off');
   set(plotdata.Material(k).yNNOpanel,'Visible','off');
   set(plotdata.Material(k).xNSMpanel,'Visible','off');
   set(plotdata.Material(k).xPSTpanel,'Visible','off');
   set(plotdata.Material(k).xPZTpanel,'Visible','off');
   set(plotdata.Material(k).xYTmIGpanel,'Visible','off');
   set(plotdata.Material(k).xZMOpanel,'Visible','off');
   set(plotdata.Material(k).TBPolarization,'Visible','off');
   set(plotdata.Material(k).OrientationPanel,'Visible','off');
   set(plotdata.Material(k).choosepcOrientation,'Visible','off');
   set(plotdata.Material(k).choosehOrientation,'Visible','off');
 
   if strcmp(plotdata.Material(k).Type,'none'),
       plotdata.Material(k).N=0; plotdata.Material(k).d=0;
       plotdata.Material(k).Orientation='none';
   
   elseif any(strcmp(plotdata.Material(k).Type,{'AlO2','BaO','LaO','MnO','MnO2','NiO2','NdO','PbO','RuO2','SrO','SrO2','TiO2','VO2'})),
       %AO-type or BO2-type monolayer CAREFULL: THIS IS ABSOLUTELY "ARTIFICIAL"
       plotdata.Material(k).d=2; %distance initiale entre 2 plans
       
   elseif strcmp(plotdata.Material(k).Type,'BaBiO3'),
       %perovskite
       plotdata.Material(k).dpc001=4.43647; %initial distance between two planes
       plotdata.Material(k).dpc111=4.43647/sqrt(3); %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       
   elseif strcmp(plotdata.Material(k).Type,'BaSnO3'),
       %perovskite
       plotdata.Material(k).dpc001=4.188634; %initial distance between two planes
       plotdata.Material(k).dpc111=4.188634/sqrt(3); %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       
   elseif strcmp(plotdata.Material(k).Type,'(Ba_x,Sr_{1-x})TiO3'),
       %perovskite
       plotdata.Material(k).dpc001=3.905; %initial distance between two planes (default x=0 i.e. SrTiO3)
       plotdata.Material(k).dpc111=3.905/sqrt(3); %initial distance between two planes (default x=0 i.e. SrTiO3)
       plotdata.Material(k).d=plotdata.Material(k).dpc001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xBSTpanel,'Visible','on');
 
    elseif strcmp(plotdata.Material(k).Type,'BaTiO3'),
       %perovskite
       plotdata.Material(k).dpc001=4.036; %initial distance between two planes
       plotdata.Material(k).dpc111=4.036/sqrt(3); %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
    
    elseif strcmp(plotdata.Material(k).Type,'BiFeO3'),
       %perovskite
       plotdata.Material(k).dpc001=4.924; %initial distance between two planes
       plotdata.Material(k).dpc111=4.924/sqrt(3); %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
   
   elseif strcmp(plotdata.Material(k).Type,'Bi2Te3'),
       % Bi2Te3: hexagonal with primitive unit cell Bi6Te9 a=b=4.44949, c=31.150873, alpha=beta=90?, gamma=120?
       plotdata.Material(k).dh0001=31.150873; %initial distance between 2 planes
       plotdata.Material(k).dh11bar20=4.44949/2;
       plotdata.Material(k).dh10bar10=4.44949*sqrt(3)/2;
       plotdata.Material(k).d=plotdata.Material(k).dh0001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosehOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='h(0001)';       
       
    elseif strcmp(plotdata.Material(k).Type,{'CaCuO2'}),
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       plotdata.Material(k).d=3.20546; %initial distance between two planes
    
    elseif strcmp(plotdata.Material(k).Type,{'Ca2RuO4'}),
       %pseudotetragonal
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       plotdata.Material(k).d=11.968145; %initial distance between two planes
          
   elseif strcmp(plotdata.Material(k).Type,'(Ca_{3-x},Sr_x)Al2O6'),
       %perovskite
       plotdata.Material(k).d=7.999957; %initial distance between two planes (default x=3 i.e. Sr3Al2O6) https://materialsproject.org/materials/mp-3393/
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       set(plotdata.Material(k).xCSApanel,'Visible','on');
       
    elseif strcmp(plotdata.Material(k).Type,'CaTiO3'),
       %perovskite
       plotdata.Material(k).dpc001=3.889; %initial distance between two planes
       plotdata.Material(k).dpc111=3.889/sqrt(3); %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       
    elseif strcmp(plotdata.Material(k).Type,'CaVO3'),
       %perovskite
       plotdata.Material(k).dpc001=3.830; %initial distance between two planes
       plotdata.Material(k).dpc111=3.830/sqrt(3); %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
          
    elseif strcmp(plotdata.Material(k).Type,'CoO'),
       plotdata.Material(k).dpc001=4.263; %initial distance between two planes
       plotdata.Material(k).dpc111=4.263/sqrt(3); %initial distance between two planes 
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       
    elseif strcmp(plotdata.Material(k).Type,'mu-Fe2O3'),
       plotdata.Material(k).d=3.16; %initial distance between 2 (010)m planes   
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosemOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='m(010)';
       
    elseif strcmp(plotdata.Material(k).Type,'beta-Ga2O3'),
       plotdata.Material(k).d=3.0371; %initial distance between 2 (010)m planes   
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosemOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='m(010)';
   
   elseif strcmp(plotdata.Material(k).Type,'GeTe'),
       % GeTe: hexagonal with primitive unit cell Ge3Te3 a=b= 4.23067, c= 10.88958, alpha=beta=90, gamma=120
       plotdata.Material(k).dh0001=10.88958; %initial distance between 2 planes
       plotdata.Material(k).dh11bar20=4.23067/2;
       plotdata.Material(k).dh10bar10=4.23067*sqrt(3)/2;
       plotdata.Material(k).d=plotdata.Material(k).dh0001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosehOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='h(0001)'; 
       
    elseif strcmp(plotdata.Material(k).Type,'(Hf_{1-x},Zr_x)O2'),
       plotdata.Material(k).dpc001=5.1; %initial distance between 2 planes - default x=0: HfO2
       plotdata.Material(k).dpc111=5.1/sqrt(3)*6/2; %initial distance between 2 planes - default x=0: HfO2
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xHZOpanel,'Visible','on'); 
       
    elseif strcmp(plotdata.Material(k).Type,'In2O3'),
       plotdata.Material(k).dpc001=10.29956; %initial distance between two planes        
       plotdata.Material(k).dpc111=10.29956/sqrt(3); %initial distance between two planes 
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       
    elseif strcmp(plotdata.Material(k).Type,'LaAlO3'),
       plotdata.Material(k).dpc001=3.787; %initial distance between two planes        
       plotdata.Material(k).dpc111=3.787/sqrt(3); %initial distance between two planes 
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
 
    elseif strcmp(plotdata.Material(k).Type,'LaCoO3'),
       plotdata.Material(k).dpc001=3.816; %initial distance between two planes
       plotdata.Material(k).dpc111=3.816/sqrt(3); %initial distance between two planes 
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       
     elseif strcmp(plotdata.Material(k).Type,'LaCrO3'),
       plotdata.Material(k).dpc001=3.885; %initial distance between two planes
       plotdata.Material(k).dpc111=3.885/sqrt(3); %initial distance between two planes 
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
    
    elseif strcmp(plotdata.Material(k).Type,'La2CuO4'), %Ruddlesden-Popper phase
       plotdata.Material(k).d=5.398; %initial distance between two planes  
       set(plotdata.Material(k).OrientationPanel,'visible','off');
 
    elseif strcmp(plotdata.Material(k).Type,'LaFeO3'),
       plotdata.Material(k).dpc001=3.959; %initial distance between two planes 
       plotdata.Material(k).dpc111=3.959/sqrt(3); %initial distance between two planes 
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';   
       
    elseif strcmp(plotdata.Material(k).Type,'LaMnO3'),
       plotdata.Material(k).dpc001=3.945; %initial distance between two planes
       plotdata.Material(k).dpc111=3.945/sqrt(3); %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';         
  
    elseif strcmp(plotdata.Material(k).Type,'LaNiO3'),
       plotdata.Material(k).dpc001=3.857; %initial distance between two planes
       plotdata.Material(k).dpc111=3.857/sqrt(3); %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';  
    
    elseif strcmp(plotdata.Material(k).Type,'La_xNi_yO3'),%off-stoichiometry
       plotdata.Material(k).dpc001=3.857; %initial distance between 2 planes
       plotdata.Material(k).dpc111=3.857/sqrt(3); %initial distance between 2 planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes       
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xLNOpanel,'Visible','on');
       set(plotdata.Material(k).yLNOpanel,'Visible','on');
       
    elseif strcmp(plotdata.Material(k).Type,'La2NiMnO6'),%doubleperovskite
       plotdata.Material(k).dpc001=3.871; %initial distance between two planes
       plotdata.Material(k).dpc111=3.871/sqrt(3)*2; %initial distance between two planes 
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
    
    elseif strcmp(plotdata.Material(k).Type,'La2NiO4'), %Ruddlesden-Popper phase
       plotdata.Material(k).d=12.42323; %initial distance between two planes  
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       
    elseif strcmp(plotdata.Material(k).Type,'(La_x,Sr_{1-x})CoO3'),
       plotdata.Material(k).dpc001=3.855; %initial distance between 2 planes - default x=0: SrCoO3 a=3.860; b=3.853; c=b;
       plotdata.Material(k).dpc111=3.855/sqrt(3); %initial distance between 2 planes - default x=0: SrTiO3
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xLSCOpanel,'Visible','on');
       
    elseif strcmp(plotdata.Material(k).Type,'LaTiO3'),
       plotdata.Material(k).dpc001=3.959; %initial distance between two planes https://materialsproject.org/materials/mp-8020/
       plotdata.Material(k).dpc111=3.959/sqrt(3); %initial distance between two planes https://materialsproject.org/materials/mp-8020/
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)'; 
 
    elseif strcmp(plotdata.Material(k).Type,'LaVO3'),
       plotdata.Material(k).dpc001=3.951; %initial distance between two planes https://materialsproject.org/materials/mp-19053/
       plotdata.Material(k).dpc111=3.951/sqrt(3); %initial distance between two planes https://materialsproject.org/materials/mp-19053/
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)'; 
       
    elseif strcmp(plotdata.Material(k).Type,'LSMO'),
       plotdata.Material(k).dpc001=3.873; %initial distance between two planes http://ematweb.cmi.ua.ac.be/emat/pdf/1214.pdf
       plotdata.Material(k).dpc111=3.873/sqrt(3); %initial distance between two planes http://ematweb.cmi.ua.ac.be/emat/pdf/1214.pdf
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)'; 
       
    elseif strcmp(plotdata.Material(k).Type,'Mg3N2'),
       plotdata.Material(k).dpc001=9.9528; %initial distance between two planes
       plotdata.Material(k).dpc111=9.9528/sqrt(3)/2; %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';  
    
    elseif strcmp(plotdata.Material(k).Type,'MgO'),
       plotdata.Material(k).dpc001=4.212; %initial distance between two planes
       plotdata.Material(k).dpc111=4.212/sqrt(3); %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';  
       
    elseif strcmp(plotdata.Material(k).Type,'MnTiO3'),
       plotdata.Material(k).dpc001=3.832; %initial distance between two planes https://materialsproject.org/materials/mp-19082/ V=112.583 for double unit cell
       plotdata.Material(k).dpc111=3.832/sqrt(3); %initial distance between two planes https://materialsproject.org/materials/mp-19082/ V=112.583 for double unit cell
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)'; 
 
     elseif strcmp(plotdata.Material(k).Type,'MoS2(P.63/m.m.c)'),
       plotdata.Material(k).dh0001=12.295; %initial distance between 2 planes
       plotdata.Material(k).d=plotdata.Material(k).dh0001;
       plotdata.Material(k).dh11bar20=3.16040/2;
       plotdata.Material(k).dh10bar10=3.16040*sqrt(3)/2;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosehOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='h(0001)';
 
     elseif strcmp(plotdata.Material(k).Type,'MoS2(P.-3.m.1)'),
       plotdata.Material(k).dh0001=5.94500; %initial distance between 2 planes
       plotdata.Material(k).d=plotdata.Material(k).dh0001;
       plotdata.Material(k).dh11bar20=3.19000/2;
       plotdata.Material(k).dh10bar10=3.19000*sqrt(3)/2;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosehOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='h(0001)';
       
     elseif strcmp(plotdata.Material(k).Type,'(Nd_{1-x},Ca_x)MnO3'),
       plotdata.Material(k).dpc001=3.834; %initial distance between 2 planes - default x=0: NdMnO3
       plotdata.Material(k).dpc111=3.834/sqrt(3); %initial distance between 2 planes - default x=0: NdMniO3
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xNCMpanel,'Visible','on');    
   
    elseif strcmp(plotdata.Material(k).Type,'NdCaMn2O6'),%doubleperovskite
       plotdata.Material(k).dpc001=(231.502/4)^(1/3); %initial distance between two planes https://materialsproject.org/materials/mp-1227144/;
       plotdata.Material(k).dpc111=(231.502/4)^(1/3)/sqrt(3)*2; %initial distance between two planes https://materialsproject.org/materials/mp-1227144/;
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       
    elseif strcmp(plotdata.Material(k).Type,'(Nd_x,La_{1-x})NiO3'),
       plotdata.Material(k).dpc001=3.857; %initial distance between 2 planes - default x=0: LaNiO3
       plotdata.Material(k).dpc111=3.857/sqrt(3); %initial distance between 2 planes - default x=0: LaNiO3
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xNLNpanel,'Visible','on');
      
   elseif strcmp(plotdata.Material(k).Type,'NdNiO2'),
       % tetragonal a=b=3.962, c=3.268 https://materialsproject.org/materials/mp-31063/
       plotdata.Material(k).d=3.268; %initial distance between 2 planes
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       
    elseif strcmp(plotdata.Material(k).Type,'NdNiO3'),
       plotdata.Material(k).dpc001=3.861; %initial distance between 2 planes - NdNiO3: pseudocubic a=3.861 cubicroot(230.143/4) https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
       plotdata.Material(k).dpc111=3.861/sqrt(3); %initial distance between 2 planes - NdNiO3: pseudocubic a=3.861 cubicroot(230.143/4) https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes       
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
 
    elseif strcmp(plotdata.Material(k).Type,'Nd_xNi_yO3'),%off-stoichiometry
       plotdata.Material(k).dpc001=3.861; %initial distance between 2 planes - NdNiO3: pseudocubic a=3.861 cubicroot(230.143/4) https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
       plotdata.Material(k).dpc111=3.861/sqrt(3); %initial distance between 2 planes - NdNiO3: pseudocubic a=3.861 cubicroot(230.143/4) https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes       
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xNNOpanel,'Visible','on');
       set(plotdata.Material(k).yNNOpanel,'Visible','on');
       
    elseif strcmp(plotdata.Material(k).Type,'Nd2NiMnO6'),%doubleperovskite
       plotdata.Material(k).dpc001=3.851; %initial distance between 2 planes https://aip.scitation.org/doi/full/10.1063/1.4906989?ver=pdfcov a=5.4162, b=5.4863, c=7.6718 V=a*b*c 4 unit cells 
       plotdata.Material(k).dpc111=3.851/sqrt(3)*2; %initial distance between 2 planes https://aip.scitation.org/doi/full/10.1063/1.4906989?ver=pdfcov a=5.4162, b=5.4863, c=7.6718 V=a*b*c 4 unit cells 
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes https://aip.scitation.org/doi/full/10.1063/1.4906989?ver=pdfcov a=5.4162, b=5.4863, c=7.6718 V=a*b*c 4 unit cells 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
            
   elseif strcmp(plotdata.Material(k).Type,'(Nd_{1-x},Sr_x)MnO3'),
       plotdata.Material(k).dpc001=3.834; %initial distance between 2 planes - default x=0: NdMnO3
       plotdata.Material(k).dpc111=3.834/sqrt(3); %initial distance between 2 planes - default x=0: NdMniO3
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xNSMpanel,'Visible','on');
   
    elseif strcmp(plotdata.Material(k).Type,'PbNiO3'),
       plotdata.Material(k).dpc001=3.811; %initial distance between 2 planes - https://materialsproject.org/materials/mp-974108/ V=55.330
       plotdata.Material(k).dpc111=3.811/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-974108/ V=55.330
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       
    elseif strcmp(plotdata.Material(k).Type,'(Pb_x,Sr_{1-x})TiO3'),
       plotdata.Material(k).dpc001=3.905; %initial distance between 2 planes - default x=0: SrTiO3
       plotdata.Material(k).dpc111=3.905/sqrt(3); %initial distance between 2 planes - default x=0: SrTiO3
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xPSTpanel,'Visible','on');

    elseif strcmp(plotdata.Material(k).Type,'PbTiO3'),
       plotdata.Material(k).dpc001=4.152; %initial distance between 2 planes - known bulk value for PbTiO3 with a=b=3.904
       plotdata.Material(k).dpc111=4.152/sqrt(3); %initial distance between 2 planes - known bulk value for PbTiO3 with a=b=3.904
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).TBPolarization,'Visible','on');
       
    elseif strcmp(plotdata.Material(k).Type,'Pb(Zr_x,Ti_{1-x})O3'),
       plotdata.Material(k).dpc001=4.152; %initial distance between 2 planes - default x=0: PbTiO3
       plotdata.Material(k).dpc111=4.152/sqrt(3); %initial distance between 2 planes - default x=0: PbTiO3
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xPZTpanel,'Visible','on');

    elseif strcmp(plotdata.Material(k).Type,{'PrBa2Cu3O7'}),
       plotdata.Material(k).d=11.916; %initial distance between 2 planes
       set(plotdata.Material(k).OrientationPanel,'visible','off');   
   
   elseif strcmp(plotdata.Material(k).Type,'PrNiO2'),
       % using the lattice parameters of SrCuO2 tetragonal a=b=3.948, c=3.485 https://materialsproject.org/materials/mp-37514/
       plotdata.Material(k).d=3.485; %initial distance between 2 planes
       set(plotdata.Material(k).OrientationPanel,'visible','off');
        
    elseif strcmp(plotdata.Material(k).Type,'PrNiO3'),
       plotdata.Material(k).dpc001=3.872; %initial distance between 2 planes - https://materialsproject.org/materials/mp-22280/ V=232.262 4 unit cells
       plotdata.Material(k).dpc111=3.872/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-22280/ V=232.262 4 unit cells
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';   
       
    elseif strcmp(plotdata.Material(k).Type,'PrVO3'),
       plotdata.Material(k).dpc001=3.936; %initial distance between 2 planes - https://materialsproject.org/materials/mp-1069346/
       plotdata.Material(k).dpc111=3.936/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-1069346/
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';   
       
    elseif strcmp(plotdata.Material(k).Type,'SmNiO3'),
       plotdata.Material(k).dpc001=3.794; %initial distance between 2 planes - https://materialsproject.org/materials/mp-1099668/
       plotdata.Material(k).dpc111=3.794/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-1099668/
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';  
       
   elseif strcmp(plotdata.Material(k).Type,{'Sr3Al2O6'}),
       plotdata.Material(k).d=7.999957; %initial distance between 2 planes https://materialsproject.org/materials/mp-3393/  
       set(plotdata.Material(k).OrientationPanel,'visible','off');
    
   elseif strcmp(plotdata.Material(k).Type,{'SrCoO2.5'}),
       plotdata.Material(k).d=7.8194; %initial distance between 2 planes    
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       
   elseif strcmp(plotdata.Material(k).Type,'SrCoO3'),
       %https://materialsproject.org/materials/mp-505766/ alpha=beta=gamma=90
       %a=3.860; b=3.853; c=b; d=cubicroot(a*b*c)
       plotdata.Material(k).dpc001=3.855; %initial distance between 2 planes
       plotdata.Material(k).dpc111=3.855/sqrt(3); %initial distance between 2 planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)'; 
   
   elseif strcmp(plotdata.Material(k).Type,'SrCrO3'),
       plotdata.Material(k).dpc001=3.8185; %initial distance between 2 planes
       plotdata.Material(k).dpc111=3.8185/sqrt(3); %initial distance between 2 planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)'; 
       
   elseif strcmp(plotdata.Material(k).Type,'SrCuO2'),
       % tetragonal a=b=3.948, c=3.485 https://materialsproject.org/materials/mp-37514/
       plotdata.Material(k).d=3.485; %initial distance between 2 planes
       set(plotdata.Material(k).OrientationPanel,'visible','off');
  
   elseif strcmp(plotdata.Material(k).Type,'SrIrO3'),
       plotdata.Material(k).dpc001=3.998; %initial distance between 2 planes - https://materialsproject.org/materials/mp-1016848/
       plotdata.Material(k).dpc111=3.998/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-1016848/
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)'; 
       
    elseif strcmp(plotdata.Material(k).Type,'SrMoO3'),
       plotdata.Material(k).dpc001=4.082; %initial distance between 2 planes - https://materialsproject.org/materials/mp-18747/
       plotdata.Material(k).dpc111=4.082/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-18747/
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)'; 
       
    elseif strcmp(plotdata.Material(k).Type,'SrRuO3'),
       plotdata.Material(k).dpc001=3.985; %initial distance between 2 planes - https://materialsproject.org/materials/mp-4346/
       plotdata.Material(k).dpc111=plotdata.Material(k).dpc001/sqrt(3); 
       plotdata.Material(k).d=plotdata.Material(k).dpc001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)'; 
       
    elseif strcmp(plotdata.Material(k).Type,'SrTiO3'),
       plotdata.Material(k).dpc001=3.905; %initial distance between 2 planes
       plotdata.Material(k).dpc111=plotdata.Material(k).dpc001/sqrt(3); 
       plotdata.Material(k).d=plotdata.Material(k).dpc001;      
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)'; 
       
     elseif strcmp(plotdata.Material(k).Type,'SrVO3'),
       plotdata.Material(k).dpc001=3.901; %initial distance between 2 planes - cubic - https://materialsproject.org/materials/mp-18717/
       plotdata.Material(k).dpc111=plotdata.Material(k).dpc001/sqrt(3); 
       plotdata.Material(k).d=plotdata.Material(k).dpc001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';  
       
     elseif strcmp(plotdata.Material(k).Type,'Tm3Fe5O12'),
       plotdata.Material(k).dpc001=12.2325/4; %initial distance between 2 planes - known bulk value for Tm3Fe5O12 with a=b=c=12.2325 (Landolt/Bornstein)
       plotdata.Material(k).dpc111=12.2325/4/sqrt(3);
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
  
     elseif strcmp(plotdata.Material(k).Type,'WS2'),
       plotdata.Material(k).dh0001=12.32300; %initial distance between 2 planes
       plotdata.Material(k).d=plotdata.Material(k).dh0001;
       %set(plotdata.Material(k).OrientationPanel,'visible','on');
       %set(plotdata.Material(k).choosehOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='h(0001)';
       
     elseif strcmp(plotdata.Material(k).Type,{'YBa2Cu3O7'}),
       plotdata.Material(k).d=11.824; %initial distance between 2 planes - tetragonal - https://materialsproject.org/materials/mp-20674/
       set(plotdata.Material(k).OrientationPanel,'visible','off');
    
     elseif strcmp(plotdata.Material(k).Type,'YBiO3'),
       %perovskite
       %https://materialsproject.org/materials/mvc-13598/
       plotdata.Material(k).dpc001=4.408542; %initial distance between two planes
       plotdata.Material(k).dpc111=4.408542/sqrt(3); %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       
   elseif strcmp(plotdata.Material(k).Type,'Y3Fe5O12'),
       plotdata.Material(k).dpc001=12.376/4; %initial distance between 2 planes - known bulk value for Y3Fe5O12 with a=b=c=12.376 (Landolt/Bornstein)
       plotdata.Material(k).dpc111=12.376/4/sqrt(3);
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
   
   elseif strcmp(plotdata.Material(k).Type,'(Y_xTm_{3-x})Fe5O12'),
       plotdata.Material(k).dpc001=12.2325/4; %initial distance between two planes (default x=0 i.e. Tm_3Fe5O12)
       plotdata.Material(k).dpc111=12.2325/4/sqrt(3); %initial distance between two planes (default x=0 i.e. Tm_3Fe5O12)
       plotdata.Material(k).d=plotdata.Material(k).dpc001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xYTmIGpanel,'Visible','on');
       
   elseif strcmp(plotdata.Material(k).Type,'YNiO3'),
       %https://materialsproject.org/materials/mvc-15448/
       %a=3.753 b=c=3.759 alpha=89.967 beta=gamma=90 V=53.025
       %d=cubicroot(53.025)=3.757
       plotdata.Material(k).d=3.757; %initial distance between 2 planes
       plotdata.Material(k).dpc111=plotdata.Material(k).dpc001/sqrt(3); 
       plotdata.Material(k).d=plotdata.Material(k).dpc001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';    
       
   elseif strcmp(plotdata.Material(k).Type,'(Zn_x,Mg_{1-x})O'),
       % same structure as ZnO but with Mg substituting Zn for Mg<30%
       % ZnO: hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
       % substituting with Mg increases out-of-plane lattice parameter
       plotdata.Material(k).dh0001=5.20380; %initial distance between 2 planes
       plotdata.Material(k).dh11bar20=3.2494/2;
       plotdata.Material(k).dh10bar10=3.2494*sqrt(3)/2;
       plotdata.Material(k).d=plotdata.Material(k).dh0001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosehOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='h(0001)';
       set(plotdata.Material(k).xZMOpanel,'Visible','on');
            
   elseif strcmp(plotdata.Material(k).Type,'Zn3N2'),
       % Ia-3 a=b=c=9.76910 alpha=beta=gamma=90
       plotdata.Material(k).dpc001=9.76910; %initial distance between two planes
       plotdata.Material(k).dpc111=plotdata.Material(k).dpc001/sqrt(3)/2; %initial distance between two planes
       plotdata.Material(k).d=plotdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='pc(001)';
       
   elseif strcmp(plotdata.Material(k).Type,'ZnO'),
       plotdata.Material(k).dh0001=5.20380; %initial distance between 2 planes
       plotdata.Material(k).dh11bar20=3.2494/2;
       plotdata.Material(k).dh10bar10=3.2494*sqrt(3)/2;
       plotdata.Material(k).d=plotdata.Material(k).dh0001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosehOrientation,'visible','on','value',1);
       plotdata.Material(k).Orientation='h(0001)';
   end;
   plotdata.Material(k).d=round(plotdata.Material(k).d,4);
   set(plotdata.Material(k).dWrite,'String',plotdata.Material(k).d);
   set(plotdata.Material(k).dSlide,'Value',plotdata.Material(k).d);
   ProgcNSimu
   
   %%   
case {'choosehOrientation(1)','choosehOrientation(2)','choosehOrientation(3)','choosehOrientation(4)','choosehOrientation(5)','choosehOrientation(6)'}
   k=str2num(entry(20));
   str=get(plotdata.Material(k).choosehOrientation,'String');
   entry=get(plotdata.Material(k).choosehOrientation,'Value');
   MaterialhOrientation=strsplit(str(entry,:));
   plotdata.Material(k).Orientation=MaterialhOrientation(1);
   if strcmp(plotdata.Material(k).Orientation,'h(0001)'),
              plotdata.Material(k).d=round(plotdata.Material(k).dh0001,4);
   elseif strcmp(plotdata.Material(k).Orientation,'h(11-20)'),
              plotdata.Material(k).d=round(plotdata.Material(k).dh11bar20,4);
   elseif strcmp(plotdata.Material(k).Orientation,'h(10-10)'),
              plotdata.Material(k).d=round(plotdata.Material(k).dh10bar10,4);
   else
        msgbox('Error with orientation');
   end;
   set(plotdata.Material(k).dWrite,'String',plotdata.Material(k).d);
   set(plotdata.Material(k).dSlide,'Value',plotdata.Material(k).d);
   ProgcNSimu
%%
    case {'choosemOrientation(1)','choosemOrientation(2)','choosemOrientation(3)','choosemOrientation(4)','choosemOrientation(5)','choosemOrientation(6)'}
   k=str2num(entry(20));
   str=get(plotdata.Material(k).choosemOrientation,'String');
   entry=get(plotdata.Material(k).choosemOrientation,'Value');
   MaterialmOrientation=strsplit(str(entry,:));
   plotdata.Material(k).Orientation=MaterialmOrientation(1);
%%
    case {'choosepcOrientation(1)','choosepcOrientation(2)','choosepcOrientation(3)','choosepcOrientation(4)','choosepcOrientation(5)','choosepcOrientation(6)'}
   k=str2num(entry(21));
   str=get(plotdata.Material(k).choosepcOrientation,'String');
   entry=get(plotdata.Material(k).choosepcOrientation,'Value');
   MaterialpcOrientation=strsplit(str(entry,:));
   plotdata.Material(k).Orientation=MaterialpcOrientation(1);
   if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
       plotdata.Material(k).d=round(plotdata.Material(k).dpc001,4);
   elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
       plotdata.Material(k).d=round(plotdata.Material(k).dpc111,4);
   else
        msgbox('Error with orientation');
   end;
   set(plotdata.Material(k).dWrite,'String',plotdata.Material(k).d);
   set(plotdata.Material(k).dSlide,'Value',plotdata.Material(k).d);
   ProgcNSimu
%%

case {'chooseddistrib(1)','chooseddistrib(2)','chooseddistrib(3)','chooseddistrib(4)','chooseddistrib(5)','chooseddistrib(6)'}
   k=str2num(entry(16));
   
   if plotdata.Material(k).ddistrib.Value==1,
    plotdata.Material(k).dDistribution='constant';
    set(plotdata.Material(k).dconst,'visible','on');
    set(plotdata.Material(k).dexpdistrib,'visible','off');
    set(plotdata.Material(k).dJacobidistrib,'visible','off');
    set(plotdata.Material(k).dupload,'visible','off');
   elseif plotdata.Material(k).ddistrib.Value==2,
    plotdata.Material(k).dDistribution='exp';
    set(plotdata.Material(k).dconst,'visible','off');
    set(plotdata.Material(k).dexpdistrib,'visible','on');
    set(plotdata.Material(k).dJacobidistrib,'visible','off');    
    set(plotdata.Material(k).dupload,'visible','off');
   elseif plotdata.Material(k).ddistrib.Value==3,
    plotdata.Material(k).dDistribution='Jacobi';
    set(plotdata.Material(k).dconst,'visible','off');
    set(plotdata.Material(k).dexpdistrib,'visible','off');
    set(plotdata.Material(k).dJacobidistrib,'visible','on');
    set(plotdata.Material(k).dupload,'visible','off');
   elseif plotdata.Material(k).ddistrib.Value==4,
    plotdata.Material(k).dDistribution='upload';
    set(plotdata.Material(k).dconst,'visible','off');
    set(plotdata.Material(k).dexpdistrib,'visible','off');
    set(plotdata.Material(k).dJacobidistrib,'visible','off');
    set(plotdata.Material(k).dupload,'visible','on');
else
    warning('Error in toggle function for dDistribution')
end
ProgcNSimu;
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
      plotdata.LMax=2*plotdata.Substrate.d/plotdata.lambda;
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
   
end
%===A utility to center the window on the screen============
function pos = centerfig(width,height)

% Find the screen size in pixels
screen_s = get(0,'ScreenSize');
pos = [screen_s(3)/2 - width/2, screen_s(4)/2 - height/2, width, height];

function toggle(source,eventdata)
global plotdata;
if strcmp(get(get(source,'SelectedObject'),'String'),'L')
    plotdata.FitDisplay.haxis='L';
    plotdata.LMin=2*plotdata.Substrate.d*sin(plotdata.TThetaMin/2*pi/180)/plotdata.lambda;
    plotdata.LMax=2*plotdata.Substrate.d*sin(plotdata.TThetaMax/2*pi/180)/plotdata.lambda;
    set(plotdata.chooseXMin,'String',plotdata.LMin);
    set(plotdata.chooseXMax,'String',plotdata.LMax);
    
elseif strcmp(get(get(source,'SelectedObject'),'String'),'2Theta')
    plotdata.FitDisplay.haxis='2Theta';   
    set(plotdata.chooseXMin,'String',plotdata.TThetaMin);
    set(plotdata.chooseXMax,'String',plotdata.TThetaMax);

else
    warning('Error in Toggle function')
end
set(plotdata.chooseXMinText,'String',[plotdata.FitDisplay.haxis,' min:']);
set(plotdata.chooseXMaxText,'String',[plotdata.FitDisplay.haxis,' max:']);
PlotFitAndData
%%
function toggleSubstrateOrientation(source,eventdata)
global plotdata;
plotdata.Substrate.Orientation=get(get(source,'SelectedObject'),'String');
if and(strcmp(char(plotdata.Substrate.Type),'Al2O3'),strcmp(char(plotdata.Substrate.Orientation),'(11-20)h')),
    SubstrateWarning(plotdata.Substrate.Orientation,plotdata.Substrate.Type);
    plotdata.Substrate.Orientation='(0001)h';
    plotdata.Substrate.TBh0001.Value=true;
elseif and(strcmp(char(plotdata.Substrate.Type),'Al2O3'),strcmp(char(plotdata.Substrate.Orientation),'(10-10)h')),
    SubstrateWarning(plotdata.Substrate.Orientation,plotdata.Substrate.Type);
    plotdata.Substrate.Orientation='(0001)h';
    plotdata.Substrate.TBh0001.Value=true;
elseif and(strcmp(char(plotdata.Substrate.Type),'NdGaO3'),strcmp(char(plotdata.Substrate.Orientation),'(001)o')),
    SubstrateWarning(plotdata.Substrate.Orientation,plotdata.Substrate.Type);
    plotdata.Substrate.Orientation='(110)o';
    plotdata.Substrate.TBo110.Value=true;  
elseif and(strcmp(char(plotdata.Substrate.Type),'TbScO3'),strcmp(char(plotdata.Substrate.Orientation),'(110)o')),
    SubstrateWarning(plotdata.Substrate.Orientation,plotdata.Substrate.Type);
    plotdata.Substrate.Orientation='(001)o';
    plotdata.Substrate.TBo001.Value=true;
elseif and(strcmp(char(plotdata.Substrate.Type),'TbScO3'),strcmp(char(plotdata.Substrate.Orientation),'(101)o')),
    SubstrateWarning(plotdata.Substrate.Orientation,plotdata.Substrate.Type);
    plotdata.Substrate.Orientation='(001)o';
    plotdata.Substrate.TBo001.Value=true; 
end
Substrate;
ProgcNSimu;
%%
function SubstrateWarning(Orientation,Name)
    dlg = dialog('Position',[300 300 250 300],'Name','Warning');

    txt = uicontrol('Parent',dlg,...
               'Style','text',...
               'Position',[20 80 210 200],...
               'String',['The ' Orientation ' orientation for ' Name ' is not available yet. Please contact Celine.Lichtensteiger@unige.ch if needed.']);

    btn = uicontrol('Parent',dlg,...
               'Position',[85 20 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');
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