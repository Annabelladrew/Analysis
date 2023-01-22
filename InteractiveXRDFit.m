% Created by Celine Lichtensteiger
% Journal of Applied Crystallography - Computer programs - 51(6)p.1745 (2018)
% https://doi.org/10.1107/S1600576718012840
% This is the heart of the program
% It initializes all the variables and creates the Graphical User Interface
%*********************************

function InteractiveXRDFit(entry)

global fitdata plotdata;
00
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
   fitdata.Substrate.d=dsubstrateSrTiO3c001;                              
   plotdata.Substrate.g=gsubstrateSrTiO3c001;
   plotdata.Substrate.thickness=2e4*3.905; %!!!all substrates have to be recalculated if this value is changed!!! 10e4*3.905;%2e4*3.905;%reference thickness in [A] for substrate calculations
   fitdata.Substrate.Type='SrTiO3';                       
   % Bottom Layer
   fitdata.Material(1).Aexp=0;
   fitdata.Material(1).Bexp=1;
   fitdata.Material(1).Cexp=4;
   fitdata.Material(1).d=4;
   fitdata.Material(1).dBorder=2;
   fitdata.Material(1).dCenter=4;
   fitdata.Material(1).dDistribution='constant';
   %fitdata.Material(1).d=2.309;%4/sqrt(3);
   fitdata.Material(1).m=1; %modulus of the Jacobi Elliptic function
   fitdata.Material(1).Meandexp=4;
   fitdata.Material(1).N=0;
   fitdata.Material(1).Orientation='none';
   fitdata.Material(1).Polarization=0;
   fitdata.Material(1).Type='none';
   fitdata.Material(1).xBST=0;
   fitdata.Material(1).xNLN=0;
   fitdata.Material(1).xPST=0;
   fitdata.Material(1).xPZT=0;
   fitdata.Material(1).xYTmIG=0;
   fitdata.Material(1).xZMO=1;
   plotdata.Material(1).V=4^2; %unit cell volume in [A]
   % Min-Max:
   plotdata.Material(1).AexpMax=10;      
   plotdata.Material(1).AexpMin=-10;     
   plotdata.Material(1).BexpMax=1000;    
   plotdata.Material(1).BexpMin=-1000;   
   plotdata.Material(1).MeandexpMax=14;  
   plotdata.Material(1).MeandexpMin=0;   
   plotdata.Material(1).dMax=14;         
   plotdata.Material(1).dMin=0;          
   plotdata.Material(1).NMax=1000;       
   plotdata.Material(1).NMin=0;   
   plotdata.Material(1).xBSTMax=1;
   plotdata.Material(1).xBSTMin=0;
   plotdata.Material(1).xNLNMax=1;
   plotdata.Material(1).xNLNMin=0;
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
       fitdata.Material(k)=fitdata.Material(1);
       plotdata.Material(k)=plotdata.Material(1);
   end;
   % dDisplay
   plotdata.dDisplay.Width=660;
   plotdata.dDisplay.Height=400;
   plotdata.dDisplay.centerfig=centerfig(plotdata.dDisplay.Width,plotdata.dDisplay.Height)+[400,-250,0,0];
   % RMS
   plotdata.RMS.hVerticalLines=[];
   plotdata.RMS.LMax=[];
   plotdata.RMS.LMin=[];
   plotdata.RMS.TThetaMax=[];
   plotdata.RMS.TThetaMin=[];
   % FitDisplay
   fitdata.FitDisplay.haxis='2Theta';
   plotdata.FitDisplay.Width=660;
   plotdata.FitDisplay.Height=400;
   plotdata.FitDisplay.centerfig=centerfig(plotdata.FitDisplay.Width,plotdata.FitDisplay.Height)-[400,0,0,0];
   plotdata.LMin=0;
   plotdata.LMax=2*fitdata.Substrate.d/plotdata.lambda; 
   fitdata.ScalingIntensity=1;
   plotdata.TThetaMin=0;
   plotdata.TThetaMax=180;
   % Main Panel
   plotdata.MainPanel.Width=792;    
   plotdata.MainPanel.Height=420;   
   plotdata.MainPanel.centerfig=centerfig(plotdata.MainPanel.Width,plotdata.MainPanel.Height)+[400,+250,0,0];
   fitdata.Repetition.N=1;
   plotdata.Repetition.Max=50;
   plotdata.Repetition.Min=0;
   % Misc
   plotdata.fit.compare.x=[];
   plotdata.fit.compare.y=[];
   plotdata.fit.x.TTheta=[0:0.01:180]';
   plotdata.fit.x.L=2*fitdata.Substrate.d*sin(plotdata.fit.x.TTheta/2*pi/180)/plotdata.lambda;
   plotdata.fit.y=10^10+zeros(size([0:0.01:180]'));
   plotdata.mMax=1;
   plotdata.mMin=0;
   plotdata.measure.compare.y=[];
   plotdata.measure.x=[];
   plotdata.measure.y=[];
   plotdata.measure.filename='none';
   plotdata.mu=1.5e4; %Penetration depth
   plotdata.prog = mfilename;
   plotdata.Q=2*pi*(2*sin([0:0.01:180]'/2*pi/180)/1.5406)'; %Momentum transfert
   plotdata.QQ=plotdata.Q.^2;
   plotdata.TotalThickness=plotdata.Substrate.thickness+100000; %in [A]
   plotdata.d=[];
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
      'String','none|DyScO3|Gd3Ga5O12|GdScO3|KTaO3|LaAlO3|LaGaO3|LSAT|MgO|Nb:SrTiO30.5%wt|NdAlO3|NdGaO3|PMN-PT|Si|SrLaAlO4|SrLaGaO4|SrTiO3|TbScO3|TiO2|YAlO3|YSZ|ZnO',...
      'value',17,...
      'Units','normalized',...
      'Position',[0,0.47,1,0.5],...
      'visible','on',...
      'CallBack',[plotdata.prog,' chooseSubstrate'] ); 

%%DyScO3
plotdata.Substrate.DyScO3.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','DyScO3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.85]);

uicontrol(plotdata.Substrate.DyScO3.TBorientation,'style','text',...
       'string',['Orthorhombic a=',num2str(round(DyScO3o110.a,4)),'A b=',num2str(round(DyScO3o110.b,4)),'A c=',num2str(round(DyScO3o110.c,4)),'A'],...% Meley, APL Materials, 6, 046102 (2018) 1/d=sqrt(h^2/a^2+k^2/b^2+l^2/v^2)',...
       'unit','normalized','position',[0,0.8,1,0.2]);

% Create three buttons in the button group.
plotdata.Substrate.DyScO3.TBortho001 = uicontrol('Style','radiobutton','String','(001)o',...
    'Units','normalized','pos',[0,0.70,1,0.15],'parent',plotdata.Substrate.DyScO3.TBorientation);

plotdata.Substrate.DyScO3.TBpc001 = uicontrol('Style','radiobutton','String','(110)o',...
    'Units','normalized','pos',[0,0.55,1,0.15],'parent',plotdata.Substrate.DyScO3.TBorientation);
uicontrol(plotdata.Substrate.DyScO3.TBorientation,'style','text','string',['equivalent to (001)pc, d=',num2str(round(dsubstrateDyScO3pc001,4)),'A'],'Units','normalized','pos',[0,0.35,0.8,0.2]);

plotdata.Substrate.DyScO3.TBpc111 = uicontrol('Style','radiobutton','String','(101)o',...
    'Units','normalized','pos',[0,0.2,1,0.15],'parent',plotdata.Substrate.DyScO3.TBorientation);
uicontrol(plotdata.Substrate.DyScO3.TBorientation,'style','text','string',['equivalent to (111)pc, d=',num2str(round(dsubstrateDyScO3pc111,4)),'A'],'Units','normalized','pos',[0,0,0.8,0.2]);

set(plotdata.Substrate.DyScO3.TBorientation,'SelectionChangeFcn',@toggleSubstrateDyScO3orientation);

%%Gd3Ga5O12
plotdata.Substrate.Gd3Ga5O12.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','Gd3Ga5O12 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create two buttons in the button group.
plotdata.Substrate.Gd3Ga5O12.TBc001 = uicontrol('Style','radiobutton','String','cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.Gd3Ga5O12.TBorientation);
uicontrol(plotdata.Substrate.Gd3Ga5O12.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateGd3Ga5O12c001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.Gd3Ga5O12.TBc111 = uicontrol('Style','radiobutton','String','cubic (111)',...
    'Units','normalized','pos',[0,0.15,1,0.2],'parent',plotdata.Substrate.Gd3Ga5O12.TBorientation);
uicontrol(plotdata.Substrate.Gd3Ga5O12.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateGd3Ga5O12c111,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.Gd3Ga5O12.TBorientation,'SelectionChangeFcn',@toggleSubstrateGd3Ga5O12orientation);

%%GdScO3
plotdata.Substrate.GdScO3.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','GdScO3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create two buttons in the button group.
plotdata.Substrate.GdScO3.TBpc001 = uicontrol('Style','radiobutton','String','pseudo-cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.GdScO3.TBorientation);
uicontrol(plotdata.Substrate.GdScO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateGdScO3pc001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.GdScO3.TBpc111 = uicontrol('Style','radiobutton','String','pseudo-cubic (111)',...
    'Units','normalized','pos',[0,0.15,1,0.2],'parent',plotdata.Substrate.GdScO3.TBorientation);
uicontrol(plotdata.Substrate.GdScO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateGdScO3pc111,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.GdScO3.TBorientation,'SelectionChangeFcn',@toggleSubstrateGdScO3orientation);

%%KTaO3
plotdata.Substrate.KTaO3.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','KTaO3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create two buttons in the button group.
plotdata.Substrate.KTaO3.TBc001 = uicontrol('Style','radiobutton','String','cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.KTaO3.TBorientation);
uicontrol(plotdata.Substrate.KTaO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateKTaO3c001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.KTaO3.TBc111 = uicontrol('Style','radiobutton','String','cubic (111)',...
    'Units','normalized','pos',[0,0.15,1,0.2],'parent',plotdata.Substrate.KTaO3.TBorientation);
uicontrol(plotdata.Substrate.KTaO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateKTaO3c111,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.KTaO3.TBorientation,'SelectionChangeFcn',@toggleSubstrateKTaO3orientation);

%%LaAlO3
plotdata.Substrate.LaAlO3.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','LaAlO3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create two buttons in the button group.
plotdata.Substrate.LaAlO3.TBpc001 = uicontrol('Style','radiobutton','String','pseudo-cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.LaAlO3.TBorientation);
uicontrol(plotdata.Substrate.LaAlO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateLaAlO3pc001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.LaAlO3.TBpc111 = uicontrol('Style','radiobutton','String','pseudo-cubic (111)',...
    'Units','normalized','pos',[0,0.15,1,0.2],'parent',plotdata.Substrate.LaAlO3.TBorientation);
uicontrol(plotdata.Substrate.LaAlO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateLaAlO3pc111,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.LaAlO3.TBorientation,'SelectionChangeFcn',@toggleSubstrateLaAlO3orientation);

%%LaGaO3
plotdata.Substrate.LaGaO3.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','LaGaO3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.85]);

uicontrol(plotdata.Substrate.LaGaO3.TBorientation,'style','text',...
       'string',['Orthorhombic a=',num2str(LaGaO3o110.a),'A b=',num2str(round(LaGaO3o110.b,4)),'A d=',num2str(round(LaGaO3o110.c,4)),'A'],...% (Mathews, PRB80,064408 2009) 1/d=sqrt(h^2/a^2+k^2/b^2+l^2/v^2)',...
       'unit','normalized','position',[0,0.8,1,0.2]); 
% Create three buttons in the button group.
plotdata.Substrate.LaGaO3.TBo110 = uicontrol('Style','radiobutton','String','(110)o',...
    'Units','normalized','pos',[0,0.6,1,0.15],'parent',plotdata.Substrate.LaGaO3.TBorientation);
uicontrol(plotdata.Substrate.LaGaO3.TBorientation,'style','text','string',['equivalent to (001)pc, d=',num2str(round(dsubstrateLaGaO3o110,4)),'A'],'Units','normalized','pos',[0.2,0.5,0.8,0.15]);

plotdata.Substrate.LaGaO3.TBo101 = uicontrol('Style','radiobutton','String','(101)o',...
    'Units','normalized','pos',[0,0.4,1,0.15],'parent',plotdata.Substrate.LaGaO3.TBorientation);
uicontrol(plotdata.Substrate.LaGaO3.TBorientation,'style','text','string',['equivalent to (111)pc, d=',num2str(round(dsubstrateLaGaO3o101,4)),'A'],'Units','normalized','pos',[0.2,0.3,0.8,0.15]);

plotdata.Substrate.LaGaO3.TBo001 = uicontrol('Style','radiobutton','String','(001)o',...
    'Units','normalized','pos',[0,0.2,1,0.15],'parent',plotdata.Substrate.LaGaO3.TBorientation);
uicontrol(plotdata.Substrate.LaGaO3.TBorientation,'style','text','string',['equivalent to (100)pc, d=',num2str(round(dsubstrateLaGaO3o001,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.15]);

set(plotdata.Substrate.LaGaO3.TBorientation,'SelectionChangeFcn',@toggleSubstrateLaGaO3orientation);

%%LSAT
plotdata.Substrate.LSAT.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','LSAT :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create two buttons in the button group.
plotdata.Substrate.LSAT.TBpc001 = uicontrol('Style','radiobutton','String','pseudo-cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.LSAT.TBorientation);
uicontrol(plotdata.Substrate.LSAT.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateLSATpc001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.LSAT.TBpc111 = uicontrol('Style','radiobutton','String','pseudo-cubic (111)',...
    'Units','normalized','pos',[0,0.15,1,0.2],'parent',plotdata.Substrate.LSAT.TBorientation);
uicontrol(plotdata.Substrate.LSAT.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateLSATpc111,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.LSAT.TBorientation,'SelectionChangeFcn',@toggleSubstrateLSATorientation);

%%MgO
plotdata.Substrate.MgO.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','MgO :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create two buttons in the button group.
plotdata.Substrate.MgO.TBc001 = uicontrol('Style','radiobutton','String','cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.MgO.TBorientation);
uicontrol(plotdata.Substrate.MgO.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateMgOc001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.MgO.TBc111 = uicontrol('Style','radiobutton','String','cubic (111)',...
    'Units','normalized','pos',[0,0.15,1,0.2],'parent',plotdata.Substrate.MgO.TBorientation);
uicontrol(plotdata.Substrate.MgO.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateMgOc111,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.MgO.TBorientation,'SelectionChangeFcn',@toggleSubstrateMgOorientation);

%%Nb:SrTiO30.5%wt
plotdata.Substrate.Nb05STO.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','on',...
    'Title','Nb:SrTiO3 0.5%wt :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create two buttons in the button group.
plotdata.Substrate.Nb05STO.TBc001 = uicontrol('Style','radiobutton','String','cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.Nb05STO.TBorientation);
uicontrol(plotdata.Substrate.Nb05STO.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateNb05STOc001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.Nb05STO.TBc111 = uicontrol('Style','radiobutton','String','cubic (111)',...
    'Units','normalized','pos',[0,0.15,1,0.2],'parent',plotdata.Substrate.Nb05STO.TBorientation);
uicontrol(plotdata.Substrate.Nb05STO.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateNb05STOc111,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.Nb05STO.TBorientation,'SelectionChangeFcn',@toggleSubstrateNb05STOorientation);

%%NdAlO3
plotdata.Substrate.NdAlO3.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','NdAlO3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create two buttons in the button group.
plotdata.Substrate.NdAlO3.TBpc001 = uicontrol('Style','radiobutton','String','pseudo-cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.NdAlO3.TBorientation);
uicontrol(plotdata.Substrate.NdAlO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateNdAlO3pc001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.NdAlO3.TBpc111 = uicontrol('Style','radiobutton','String','pseudo-cubic (111)',...
    'Units','normalized','pos',[0,0.15,1,0.2],'parent',plotdata.Substrate.NdAlO3.TBorientation);
uicontrol(plotdata.Substrate.NdAlO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateNdAlO3pc111,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.NdAlO3.TBorientation,'SelectionChangeFcn',@toggleSubstrateNdAlO3orientation);

%%NdGaO3
plotdata.Substrate.NdGaO3.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','NdGaO3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.85]);

uicontrol(plotdata.Substrate.NdGaO3.TBorientation,'style','text',...
       'string',['Orthorhombic a=',num2str(NdGaO3o110.a),'A b=',num2str(round(NdGaO3o110.b,4)),'A d=',num2str(round(NdGaO3o110.c,4)),'A'],...% (Mathews, PRB80,064408 2009) 1/d=sqrt(h^2/a^2+k^2/b^2+l^2/v^2)',...
       'unit','normalized','position',[0,0.7,1,0.2]);
% Create two buttons in the button group.
plotdata.Substrate.NdGaO3.TBpc001 = uicontrol('Style','radiobutton','String','(110)o',...
    'Units','normalized','pos',[0,0.55,1,0.15],'parent',plotdata.Substrate.NdGaO3.TBorientation);
uicontrol(plotdata.Substrate.NdGaO3.TBorientation,'style','text','string',['equivalent to (001)pc, d=',num2str(round(dsubstrateNdGaO3pc001,4)),'A'],'Units','normalized','pos',[0,0.35,0.8,0.2]);

plotdata.Substrate.NdGaO3.TBpc111 = uicontrol('Style','radiobutton','String','(101)o',...
    'Units','normalized','pos',[0,0.2,1,0.15],'parent',plotdata.Substrate.NdGaO3.TBorientation);
uicontrol(plotdata.Substrate.NdGaO3.TBorientation,'style','text','string',['equivalent to (111)pc, d=',num2str(round(dsubstrateNdGaO3pc111,4)),'A'],'Units','normalized','pos',[0,0,0.8,0.2]);

set(plotdata.Substrate.NdGaO3.TBorientation,'SelectionChangeFcn',@toggleSubstrateNdGaO3orientation);

%%PMN71PT29t001
plotdata.Substrate.PMN71PT29.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','PMN71PT29 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create one button in the button group.
plotdata.Substrate.PMN71PT29.TBt001 = uicontrol('Style','radiobutton','String','tetragonal (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.PMN71PT29.TBorientation);
uicontrol(plotdata.Substrate.PMN71PT29.TBorientation,'style','text','string',['d=',num2str(round(dsubstratePMN71PT29t001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);

%%Si
plotdata.Substrate.Si.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','Si :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create one button in the button group.
plotdata.Substrate.Si.TBc001 = uicontrol('Style','radiobutton','String','tetragonal (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.Si.TBorientation);
uicontrol(plotdata.Substrate.Si.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateSic001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);

%%SrLaAlO4
plotdata.Substrate.SrLaAlO4.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','SrLaAlO4 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create one button in the button group.
plotdata.Substrate.SrLaAlO4.TBt001 = uicontrol('Style','radiobutton','String','tetragonal (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.SrLaAlO4.TBorientation);
uicontrol(plotdata.Substrate.SrLaAlO4.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateSrLaAlO4t001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);

%%SrLaGaO4
plotdata.Substrate.SrLaGaO4.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','SrLaGaO4 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create one button in the button group.
plotdata.Substrate.SrLaGaO4.TBt001 = uicontrol('Style','radiobutton','String','tetragonal (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.SrLaGaO4.TBorientation);
uicontrol(plotdata.Substrate.SrLaGaO4.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateSrLaGaO4t001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);

%%SrTiO3
plotdata.Substrate.SrTiO3.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','on',...
    'Title','SrTiO3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create two buttons in the button group.
plotdata.Substrate.SrTiO3.TBc001 = uicontrol('Style','radiobutton','String','cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.SrTiO3.TBorientation);
uicontrol(plotdata.Substrate.SrTiO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateSrTiO3c001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.SrTiO3.TBc111 = uicontrol('Style','radiobutton','String','cubic (111)',...
    'Units','normalized','pos',[0,0.15,1,0.2],'parent',plotdata.Substrate.SrTiO3.TBorientation);
uicontrol(plotdata.Substrate.SrTiO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateSrTiO3c111,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.SrTiO3.TBorientation,'SelectionChangeFcn',@toggleSubstrateSrTiO3orientation);

%%TbScO3
plotdata.Substrate.TbScO3.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','TbScO3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create one button in the button group.
plotdata.Substrate.TbScO3.TBortho001 = uicontrol('Style','radiobutton','String','orthorhombic (001)',...
    'Units','normalized','pos',[0,0.75,1,0.2],'parent',plotdata.Substrate.TbScO3.TBorientation);
uicontrol(plotdata.Substrate.TbScO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateTbScO3o001,4)),'A'],'Units','normalized','pos',[0.2,0.7,0.8,0.1]);

%%TiO2
plotdata.Substrate.TiO2.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','TiO2 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create one button in the button group.
plotdata.Substrate.TiO2.TBt001 = uicontrol('Style','radiobutton','String','tetragonal (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.TiO2.TBorientation);
uicontrol(plotdata.Substrate.TiO2.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateTiO2t001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);

%%YAlO3
plotdata.Substrate.YAlO3.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','YAlO3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create two buttons in the button group.
plotdata.Substrate.YAlO3.TBpc001 = uicontrol('Style','radiobutton','String','pseudo-cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.YAlO3.TBorientation);
uicontrol(plotdata.Substrate.YAlO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateYAlO3pc001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);
plotdata.Substrate.YAlO3.TBpc111 = uicontrol('Style','radiobutton','String','pseudo-cubic (111)',...
    'Units','normalized','pos',[0,0.15,1,0.2],'parent',plotdata.Substrate.YAlO3.TBorientation);
uicontrol(plotdata.Substrate.YAlO3.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateYAlO3pc111,4)),'A'],'Units','normalized','pos',[0.2,0.1,0.8,0.1]);
set(plotdata.Substrate.YAlO3.TBorientation,'SelectionChangeFcn',@toggleSubstrateYAlO3orientation);

%%YSZ (yttria stabilized zirconia ZrO2 including 9.5mol% Y2O3)
plotdata.Substrate.YSZ.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','YSZ (yttria stabilized zirconia ZrO2 incl. 9.5mol%Y2O3 :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create one button in the button group.
plotdata.Substrate.YSZ.TBc001 = uicontrol('Style','radiobutton','String','cubic (001)',...
    'Units','normalized','pos',[0,0.65,1,0.2],'parent',plotdata.Substrate.YSZ.TBorientation);
uicontrol(plotdata.Substrate.YSZ.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateYSZc001,4)),'A'],'Units','normalized','pos',[0.2,0.6,0.8,0.1]);

%%ZnO
plotdata.Substrate.ZnO.TBorientation=uibuttongroup(plotdata.Substrate.Panel,'visible','off',...
    'Title','ZnO :',...
    'BorderType','etchedout',...
    'position',[0,0,1,0.75]);  
% Create three buttons in the button group.
plotdata.Substrate.ZnO.TBh0001 = uicontrol('Style','radiobutton','String','hexagonal (0001)',...
    'Units','normalized','pos',[0,0.8,1,0.15],'parent',plotdata.Substrate.ZnO.TBorientation);
uicontrol(plotdata.Substrate.ZnO.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateZnOh0001,4)),'A'],'Units','normalized','pos',[0.2,0.7,0.8,0.15]);
plotdata.Substrate.ZnO.TBh11bar20 = uicontrol('Style','radiobutton','String','hexagonal (11-20)',...
    'Units','normalized','pos',[0,0.45,1,0.15],'parent',plotdata.Substrate.ZnO.TBorientation);
uicontrol(plotdata.Substrate.ZnO.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateZnOh11bar20,4)),'A'],'Units','normalized','pos',[0.2,0.35,0.8,0.15]);
plotdata.Substrate.ZnO.TBh10bar10 = uicontrol('Style','radiobutton','String','hexagonal (10-10)',...
    'Units','normalized','pos',[0,0.1,1,0.15],'parent',plotdata.Substrate.ZnO.TBorientation);
uicontrol(plotdata.Substrate.ZnO.TBorientation,'style','text','string',['d=',num2str(round(dsubstrateZnOh10bar10,4)),'A'],'Units','normalized','pos',[0.2,0,0.8,0.15]);
set(plotdata.Substrate.ZnO.TBorientation,'SelectionChangeFcn',@toggleSubstrateZnOorientation);

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
      'String',num2str(fitdata.Repetition.N),...
      'Units','normalized',...
      'Position',[0.5,0.945,0.1,0.05],...
      'callback', [plotdata.prog,' writeRepetition']);
   plotdata.RepetitionSlide = uicontrol(plotdata.Superlattice.Panel,...
      'Style','slider',...
      'Min' ,plotdata.Repetition.Min,'Max',plotdata.Repetition.Max, ...
      'Units','normalized',...
      'Position',[0.6,0.94,0.3,0.05], ...
      'Value', fitdata.Repetition.N,...
      'SliderStep',[1/plotdata.Repetition.Max 1/plotdata.Repetition.Max], ...
      'CallBack', [plotdata.prog,' slideRepetition']);
   
  %==Material chooser and slider ========================= 
  for k=1:6,
      
   plotdata.Material(k).Panel=uipanel(plotdata.MainPanel.fig,... 
       'BorderType','etchedout',...
       'position',[0.17+(k-1)*0.135,0.01,0.13,0.84]);
    %'position',[0.17+(k-1)*0.165,0.01,0.16,0.84]);
   plotdata.Material(k).choose=uicontrol(plotdata.Material(k).Panel,'Style','popupmenu', ...
      'String','none|AlO2|BaBiO3|BaO|BaSnO3|(Ba_x,Sr_{1-x})TiO3|BaTiO3|BiFeO3|CaCuO2|Ca2RuO4|CaTiO3|CaVO3|CoO|LaAlO3|LaCoO3|La2CuO4|LaFeO3|LaMnO3|LaNiO3|La2NiMnO6|LaO|LaTiO3|LaVO3|LSMO (La0.67Sr0.33MnO3)|Mg3N2|MgO|MnO|MnO2|MnTiO3|(Nd_x,La_{1-x})NiO3|NdNiO2|NdNiO3|Nd2NiMnO6|NdO|NiO2|PbO|PbNiO3|(Pb_x,Sr_{1-x})TiO3|PbTiO3|Pb(Zr_x,Ti_{1-x})O3|PrBa2Cu3O7|PrNiO2|PrNiO3|PrVO3|RuO2|SmNiO3|Sr3Al2O6|SrCoO2.5|SrCoO3|SrCrO3|SrCuO2|SrIrO3|SrMoO3|SrO|SrO2|SrRuO3|SrTiO3|SrVO3|TiO2|Tm3Fe5O12|VO2|YBa2Cu3O7|YBiO3|Y3Fe5O12|(Y_xTm_{3-x})Fe5O12|YNiO3|(Zn_x,Mg_{1-x})O|Zn3N2|ZnO|ZrO2',...      
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
  %Case pseudo-cubic: 
  plotdata.Material(k).choosepcOrientation=uicontrol('Parent',plotdata.Material(k).OrientationPanel,...
      'Style','popupmenu',...
      'String','pc(001)|pc(111)',... 
      'value',1,...
      'visible','off',...
      'Units','normalized',...
      'Position',[0,0,1,1],...
      'CallBack',[plotdata.prog,' choosepcOrientation(',num2str(k),')'] );
  %Case hexagonal: 
  plotdata.Material(k).choosehOrientation=uicontrol('Parent',plotdata.Material(k).OrientationPanel,...
      'Style','popupmenu',...
      'String','h(0001)|h(11-20)|h(10-10)',... 
      'value',1,...
      'visible','off',...
      'Units','normalized',...
      'Position',[0,0,1,1],...
      'CallBack',[plotdata.prog,' choosehOrientation(',num2str(k),')'] );

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
      'String',num2str(fitdata.Material(k).Polarization),...
      'Units','normalized',...
      'Position',[0.7,0.5,0.3,0.6],...
      'callback', [plotdata.prog,' writePolarization(',num2str(k),')']);
  
  plotdata.Material(k).PolarizationSlide = uicontrol(plotdata.Material(k).TBPolarization,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).PolarizationMin,'Max',plotdata.Material(k).PolarizationMax, ...
      'Units','normalized',...
      'Position',[0.7,0,0.3,0.6], ...
      'Value', fitdata.Material(k).Polarization,...
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
      'String',num2str(fitdata.Material(k).xBST),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexBST(',num2str(k),')']);
  
  plotdata.Material(k).xBSTSlide = uicontrol(plotdata.Material(k).xBSTpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xBSTMin,'Max',plotdata.Material(k).xBSTMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', fitdata.Material(k).xBST,...
      'SliderStep',[0.1/(plotdata.Material(k).xBSTMax-plotdata.Material(k).xBSTMin) 0.1/(plotdata.Material(k).xBSTMax-plotdata.Material(k).xBSTMin)], ...
      'CallBack', [plotdata.prog,' slidexBST(',num2str(k),')']);

  %==NLN concentration chooser and slider =========================
  plotdata.Material(k).xNLNpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xNLNWrite = uicontrol(plotdata.Material(k).xNLNpanel,'Style','edit',...
      'String',num2str(fitdata.Material(k).xNLN),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexNLN(',num2str(k),')']);
  
  plotdata.Material(k).xNLNSlide = uicontrol(plotdata.Material(k).xNLNpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xNLNMin,'Max',plotdata.Material(k).xNLNMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', fitdata.Material(k).xNLN,...
      'SliderStep',[0.01/(plotdata.Material(k).xNLNMax-plotdata.Material(k).xNLNMin) 0.01/(plotdata.Material(k).xNLNMax-plotdata.Material(k).xNLNMin)], ...
      'CallBack', [plotdata.prog,' slidexNLN(',num2str(k),')']);
  
  %==PST concentration chooser and slider =========================
  plotdata.Material(k).xPSTpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xPSTWrite = uicontrol(plotdata.Material(k).xPSTpanel,'Style','edit',...
      'String',num2str(fitdata.Material(k).xPST),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexPST(',num2str(k),')']);
  
  plotdata.Material(k).xPSTSlide = uicontrol(plotdata.Material(k).xPSTpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xPSTMin,'Max',plotdata.Material(k).xPSTMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', fitdata.Material(k).xPST,...
      'SliderStep',[0.1/(plotdata.Material(k).xPSTMax-plotdata.Material(k).xPSTMin) 0.1/(plotdata.Material(k).xPSTMax-plotdata.Material(k).xPSTMin)], ...
      'CallBack', [plotdata.prog,' slidexPST(',num2str(k),')']);
  
  %==PZT concentration chooser and slider =========================
  plotdata.Material(k).xPZTpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xPZTWrite = uicontrol(plotdata.Material(k).xPZTpanel,'Style','edit',...
      'String',num2str(fitdata.Material(k).xPZT),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexPZT(',num2str(k),')']);
  
  plotdata.Material(k).xPZTSlide = uicontrol(plotdata.Material(k).xPZTpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xPZTMin,'Max',plotdata.Material(k).xPZTMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', fitdata.Material(k).xPZT,...
      'SliderStep',[0.1/(plotdata.Material(k).xPZTMax-plotdata.Material(k).xPZTMin) 0.1/(plotdata.Material(k).xPZTMax-plotdata.Material(k).xPZTMin)], ...
      'CallBack', [plotdata.prog,' slidexPZT(',num2str(k),')']);
 
  %==YTmIG concentration chooser and slider =========================
  plotdata.Material(k).xYTmIGpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xYTmIGWrite = uicontrol(plotdata.Material(k).xYTmIGpanel,'Style','edit',...
      'String',num2str(fitdata.Material(k).xYTmIG),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexYTmIG(',num2str(k),')']);
  
  plotdata.Material(k).xYTmIGSlide = uicontrol(plotdata.Material(k).xYTmIGpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xYTmIGMin,'Max',plotdata.Material(k).xYTmIGMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', fitdata.Material(k).xYTmIG,...
      'SliderStep',[0.1/(plotdata.Material(k).xYTmIGMax-plotdata.Material(k).xYTmIGMin) 0.1/(plotdata.Material(k).xYTmIGMax-plotdata.Material(k).xYTmIGMin)], ...
      'CallBack', [plotdata.prog,' slidexYTmIG(',num2str(k),')']);
  
    %==ZnMgO concentration chooser and slider =========================
  plotdata.Material(k).xZMOpanel=uibuttongroup(plotdata.Material(k).Panel,...
      'visible','off',...
      'Title','x:',...
      'BorderType','etchedout',...
      'Position',[0,0.7,1,0.1]);

  plotdata.Material(k).xZMOWrite = uicontrol(plotdata.Material(k).xZMOpanel,'Style','edit',...
      'String',num2str(fitdata.Material(k).xZMO),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writexZMO(',num2str(k),')']);
  
  plotdata.Material(k).xZMOSlide = uicontrol(plotdata.Material(k).xZMOpanel,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).xZMOMin,'Max',plotdata.Material(k).xZMOMax, ...
      'Units','normalized',...
      'Position',[0.35,0.3,0.65,0.5], ...
      'Value', fitdata.Material(k).xZMO,...
      'SliderStep',[0.01/(plotdata.Material(k).xZMOMax-plotdata.Material(k).xZMOMin) 0.01/(plotdata.Material(k).xZMOMax-plotdata.Material(k).xZMOMin)], ...
      'CallBack', [plotdata.prog,' slidexZMO(',num2str(k),')']);

   %slider for N
   plotdata.Material(k).chooseN=uipanel(plotdata.Material(k).Panel, ...
      'Title','N:',...
      'BorderType','etchedout',...
      'Position',[0,0.6,1,0.1]);
   plotdata.Material(k).NWrite = uicontrol(plotdata.Material(k).chooseN,'Style','edit',...
      'String',num2str(fitdata.Material(k).N),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writeN(',num2str(k),')']);
   plotdata.Material(k).NSlide = uicontrol(plotdata.Material(k).chooseN,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).NMin,'Max',plotdata.Material(k).NMax, ...
      'Units','normalized',...
      'Position',[0.35,0.35,0.65,0.5], ...
      'Value', fitdata.Material(k).N,...
      'SliderStep',[1/(plotdata.Material(k).NMax-plotdata.Material(k).NMin) 1/(plotdata.Material(k).NMax-plotdata.Material(k).NMin)], ...
      'CallBack', [plotdata.prog,' slideN(',num2str(k),')']);   

%== d: choose between constant, exponential and sinusoidal distribution =========================

plotdata.Material(k).chooseddistrib = uipanel(plotdata.Material(k).Panel,'visible','on',...
    'Title','Choose d distribution:',...
    'BorderType','etchedout',...
    'Position',[0 0 1 0.55]);

fitdata.Material(k).ddistrib = uicontrol(plotdata.Material(k).chooseddistrib,'Style','popupmenu',...
    'String','d=constant|<html>c<sub>z</sub>=Aexp<sup>z/B</sup>+C</html>',...%    'String','d=constant|<html>d<sub>z</sub>=Aexp<sup>z/B</sup>+C</html>|Jacobi elliptic function|upload|',...
    'value',1,...
    'Units','normalized',...
    'visible','on',...
    'pos',[0,0.475,1,0.5],...
    'CallBack',[plotdata.prog,' chooseddistrib(',num2str(k),')'] );
%% d=constant
fitdata.Material(k).dconst=uipanel(plotdata.Material(k).chooseddistrib,...
    'visible','on','Position',[0 0 1 0.825]);
  
uicontrol(fitdata.Material(k).dconst,...
    'Style','text','String','  d:','Units','normalized','Position',[0 0.9 1 0.1]);

fitdata.Material(k).dWrite = uicontrol(fitdata.Material(k).dconst,'Style','edit',...
      'String',num2str(fitdata.Material(k).d),...
      'Units','normalized',...
      'Position',[0.6,0.8,0.4,0.2],...
      'callback', [plotdata.prog,' writed(',num2str(k),')']);
  
fitdata.Material(k).dSlide = uicontrol(fitdata.Material(k).dconst,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).dMin,'Max',plotdata.Material(k).dMax, ...
      'Units','normalized',...
      'Position',[0,0.6,1,0.2], ...
      'Value', fitdata.Material(k).d,...
      'SliderStep',[0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin) 0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin)], ...
      'CallBack', [plotdata.prog,' slided(',num2str(k),')']);
%% d=exp
fitdata.Material(k).dexpdistrib=uipanel(plotdata.Material(k).chooseddistrib,...
      'visible','off',...
      'BorderType','etchedout',...
      'Position',[0,0,1,0.825]);
%A  
plotdata.Material(k).chooseAexp=uipanel(fitdata.Material(k).dexpdistrib,...
    'Title','A:',...
    'BorderType','etchedout',...
    'Position',[0 0.7 1 0.3]);

plotdata.Material(k).AexpWrite = uicontrol(plotdata.Material(k).chooseAexp,'Style','edit',...
      'String',num2str(fitdata.Material(k).Aexp),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writeAexp(',num2str(k),')']);
  
plotdata.Material(k).AexpSlide = uicontrol(plotdata.Material(k).chooseAexp,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).AexpMin,'Max',plotdata.Material(k).AexpMax, ...
      'Units','normalized',...
      'Position',[0.35,0.35,0.65,0.5], ...
      'Value', fitdata.Material(k).Aexp,...
      'SliderStep',[0.001*1/(plotdata.Material(k).AexpMax-plotdata.Material(k).AexpMin) 0.001*1/(plotdata.Material(k).AexpMax-plotdata.Material(k).AexpMin)], ...
      'CallBack', [plotdata.prog,' slideAexp(',num2str(k),')']);
%B
plotdata.Material(k).chooseBexp=uipanel(fitdata.Material(k).dexpdistrib,...
    'Title','B:',...
    'BorderType','etchedout',...
    'Position',[0 0.4 1 0.3]);

plotdata.Material(k).BexpWrite = uicontrol(plotdata.Material(k).chooseBexp,'Style','edit',...
      'String',num2str(fitdata.Material(k).Bexp),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writeBexp(',num2str(k),')']);
  
plotdata.Material(k).BexpSlide = uicontrol(plotdata.Material(k).chooseBexp,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).BexpMin,'Max',plotdata.Material(k).BexpMax,...
      'Units','normalized',...
      'Position',[0.35,0.35,0.65,0.5], ...
      'Value', fitdata.Material(k).Bexp,...
      'SliderStep',[1/(plotdata.Material(k).BexpMax-plotdata.Material(k).BexpMin) 1/(plotdata.Material(k).BexpMax-plotdata.Material(k).BexpMin)], ...
      'CallBack', [plotdata.prog,' slideBexp(',num2str(k),')']);
% average d
plotdata.Material(k).chooseAveraged=uipanel(fitdata.Material(k).dexpdistrib,...
    'Title','Average d:',...
    'BorderType','etchedout',...
    'Position',[0 0.1 1 0.3]);

plotdata.Material(k).MeandexpWrite = uicontrol(plotdata.Material(k).chooseAveraged,'Style','edit',...
      'String',num2str(fitdata.Material(k).Meandexp),...
      'Units','normalized',...
      'Position',[0,0.1,0.3,0.8],...
      'callback', [plotdata.prog,' writeMeandexp(',num2str(k),')']);
  
plotdata.Material(k).MeandexpSlide = uicontrol(plotdata.Material(k).chooseAveraged,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).MeandexpMin,'Max',plotdata.Material(k).MeandexpMax, ...
      'Units','normalized',...
      'Position',[0.35,0.35,0.65,0.5], ...
      'Value', fitdata.Material(k).Meandexp,...
      'SliderStep',[0.001*1/(plotdata.Material(k).MeandexpMax-plotdata.Material(k).MeandexpMin) 0.001*1/(plotdata.Material(k).MeandexpMax-plotdata.Material(k).MeandexpMin)], ...
      'CallBack', [plotdata.prog,' slideMeandexp(',num2str(k),')']);
  
uicontrol(fitdata.Material(k).dexpdistrib,'Style','text', ...
      'String','   C:',...
      'Position',[0,0,22,12]);
  
plotdata.Material(k).CexpWrite=uicontrol(fitdata.Material(k).dexpdistrib,'Style','text', ...
      'String',num2str(fitdata.Material(k).Cexp),...
      'Position',[25,0,50,12]);
%% Jacobi elliptic function
fitdata.Material(k).dJacobidistrib=uipanel(plotdata.Material(k).chooseddistrib,...
    'visible','off','Position',[0 0 1 0.825]);
  
uicontrol(fitdata.Material(k).dJacobidistrib,...
    'Style','text','String','d @center:','Units','normalized','Position',[0 0.8 0.7 0.15]);%[0 0.8 1 0.15]);

fitdata.Material(k).dCenterWrite = uicontrol(fitdata.Material(k).dJacobidistrib,'Style','edit',...
      'String',num2str(fitdata.Material(k).dCenter),...
      'Units','normalized',...
      'Position',[0.6,0.8,0.4,0.2],...
      'callback', [plotdata.prog,' writedCenter(',num2str(k),')']);
  
fitdata.Material(k).dCenterSlide = uicontrol(fitdata.Material(k).dJacobidistrib,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).dMin,'Max',plotdata.Material(k).dMax, ...
      'Units','normalized',...
      'Position',[0,0.6,1,0.2], ...
      'Value', fitdata.Material(k).dCenter,...
      'SliderStep',[0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin) 0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin)], ...
      'CallBack', [plotdata.prog,' slidedCenter(',num2str(k),')']);
  
uicontrol(fitdata.Material(k).dJacobidistrib,...
    'Style','text','String','d @border:','Units','normalized','Position',[0 0.45 0.7 0.15]);%[0 0.8 1 0.15]);

fitdata.Material(k).dBorderWrite = uicontrol(fitdata.Material(k).dJacobidistrib,'Style','edit',...
      'String',num2str(fitdata.Material(k).dBorder),...
      'Units','normalized',...
      'Position',[0.6,0.45,0.4,0.2],...
      'callback', [plotdata.prog,' writedBorder(',num2str(k),')']);
  
fitdata.Material(k).dBorderSlide = uicontrol(fitdata.Material(k).dJacobidistrib,...
      'Style','slider',...
      'Min' ,plotdata.Material(k).dMin,'Max',plotdata.Material(k).dMax, ...
      'Units','normalized',...
      'Position',[0,0.25,1,0.2], ...
      'Value', fitdata.Material(k).dBorder,...
      'SliderStep',[0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin) 0.001*1/(plotdata.Material(k).dMax-plotdata.Material(k).dMin)], ...
      'CallBack', [plotdata.prog,' slidedBorder(',num2str(k),')']);
  
uicontrol(fitdata.Material(k).dJacobidistrib,...
    'Style','text','String','modulus m:','Units','normalized','Position',[0 0.1 0.7 0.15]);

plotdata.Material(k).mWrite = uicontrol(fitdata.Material(k).dJacobidistrib,'Style','edit',...
      'String',num2str(fitdata.Material(k).m),...
      'Units','normalized',...
      'Position',[0.6,0.1,0.4,0.2],...
      'callback', [plotdata.prog,' writem(',num2str(k),')']);
  
plotdata.Material(k).mSlide = uicontrol(fitdata.Material(k).dJacobidistrib,...
      'Style','slider',...
      'Min' ,plotdata.mMin,'Max',plotdata.mMax, ...
      'Units','normalized',...
      'Position',[0,-0.1,1,0.2], ...
      'Value', fitdata.Material(k).m,...
      'SliderStep',[0.001*1/(plotdata.mMax-plotdata.mMin) 0.001*1/(plotdata.mMax-plotdata.mMin)], ...
      'CallBack', [plotdata.prog,' slidem(',num2str(k),')']);
  
  %% upload d-spacing in Angströms
fitdata.Material(k).dupload=uipanel(plotdata.Material(k).chooseddistrib,...
    'visible','off','Position',[0 0 1 0.825]);
  
uicontrol(fitdata.Material(k).dupload,...
    'Style','text','String','Load file containing d spacing in Angströms, number of rows = N, first row = bottom layer:','Units','normalized',...
    'Position',[0 0.5 1 0.5]);%[0 0.8 1 0.15]);

   uicontrol(fitdata.Material(k).dupload,...
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
       'string',[fitdata.FitDisplay.haxis,' min:'],... 
       'position',[0,0,100,20]);
   plotdata.chooseXMin=uicontrol('parent',plotdata.FitDisplay.fig,'Style','edit', ...
      'String',num2str(plotdata.TThetaMin),... 
      'Position',[80,3,40,20],...
      'CallBack',[plotdata.prog,' chooseXMin'] );
  
  %==X max=========================
   plotdata.chooseXMaxText=uicontrol(plotdata.FitDisplay.fig,'style','text',...
       'string',[fitdata.FitDisplay.haxis,' max:'],... 
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
      'String',num2str(fitdata.ScalingIntensity),... 
      'Position',[200,3,40,20],...
      'CallBack',[plotdata.prog,' chooseScalingIntensity'] );

  %==RMS min=========================
   plotdata.chooseRMSTThetaMinText=uicontrol('parent',plotdata.FitDisplay.fig,'style','text',...
       'visible','off',...
       'string',' RMS TThetamin:',... 
       'position',[110,25,100,25]);
%       'position',[plotdata.FitDisplay.Width/2,plotdata.FitDisplay.Height-27,100,20]);
   plotdata.chooseRMSTThetaMin=uicontrol('parent',plotdata.FitDisplay.fig,'Style','edit', ...
       'visible','off',...
       'String',num2str(plotdata.RMS.TThetaMin),... 
       'position',[200,35,40,20],...
       'CallBack',[plotdata.prog,' chooseRMSMinMax'] );
      %'Position',[plotdata.FitDisplay.Width/2+90,plotdata.FitDisplay.Height-24,40,20],...
   plotdata.chooseRMSLMinText=uicontrol('parent',plotdata.FitDisplay.fig,'style','text',...
       'visible','off',...
       'string',' RMS Lmin:',...       
       'position',[110,25,100,25]);
       %'position',[plotdata.FitDisplay.Width/2,plotdata.FitDisplay.Height-27,100,20]);
   plotdata.chooseRMSLMin=uicontrol('parent',plotdata.FitDisplay.fig,'Style','edit', ...
       'visible','off',...
       'String',num2str(plotdata.RMS.LMin),... 
       'position',[200,35,40,20],...
       'CallBack',[plotdata.prog,' chooseRMSMinMax'] );
  %==RMS max=========================
   plotdata.chooseRMSTThetaMaxText=uicontrol('parent',plotdata.FitDisplay.fig,'style','text',...
       'string',' RMS TThetamax:',... 
       'visible','off',...
       'position',[250,25,100,25]);
   plotdata.chooseRMSTThetaMax=uicontrol('parent',plotdata.FitDisplay.fig,'Style','edit', ...
      'String',num2str(plotdata.RMS.TThetaMax),... 
      'visible','off',...
       'position',[340,35,40,20],...
      'CallBack',[plotdata.prog,' chooseRMSMinMax'] );
   plotdata.chooseRMSLMaxText=uicontrol('parent',plotdata.FitDisplay.fig,'style','text',...
       'string',' RMS Lmax:',... 
       'visible','off',...
       'position',[250,25,100,25]);
   plotdata.chooseRMSLMax=uicontrol('parent',plotdata.FitDisplay.fig,'Style','edit', ...
      'String',num2str(plotdata.RMS.LMax),... 
      'visible','off',...
      'position',[340,35,40,20],...
      'CallBack',[plotdata.prog,' chooseRMSMinMax'] );
     
   I=plotdata.Substrate.g.*conj(plotdata.Substrate.g);
   plotdata.fit.y=I/max(I);
   semilogy(plotdata.fit.x.TTheta,plotdata.fit.y*fitdata.ScalingIntensity,'b');    
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
   end
   close all;
   clear all;
%%   
case 'export' 
   if strcmp(fitdata.FitDisplay.haxis,'L')
      x=plotdata.fit.x.L; 
   elseif strcmp(fitdata.FitDisplay.haxis,'2Theta')
      x=plotdata.fit.x.TTheta; 
   else
    warning('Error in export')
   end
   y=plotdata.fit.y*fitdata.ScalingIntensity;
   
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
    save([startfolder,'parameters'],'fitdata');
    %save important parameters
    parameters=[];
    %Ref to measurement file
    if ~strcmp(plotdata.measure.filename,'none')
      parameters=[parameters;{'Measurement file name: '},plotdata.measure.filename];
    end
    %Substrate
    parameters=[parameters;{'Substrate type: '},{char(fitdata.Substrate.Type)};{'Substrate d[A]: '},{fitdata.Substrate.d};{'Substrate N: '},{fitdata.Substrate.N}];
    %Bottom layer
    if ~strcmp(fitdata.Material(1).Type,'none')
        parameters=[parameters;{'Bottom layer type: '},{char(fitdata.Material(1).Type)}];
        parameters=[parameters;{'Bottom layer orientation: '},{fitdata.Material(1).Orientation}];
        if strcmp(fitdata.Material(1).Type,'(Ba_x,Sr_{1-x})TiO3'),
            parameters=[parameters;{'Bottom layer Ba atomic concentration x: '},{fitdata.Material(1).xBST}];
        elseif strcmp(fitdata.Material(1).Type,'(Nd_x,La_{1-x})NiO3'),
            parameters=[parameters;{'Bottom layer Nd atomic concentration x: '},{fitdata.Material(1).xNLN}];            
        elseif strcmp(fitdata.Material(1).Type,'(Pb_x,Sr_{1-x})TiO3'),
            parameters=[parameters;{'Bottom layer Pb atomic concentration x: '},{fitdata.Material(1).xPST}];
        elseif strcmp(fitdata.Material(1).Type,'Pb(Zr_x,Ti_{1-x})O3'),
            parameters=[parameters;{'Bottom layer Zr atomic concentration x: '},{fitdata.Material(1).xPZT}];
        elseif strcmp(fitdata.Material(1).Type,'(Y_xTm_{3-x})Fe5O12'),
            parameters=[parameters;{'Bottom layer Y atomic concentration x: '},{fitdata.Material(1).xYTmIG}];
        elseif strcmp(fitdata.Material(1).Type,'(Zn_x,Mg_{1-x})O'),
            parameters=[parameters;{'Bottom layer Zn atomic concentration x: '},{fitdata.Material(1).xZMO}];
        end;
        if strcmp(fitdata.Material(1).dDistribution,'constant')
                parameters=[parameters;{'Bottom layer d[A]: '},{round(fitdata.Material(1).d,4)}];
        elseif strcmp(fitdata.Material(1).dDistribution,'exp')
            parameters=[parameters;{'Bottom layer d(z)[A]: '},{[num2str(fitdata.Material(1).Aexp),' x exp(z/',num2str(fitdata.Material(1).Bexp),')+',num2str(fitdata.Material(1).Cexp)]};...
                {'Bottom layer <d>(average)[A]: '},fitdata.Material(1).Meandexp];
        else warning('Error in fitdata.Material(1).dDistribution in Export case');
        end;
        parameters=[parameters;{'Bottom layer N: '},{fitdata.Material(1).N}];
    end
    %Number of repetitions (superlattice) 
    if fitdata.Repetition.N~=1,
        parameters=[parameters;{'Number of superlattice repetitions [Material 1+2+3]: '},{fitdata.Repetition.N}];
    end;
    %Material 1+2+3+4
    for k=2:5,
        if ~strcmp(fitdata.Material(k).Type,'none')
            parameters=[parameters;{['Material ',num2str(k-1),' type : ']},{char(fitdata.Material(k).Type)}];
            parameters=[parameters;{['Material ',num2str(k-1),' orientation : ']},{char(fitdata.Material(k).Orientation)}];
            if strcmp(fitdata.Material(k).Type,'(Ba_x,Sr_{1-x})TiO3'),
                parameters=[parameters;{['Material ',num2str(k-1),' Ba atomic concentration x : ']},{fitdata.Material(k).xBST}];
            elseif strcmp(fitdata.Material(k).Type,'(Nd_x,La_{1-x})NiO3'),
            parameters=[parameters;{['Material ',num2str(k-1),' Nd atomic concentration x: ']},{fitdata.Material(k).xNLN}];
            elseif strcmp(fitdata.Material(k).Type,'(Pb_x,Sr_{1-x})TiO3'),
                parameters=[parameters;{['Material ',num2str(k-1),' Pb atomic concentration x : ']},{fitdata.Material(k).xPST}];
            elseif strcmp(fitdata.Material(k).Type,'Pb(Zr_x,Ti_{1-x})O3'),
                parameters=[parameters;{['Material ',num2str(k-1),' Zr atomic concentration x : ']},{fitdata.Material(k).xPZT}];
            elseif strcmp(fitdata.Material(k).Type,'(Y_xTm_{3-x})Fe5O12'),
                parameters=[parameters;{['Material ',num2str(k-1),' Y atomic concentration x: ']},{fitdata.Material(k).xYTmIG}];
            elseif strcmp(fitdata.Material(k).Type,'(Zn_x,Mg_{1-x})O'),
                parameters=[parameters;{['Material ',num2str(k-1),' Zn atomic concentration x: ']},{fitdata.Material(k).xZMO}];
            end;
            if strcmp(fitdata.Material(k).dDistribution,'constant')
                    parameters=[parameters;{['Material ',num2str(k-1),' d[A]: ']},{round(fitdata.Material(k).d,4)}];
        elseif strcmp(fitdata.Material(k).dDistribution,'exp')
            parameters=[parameters;{['Material ',num2str(k-1),' d(z)[A]: ']},{[num2str(fitdata.Material(k).Aexp),' x exp(z/',num2str(fitdata.Material(k).Bexp),')+',num2str(fitdata.Material(k).Cexp)]};...
                {['Material ',num2str(k-1),' <d>(average)[A]: ']},fitdata.Material(k).Meandexp];
        else warning('Error in fitdata.Material(k).dDistribution in Export case');
        end;
        parameters=[parameters;{['Material ',num2str(k-1),' N: ']},{fitdata.Material(k).N}];
        end
    end;
    %Top layer
    if ~strcmp(fitdata.Material(6).Type,'none')
        parameters=[parameters;{'Top layer type: '},{char(fitdata.Material(6).Type)}];
        parameters=[parameters;{'Top layer orientation: '},{char(fitdata.Material(6).Orientation)}];
        if strcmp(fitdata.Material(6).Type,'(Ba_x,Sr_{1-x})TiO3'),
            parameters=[parameters;{'Top layer Ba atomic concentration x: '},{fitdata.Material(6).xBST}];
        elseif strcmp(fitdata.Material(6).Type,'(Nd_x,La_{1-x})NiO3'),
            parameters=[parameters;{'Top layer Nd atomic concentration x: '},{fitdata.Material(6).xNLN}];            
        elseif strcmp(fitdata.Material(6).Type,'(Pb_x,Sr_{1-x})TiO3'),
            parameters=[parameters;{'Top layer Pb atomic concentration x: '},{fitdata.Material(6).xPST}];
        elseif strcmp(fitdata.Material(6).Type,'Pb(Zr_x,Ti_{1-x})O3'),
            parameters=[parameters;{'Top layer Zr atomic concentration x: '},{fitdata.Material(6).xPZT}];
        elseif strcmp(fitdata.Material(6).Type,'(Y_xTm_{3-x})Fe5O12'),
            parameters=[parameters;{'Top layer Y atomic concentration x: '},{fitdata.Material(6).xYTmIG}]; 
        elseif strcmp(fitdata.Material(6).Type,'(Zn_x,Mg_{1-x})O'),
            parameters=[parameters;{'Top layer Zn atomic concentration x: '},{fitdata.Material(6).xZMO}];
        end;
        if strcmp(fitdata.Material(6).dDistribution,'constant')
                parameters=[parameters;{'Top layer d[A]: '},{round(fitdata.Material(6).d,4)}];
        elseif strcmp(fitdata.Material(6).dDistribution,'exp')
            parameters=[parameters;{'Top layer d(z)[A]: '},{[num2str(fitdata.Material(6).Aexp),' x exp(z/',num2str(fitdata.Material(6).Bexp),')+',num2str(fitdata.Material(6).Cexp)]};...
                {'Top layer <d>(average)[A]: '},fitdata.Material(6).Meandexp];
        else warning('Error in fitdata.Material(6).dDistribution in Export case');
        end;
        parameters=[parameters;{'Top layer N: '},{fitdata.Material(6).N}];
    end
    parameters=cell2table(parameters,'VariableNames',{'Parameters','Values'});
    writetable(parameters,[startfolder,'parameters.csv']);
   end
%%   
case {'slidePolarization(1)','slidePolarization(2)','slidePolarization(3)','slidePolarization(4)','slidePolarization(5)','slidePolarization(6)'}
   k=str2num(entry(19));
   fitdata.Material(k).Polarization=get(plotdata.Material(k).PolarizationSlide,'Value');
   set(plotdata.Material(k).PolarizationWrite,'String',fitdata.Material(k).Polarization);
   set(plotdata.Material(k).TBPolarization,'SelectedObject',plotdata.Material(k).TBPolarizationYes);
   ProgcNSimu
%%     
case {'writePolarization(1)','writePolarization(2)','writePolarization(3)','writePolarization(4)','writePolarization(5)','writePolarization(6)'}
   k=str2num(entry(19));
   fitdata.Material(k).Polarization=str2num(get(plotdata.Material(k).PolarizationWrite,'string'));
   if fitdata.Material(k).Polarization>plotdata.Material(k).PolarizationMax,
       fitdata.Material(k).Polarization=plotdata.Material(k).PolarizationMax;
       set(plotdata.Material(k).PolarizationWrite,'String',fitdata.Material(k).Polarization);
   elseif fitdata.Material(k).Polarization<plotdata.Material(k).PolarizationMin,
       fitdata.Material(k).Polarization=plotdata.Material(k).PolarizationMin;
       set(plotdata.Material(k).PolarizationWrite,'String',fitdata.Material(k).Polarization);
   end;
   set(plotdata.Material(k).PolarizationSlide,'value',fitdata.Material(k).Polarization);
   set(plotdata.Material(k).TBPolarization,'SelectedObject',plotdata.Material(k).TBPolarizationYes);
   ProgcNSimu
%%   
case {'slidexBST(1)','slidexBST(2)','slidexBST(3)','slidexBST(4)','slidexBST(5)','slidexBST(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).xBST=get(plotdata.Material(k).xBSTSlide,'Value');
   set(plotdata.Material(k).xBSTWrite,'String',fitdata.Material(k).xBST);
   ProgcNSimu
%%     
case {'writexBST(1)','writexBST(2)','writexBST(3)','writexBST(4)','writexBST(5)','writexBST(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).xBST=str2num(get(plotdata.Material(k).xBSTWrite,'string'));
   if fitdata.Material(k).xBST>plotdata.Material(k).xBSTMax,
       fitdata.Material(k).xBST=plotdata.Material(k).xBSTMax;
       set(plotdata.Material(k).xBSTWrite,'String',fitdata.Material(k).xBST);
   elseif fitdata.Material(k).xBST<plotdata.Material(k).xBSTMin,
       fitdata.Material(k).xBST=plotdata.Material(k).xBSTMin;
       set(plotdata.Material(k).xBSTWrite,'String',fitdata.Material(k).xBST);
   end;
   set(plotdata.Material(k).xBSTSlide,'value',fitdata.Material(k).xBST);
   ProgcNSimu
%%  
case {'slidexNLN(1)','slidexNLN(2)','slidexNLN(3)','slidexNLN(4)','slidexNLN(5)','slidexNLN(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).xNLN=get(plotdata.Material(k).xNLNSlide,'Value');
   set(plotdata.Material(k).xNLNWrite,'String',fitdata.Material(k).xNLN);
   ProgcNSimu
%%     
case {'writexNLN(1)','writexNLN(2)','writexNLN(3)','writexNLN(4)','writexNLN(5)','writexNLN(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).xNLN=str2num(get(plotdata.Material(k).xNLNWrite,'string'));
   if fitdata.Material(k).xNLN>plotdata.Material(k).xNLNMax,
       fitdata.Material(k).xNLN=plotdata.Material(k).xNLNMax;
       set(plotdata.Material(k).xNLNWrite,'String',fitdata.Material(k).xNLN);
   elseif fitdata.Material(k).xNLN<plotdata.Material(k).xNLNMin,
       fitdata.Material(k).xNLN=plotdata.Material(k).xNLNMin;
       set(plotdata.Material(k).xNLNWrite,'String',fitdata.Material(k).xNLN);
   end;
   set(plotdata.Material(k).xNLNSlide,'value',fitdata.Material(k).xNLN);
   ProgcNSimu
%%   
case {'slidexPST(1)','slidexPST(2)','slidexPST(3)','slidexPST(4)','slidexPST(5)','slidexPST(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).xPST=get(plotdata.Material(k).xPSTSlide,'Value');
   set(plotdata.Material(k).xPSTWrite,'String',fitdata.Material(k).xPST);
   ProgcNSimu
%%     
case {'writexPST(1)','writexPST(2)','writexPST(3)','writexPST(4)','writexPST(5)','writexPST(5)'}
   k=str2num(entry(11));
   fitdata.Material(k).xPST=str2num(get(plotdata.Material(k).xPSTWrite,'string'));
   if fitdata.Material(k).xPST>plotdata.Material(k).xPSTMax,
       fitdata.Material(k).xPST=plotdata.Material(k).xPSTMax;
       set(plotdata.Material(k).xPSTWrite,'String',fitdata.Material(k).xPST);
   elseif fitdata.Material(k).xPST<plotdata.Material(k).xPSTMin,
       fitdata.Material(k).xPST=plotdata.Material(k).xPSTMin;
       set(plotdata.Material(k).xPSTWrite,'String',fitdata.Material(k).xPST);
   end;
   set(plotdata.Material(k).xPSTSlide,'value',fitdata.Material(k).xPST);
   ProgcNSimu
%%   
case {'slidexPZT(1)','slidexPZT(2)','slidexPZT(3)','slidexPZT(4)','slidexPZT(5)','slidexPZT(5)'}
   k=str2num(entry(11));
   fitdata.Material(k).xPZT=get(plotdata.Material(k).xPZTSlide,'Value');
   set(plotdata.Material(k).xPZTWrite,'String',fitdata.Material(k).xPZT);
   ProgcNSimu
%%     
case {'writexPZT(1)','writexPZT(2)','writexPZT(3)','writexPZT(4)','writexPZT(5)','writexPZT(5)'}
   k=str2num(entry(11));
   fitdata.Material(k).xPZT=str2num(get(plotdata.Material(k).xPZTWrite,'string'));
   if fitdata.Material(k).xPZT>plotdata.Material(k).xPZTMax,
       fitdata.Material(k).xPZT=plotdata.Material(k).xPZTMax;
       set(plotdata.Material(k).xPZTWrite,'String',fitdata.Material(k).xPZT);
   elseif fitdata.Material(k).xPZT<plotdata.Material(k).xPZTMin,
       fitdata.Material(k).xPZT=plotdata.Material(k).xPZTMin;
       set(plotdata.Material(k).xPZTWrite,'String',fitdata.Material(k).xPZT);
   end;
   set(plotdata.Material(k).xPZTSlide,'value',fitdata.Material(k).xPZT);
   ProgcNSimu
%%   
case {'slidexYTmIG(1)','slidexYTmIG(2)','slidexYTmIG(3)','slidexYTmIG(4)','slidexYTmIG(5)','slidexYTmIG(5)'}
   k=str2num(entry(13));
   fitdata.Material(k).xYTmIG=get(plotdata.Material(k).xYTmIGSlide,'Value');
   set(plotdata.Material(k).xYTmIGWrite,'String',fitdata.Material(k).xYTmIG);
   ProgcNSimu
%%     
case {'writexYTmIG(1)','writexYTmIG(2)','writexYTmIG(3)','writexYTmIG(4)','writexYTmIG(5)','writexYTmIG(6)'}
   k=str2num(entry(13));
   fitdata.Material(k).xYTmIG=str2num(get(plotdata.Material(k).xYTmIGWrite,'string'));
   if fitdata.Material(k).xYTmIG>plotdata.Material(k).xYTmIGMax,
       fitdata.Material(k).xYTmIG=plotdata.Material(k).xYTmIGMax;
       set(plotdata.Material(k).xYTmIGWrite,'String',fitdata.Material(k).xYTmIG);
   elseif fitdata.Material(k).xYTmIG<plotdata.Material(k).xYTmIGMin,
       fitdata.Material(k).xYTmIG=plotdata.Material(k).xYTmIGMin;
       set(plotdata.Material(k).xYTmIGWrite,'String',fitdata.Material(k).xYTmIG);
   end;
   set(plotdata.Material(k).xYTmIGSlide,'value',fitdata.Material(k).xYTmIG);
   ProgcNSimu
%%     
case {'slidexZMO(1)','slidexZMO(2)','slidexZMO(3)','slidexZMO(4)','slidexZMO(5)','slidexZMO(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).xZMO=get(plotdata.Material(k).xZMOSlide,'Value');
   set(plotdata.Material(k).xZMOWrite,'String',fitdata.Material(k).xZMO);
   ProgcNSimu
%%     
case {'writexZMO(1)','writexZMO(2)','writexZMO(3)','writexZMO(4)','writexZMO(5)','writexZMO(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).xZMO=str2num(get(plotdata.Material(k).xZMOWrite,'string'));
   if fitdata.Material(k).xZMO>plotdata.Material(k).xZMOMax,
       fitdata.Material(k).xZMO=plotdata.Material(k).xZMOMax;
       set(plotdata.Material(k).xZMOWrite,'String',fitdata.Material(k).xZMO);
   elseif fitdata.Material(k).xZMO<plotdata.Material(k).xZMOMin,
       fitdata.Material(k).xZMO=plotdata.Material(k).xZMOMin;
       set(plotdata.Material(k).xZMOWrite,'String',fitdata.Material(k).xZMO);
   end;
   set(plotdata.Material(k).xZMOSlide,'value',fitdata.Material(k).xZMO);
   ProgcNSimu
%%
case {'slideN(1)','slideN(2)','slideN(3)','slideN(4)','slideN(5)','slideN(6)'}
   k=str2num(entry(8));
   fitdata.Material(k).N=get(plotdata.Material(k).NSlide,'Value');
   set(plotdata.Material(k).NWrite,'String',fitdata.Material(k).N);
   ProgcNSimu
%%
case {'writeN(1)','writeN(2)','writeN(3)','writeN(4)','writeN(5)','writeN(6)'}
   k=str2num(entry(8));
   fitdata.Material(k).N=str2num(get(plotdata.Material(k).NWrite,'string'));
   if fitdata.Material(k).N>plotdata.Material(k).NMax,
       fitdata.Material(k).N=plotdata.Material(k).NMax;
       set(plotdata.Material(k).NWrite,'String',fitdata.Material(k).N);
   elseif fitdata.Material(k).N<plotdata.Material(k).NMin,
       fitdata.Material(k).N=plotdata.Material(k).NMin;
       set(plotdata.Material(k).NWrite,'String',fitdata.Material(k).N);
   end;
   set(plotdata.Material(k).NSlide,'value',fitdata.Material(k).N);
   ProgcNSimu
%%   
case {'slided(1)','slided(2)','slided(3)','slided(4)','slided(5)','slided(6)'}
   k=str2num(entry(8)); 
   fitdata.Material(k).d=get(fitdata.Material(k).dSlide,'Value');
   set(fitdata.Material(k).dWrite,'String',fitdata.Material(k).d);
   ProgcNSimu
%%  
case {'writed(1)','writed(2)','writed(3)','writed(4)','writed(5)','writed(6)'}
   k=str2num(entry(8)); 
   fitdata.Material(k).d=str2num(get(fitdata.Material(k).dWrite,'string'));
   if fitdata.Material(k).d>plotdata.Material(k).dMax,
       fitdata.Material(k).d=plotdata.Material(k).dMax;
       set(fitdata.Material(k).dWrite,'String',fitdata.Material(k).d);
   elseif fitdata.Material(k).d<plotdata.Material(k).dMin,
       fitdata.Material(k).d=plotdata.Material(k).dMin;
       set(fitdata.Material(k).dWrite,'String',fitdata.Material(k).d);
   end;
   set(fitdata.Material(k).dSlide,'value',fitdata.Material(k).d);
   ProgcNSimu
   %%   
case {'slidedCenter(1)','slidedCenter(2)','slidedCenter(3)','slidedCenter(4)','slidedCenter(5)','slidedCenter(6)'}
   k=str2num(entry(14)); 
   fitdata.Material(k).dCenter=get(fitdata.Material(k).dCenterSlide,'Value');
   set(fitdata.Material(k).dCenterWrite,'String',fitdata.Material(k).dCenter);
   ProgcNSimu
%%  
case {'writedCenter(1)','writedCenter(2)','writedCenter(3)','writedCenter(4)','writedCenter(5)','writedCenter(6)'}
   k=str2num(entry(14)); 
   fitdata.Material(k).dCenter=str2num(get(fitdata.Material(k).dCenterWrite,'string'));
   if fitdata.Material(k).dCenter>plotdata.Material(k).dMax,
       fitdata.Material(k).dCenter=plotdata.Material(k).dMax;
       set(fitdata.Material(k).dCenterWrite,'String',fitdata.Material(k).dCenter);
   elseif fitdata.Material(k).dCenter<plotdata.Material(k).dMin,
       fitdata.Material(k).dCenter=plotdata.Material(k).dMin;
       set(fitdata.Material(k).dCenterWrite,'String',fitdata.Material(k).dCenter);
   end;
   set(fitdata.Material(k).dCenterSlide,'value',fitdata.Material(k).dCenter);
   ProgcNSimu
%%
case {'slidedBorder(1)','slidedBorder(2)','slidedBorder(3)','slidedBorder(4)','slidedBorder(5)','slidedBorder(6)'}
   k=str2num(entry(14)); 
   fitdata.Material(k).dBorder=get(fitdata.Material(k).dBorderSlide,'Value');
   set(fitdata.Material(k).dBorderWrite,'String',fitdata.Material(k).dBorder);
   ProgcNSimu
%%  
case {'writedBorder(1)','writedBorder(2)','writedBorder(3)','writedBorder(4)','writedBorder(5)','writedBorder(6)'}
   k=str2num(entry(14)); 
   fitdata.Material(k).dBorder=str2num(get(fitdata.Material(k).dBorderWrite,'string'));
   if fitdata.Material(k).dBorder>plotdata.Material(k).dMax,
       fitdata.Material(k).dBorder=plotdata.Material(k).dMax;
       set(fitdata.Material(k).dBorderWrite,'String',fitdata.Material(k).dBorder);
   elseif fitdata.Material(k).dBorder<plotdata.Material(k).dMin,
       fitdata.Material(k).dBorder=plotdata.Material(k).dMin;
       set(fitdata.Material(k).dBorderWrite,'String',fitdata.Material(k).dBorder);
   end;
   set(fitdata.Material(k).dBorderSlide,'value',fitdata.Material(k).dBorder);
   ProgcNSimu
%%
case {'slideAexp(1)','slideAexp(2)','slideAexp(3)','slideAexp(4)','slideAexp(5)','slideAexp(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).Aexp=get(plotdata.Material(k).AexpSlide,'Value');
   set(plotdata.Material(k).AexpWrite,'String',fitdata.Material(k).Aexp);
   ProgcNSimu
   k=str2num(entry(11)); 
   set(plotdata.Material(k).CexpWrite,'String',fitdata.Material(k).Cexp);
%%   
case {'writeAexp(1)','writeAexp(2)','writeAexp(3)','writeAexp(4)','writeAexp(5)','writeAexp(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).Aexp=str2num(get(plotdata.Material(k).AexpWrite,'string'));
   if fitdata.Material(k).Aexp>plotdata.Material(k).AexpMax,
       fitdata.Material(k).Aexp=plotdata.Material(k).AexpMax;
       set(plotdata.Material(k).AexpWrite,'String',fitdata.Material(k).Aexp);
   elseif fitdata.Material(k).Aexp<plotdata.Material(k).AexpMin,
       fitdata.Material(k).Aexp=plotdata.Material(k).AexpMin;
       set(plotdata.Material(k).AexpWrite,'String',fitdata.Material(k).Aexp);
   end;
   set(plotdata.Material(k).AexpSlide,'value',fitdata.Material(k).Aexp);
   ProgcNSimu
   k=str2num(entry(11));
   set(plotdata.Material(k).CexpWrite,'String',num2str(fitdata.Material(k).Cexp));
%%   
case {'slideBexp(1)','slideBexp(2)','slideBexp(3)','slideBexp(4)','slideBexp(5)','slideBexp(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).Bexp=get(plotdata.Material(k).BexpSlide,'Value');
   set(plotdata.Material(k).BexpWrite,'String',fitdata.Material(k).Bexp);
   ProgcNSimu
   k=str2num(entry(11));
   set(plotdata.Material(k).CexpWrite,'String',num2str(fitdata.Material(k).Cexp));
%%   
case {'writeBexp(1)','writeBexp(2)','writeBexp(3)','writeBexp(4)','writeBexp(5)','writeBexp(6)'}
   k=str2num(entry(11));
   fitdata.Material(k).Bexp=str2num(get(plotdata.Material(k).BexpWrite,'string'));
   if fitdata.Material(k).Bexp>plotdata.Material(k).BexpMax,
       fitdata.Material(k).Bexp=plotdata.Material(k).BexpMax;
       set(plotdata.Material(k).BexpWrite,'String',fitdata.Material(k).Bexp);
   elseif fitdata.Material(k).Bexp<plotdata.Material(k).BexpMin,
       fitdata.Material(k).Bexp=plotdata.Material(k).BexpMin;
       set(plotdata.Material(k).BexpWrite,'String',fitdata.Material(k).Bexp);
   end;
   set(plotdata.Material(k).BexpSlide,'value',fitdata.Material(k).Bexp);
   ProgcNSimu
   k=str2num(entry(11));
   set(plotdata.Material(k).CexpWrite,'String',num2str(fitdata.Material(k).Cexp));
%%   
case {'slideMeandexp(1)','slideMeandexp(2)','slideMeandexp(3)','slideMeandexp(4)','slideMeandexp(5)','slideMeandexp(6)'}
   k=str2num(entry(15));
   fitdata.Material(k).Meandexp=get(plotdata.Material(k).MeandexpSlide,'Value');
   set(plotdata.Material(k).MeandexpWrite,'String',fitdata.Material(k).Meandexp);
   ProgcNSimu
   k=str2num(entry(15));
   set(plotdata.Material(k).CexpWrite,'String',num2str(fitdata.Material(k).Cexp));
%%   
case {'writeMeandexp(1)','writeMeandexp(2)','writeMeandexp(3)','writeMeandexp(4)','writeMeandexp(5)','writeMeandexp(6)'}
   k=str2num(entry(15));
   fitdata.Material(k).Meandexp=str2num(get(plotdata.Material(k).MeandexpWrite,'string'));
   if fitdata.Material(k).Meandexp>plotdata.Material(k).MeandexpMax,
       fitdata.Material(k).Meandexp=plotdata.Material(k).MeandexpMax;
       set(plotdata.Material(k).MeandexpWrite,'String',fitdata.Material(k).Meandexp);
   elseif fitdata.Material(k).Meandexp<plotdata.Material(k).MeandexpMin,
       fitdata.Material(k).Meandexp=plotdata.Material(k).MeandexpMin;
       set(plotdata.Material(k).MeandexpWrite,'String',fitdata.Material(k).Meandexp);
   end;
   set(plotdata.Material(k).MeandexpSlide,'value',fitdata.Material(k).Meandexp);
   ProgcNSimu
   k=str2num(entry(15));
   set(plotdata.Material(k).CexpWrite,'String',num2str(fitdata.Material(k).Cexp));
%%   
case {'slidem(1)','slidem(2)','slidem(3)','slidem(4)','slidem(5)','slidem(6)'}
   k=str2num(entry(8)); 
   fitdata.Material(k).m=get(plotdata.Material(k).mSlide,'Value');
   set(plotdata.Material(k).mWrite,'String',fitdata.Material(k).m);
   ProgcNSimu
%%  
case {'writem(1)','writem(2)','writem(3)','writem(4)','writem(5)','writem(6)'}
   k=str2num(entry(8)); 
   fitdata.Material(k).m=str2num(get(plotdata.Material(k).mWrite,'string'));
   if fitdata.Material(k).m>plotdata.mMax,
       fitdata.Material(k).m=plotdata.mMax;
       set(plotdata.Material(k).mWrite,'String',fitdata.Material(k).m);
   elseif fitdata.Material(k).m<plotdata.mMin,
       fitdata.Material(k).m=plotdata.mMin;
       set(plotdata.Material(k).mWrite,'String',fitdata.Material(k).m);
   end;
   set(plotdata.Material(k).mSlide,'value',fitdata.Material(k).m);
   ProgcNSimu
%%
case 'slideRepetition'
   fitdata.Repetition.N=get(plotdata.RepetitionSlide,'Value');
   set(plotdata.RepetitionWrite,'String',fitdata.Repetition.N);
   ProgcNSimu
%%    
case 'writeRepetition'
   fitdata.Repetition.N=str2num(get(plotdata.RepetitionWrite,'string'));
   if fitdata.Repetition.N>plotdata.Repetition.Max,
       fitdata.Repetition.N=plotdata.Repetition.Max;
       set(plotdata.RepetitionWrite,'String',fitdata.Repetition.N);
   elseif fitdata.Repetition.N<plotdata.Repetition.Min,
       fitdata.Repetition.N=plotdata.Repetition.Min;
       set(plotdata.RepetitionWrite,'String',fitdata.Repetition.N);
   end;
   set(plotdata.RepetitionSlide,'value',fitdata.Repetition.N);
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
    fitdata.Material(k).N=size(data,1);
    set(plotdata.Material(k).NWrite,'String',fitdata.Material(k).N);
    set(plotdata.Material(k).NSlide,'value',fitdata.Material(k).N);
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
    if strcmp(fitdata.FitDisplay.haxis,'L')
        helpdlg('Displays your measurement in red. The two vertical green lines delimitate the region over which the RMS (root mean square of the pair wise differences of the fit and measurement) is calculated. You can modify the limits by changing the limiting values on the Fit display panel.');
        set(plotdata.chooseRMSLMinText,'visible','on');
        set(plotdata.chooseRMSLMin,'visible','on');
        set(plotdata.chooseRMSLMaxText,'visible','on');
        set(plotdata.chooseRMSLMax,'visible','on');
        set(plotdata.chooseRMSTThetaMinText,'visible','off');
        set(plotdata.chooseRMSTThetaMin,'visible','off');
        set(plotdata.chooseRMSTThetaMaxText,'visible','off');
        set(plotdata.chooseRMSTThetaMax,'visible','off');
    elseif strcmp(fitdata.FitDisplay.haxis,'2Theta')
        helpdlg('Displays your measurement in red. The two vertical green lines delimitate the region over which the RMS (root mean square of the pair wise differences of the fit and measurement) is calculated. You can modify the limits by changing the limiting values on the Fit display panel.');
        set(plotdata.chooseRMSTThetaMinText,'visible','on');
        set(plotdata.chooseRMSTThetaMin,'visible','on');
        set(plotdata.chooseRMSTThetaMaxText,'visible','on');
        set(plotdata.chooseRMSTThetaMax,'visible','on');
        set(plotdata.chooseRMSLMinText,'visible','off');
        set(plotdata.chooseRMSLMin,'visible','off');
        set(plotdata.chooseRMSLMaxText,'visible','off');
        set(plotdata.chooseRMSLMax,'visible','off');
    else warning('Error in fitdata.FitDisplay.haxis')
    end;
   end;
          
case 'chooseSubstrate'
   load('Substrates.mat');
   str=get(plotdata.Substrate.choose,'String');
   entry=get(plotdata.Substrate.choose,'value');
   SubstrateType=strsplit(str(entry,:));
   fitdata.Substrate.Type=SubstrateType(1);
   set([plotdata.Substrate.DyScO3.TBorientation,...
       plotdata.Substrate.Gd3Ga5O12.TBorientation,...
       plotdata.Substrate.GdScO3.TBorientation,...
       plotdata.Substrate.KTaO3.TBorientation,...
       plotdata.Substrate.LaAlO3.TBorientation,...
       plotdata.Substrate.LaGaO3.TBorientation,...
       plotdata.Substrate.LSAT.TBorientation,...
       plotdata.Substrate.MgO.TBorientation,...
       plotdata.Substrate.Nb05STO.TBorientation,...
       plotdata.Substrate.NdAlO3.TBorientation,...
       plotdata.Substrate.NdGaO3.TBorientation,...
       plotdata.Substrate.PMN71PT29.TBorientation,...
       plotdata.Substrate.Si.TBorientation,...
       plotdata.Substrate.SrLaAlO4.TBorientation,...
       plotdata.Substrate.SrLaGaO4.TBorientation,...
       plotdata.Substrate.SrTiO3.TBorientation,...
       plotdata.Substrate.TbScO3.TBorientation,...
       plotdata.Substrate.TiO2.TBorientation,...
       plotdata.Substrate.YAlO3.TBorientation,...
       plotdata.Substrate.YSZ.TBorientation,...
       plotdata.Substrate.ZnO.TBorientation],...
       'visible','off');

   switch char(fitdata.Substrate.Type)
       case 'none'
              fitdata.Substrate.Orientation='none';
              fitdata.Substrate.d=0;
              plotdata.Substrate.g=gsubstrateKTaO3c001.*0;
       case 'DyScO3'
              set(plotdata.Substrate.DyScO3.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='o001';
              fitdata.Substrate.d=dsubstrateDyScO3o001;
              plotdata.Substrate.g=gsubstrateDyScO3o001;
       case 'Gd3Ga5O12'
              set(plotdata.Substrate.Gd3Ga5O12.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='c001';
              fitdata.Substrate.d=dsubstrateGd3Ga5O12c001;
              plotdata.Substrate.g=gsubstrateGd3Ga5O12c001;
       case 'GdScO3'
              set(plotdata.Substrate.GdScO3.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='pc001';
              fitdata.Substrate.d=dsubstrateGdScO3pc001;
              plotdata.Substrate.g=gsubstrateGdScO3pc001;
       case 'KTaO3'
              set(plotdata.Substrate.KTaO3.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='c001';
              fitdata.Substrate.d=dsubstrateKTaO3c001;
              plotdata.Substrate.g=gsubstrateKTaO3c001;
       case 'LaAlO3'
              set(plotdata.Substrate.LaAlO3.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='pc001';
              fitdata.Substrate.d=dsubstrateLaAlO3pc001;
              plotdata.Substrate.g=gsubstrateLaAlO3pc001;       
       case 'LaGaO3'              
              set(plotdata.Substrate.LaGaO3.TBorientation,'visible','on');              
              fitdata.Substrate.Orientation='o110';
              fitdata.Substrate.d=dsubstrateLaGaO3o110;
              plotdata.Substrate.g=gsubstrateLaGaO3o110;
       case 'LSAT'
              set(plotdata.Substrate.LSAT.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='pc001';
              fitdata.Substrate.d=dsubstrateLSATpc001;
              plotdata.Substrate.g=gsubstrateLSATpc001;
       case 'MgO'
              set(plotdata.Substrate.MgO.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='c001';
              fitdata.Substrate.d=dsubstrateMgOc001;
              plotdata.Substrate.g=gsubstrateMgOc001;
       case 'Nb:SrTiO30.5%wt'
              set(plotdata.Substrate.Nb05STO.TBorientation,'visible','on');
              fitdata.Substrate.d=dsubstrateNb05STOc001;
              plotdata.Substrate.g=gsubstrateNb05STOc001;
       case 'NdAlO3'
              set(plotdata.Substrate.NdAlO3.TBorientation,'visible','on');
              fitdata.Substrate.d=dsubstrateNdAlO3pc001;
              plotdata.Substrate.g=gsubstrateNdAlO3pc001;
       case 'NdGaO3'              
              set(plotdata.Substrate.NdGaO3.TBorientation,'visible','on');              
              fitdata.Substrate.Orientation='pc001';
              fitdata.Substrate.d=dsubstrateNdGaO3pc001;
              plotdata.Substrate.g=gsubstrateNdGaO3pc001;
       case 'PMN-PT'
              set(plotdata.Substrate.PMN71PT29.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='t001';
              fitdata.Substrate.d=dsubstratePMN71PT29t001;
              plotdata.Substrate.g=gsubstratePMN71PT29t001;
       case 'Si'
              set(plotdata.Substrate.Si.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='c001';
              fitdata.Substrate.d=dsubstrateSic001;
              plotdata.Substrate.g=gsubstrateSic001;
       case 'SrLaAlO4'
              set(plotdata.Substrate.SrLaAlO4.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='t001';
              fitdata.Substrate.d=dsubstrateSrLaAlO4t001;
              plotdata.Substrate.g=gsubstrateSrLaAlO4t001;
       case 'SrLaGaO4'
              set(plotdata.Substrate.SrLaGaO4.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='t001';
              fitdata.Substrate.d=dsubstrateSrLaGaO4t001;
              plotdata.Substrate.g=gsubstrateSrLaGaO4t001;
       case 'SrTiO3'
              set(plotdata.Substrate.SrTiO3.TBorientation,'visible','on');
              fitdata.Substrate.d=dsubstrateSrTiO3c001;
              plotdata.Substrate.g=gsubstrateSrTiO3c001;
       case 'TbScO3'
              set(plotdata.Substrate.TbScO3.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='o001';
              fitdata.Substrate.d=dsubstrateTbScO3o001;
              plotdata.Substrate.g=gsubstrateTbScO3o001;            
       case 'TiO2'
              set(plotdata.Substrate.TiO2.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='t001';
              fitdata.Substrate.d=dsubstrateTiO2t001;
              plotdata.Substrate.g=gsubstrateTiO2t001;  
       case 'YAlO3'
              set(plotdata.Substrate.YAlO3.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='pc001';
              fitdata.Substrate.d=dsubstrateYAlO3pc001;
              plotdata.Substrate.g=gsubstrateYAlO3pc001;
       case 'YSZ'
              set(plotdata.Substrate.YSZ.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='c001';
              fitdata.Substrate.d=dsubstrateYSZc001;
              plotdata.Substrate.g=gsubstrateYSZc001;
       case 'ZnO'
              set(plotdata.Substrate.ZnO.TBorientation,'visible','on');
              fitdata.Substrate.Orientation='h0001';
              fitdata.Substrate.d=dsubstrateZnOh0001;
              plotdata.Substrate.g=gsubstrateZnOh0001;
   end;
   ProgcNSimu
%%   
case {'chooseLayerType(1)','chooseLayerType(2)','chooseLayerType(3)','chooseLayerType(4)','chooseLayerType(5)','chooseLayerType(6)'}
   k=str2num(entry(17));
   str=get(plotdata.Material(k).choose,'String');
   entry=get(plotdata.Material(k).choose,'Value');
   MaterialType=strsplit(str(entry,:));
   fitdata.Material(k).Type=MaterialType(1);
   fitdata.Material(k).Polarization=0;
   set(plotdata.Material(k).xBSTpanel,'Visible','off');
   set(plotdata.Material(k).xNLNpanel,'Visible','off');
   set(plotdata.Material(k).xPSTpanel,'Visible','off');
   set(plotdata.Material(k).xPZTpanel,'Visible','off');
   set(plotdata.Material(k).xYTmIGpanel,'Visible','off');
   set(plotdata.Material(k).xZMOpanel,'Visible','off');
   set(plotdata.Material(k).TBPolarization,'Visible','off');
   set(plotdata.Material(k).OrientationPanel,'Visible','off');
   set(plotdata.Material(k).choosepcOrientation,'Visible','off');
   set(plotdata.Material(k).choosehOrientation,'Visible','off');
 
   if strcmp(fitdata.Material(k).Type,'none'),
       fitdata.Material(k).N=0; fitdata.Material(k).d=0;
       fitdata.Material(k).Orientation='none';
   
   elseif any(strcmp(fitdata.Material(k).Type,{'AlO2','BaO','LaO','MnO','MnO2','NiO2','NdO','PbO','RuO2','SrO','SrO2','TiO2','VO2','ZrO2'})),
       %AO-type or BO2-type monolayer CAREFULL: THIS IS ABSOLUTELY "ARTIFICIAL"
       fitdata.Material(k).d=2; %distance initiale entre 2 plans
       
   elseif strcmp(fitdata.Material(k).Type,'BaBiO3'),
       %perovskite
       fitdata.Material(k).dpc001=4.43647; %initial distance between two planes
       fitdata.Material(k).dpc111=4.43647/sqrt(3); %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       
   elseif strcmp(fitdata.Material(k).Type,'BaSnO3'),
       %perovskite
       fitdata.Material(k).dpc001=4.188634; %initial distance between two planes
       fitdata.Material(k).dpc111=4.188634/sqrt(3); %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       
   elseif strcmp(fitdata.Material(k).Type,'(Ba_x,Sr_{1-x})TiO3'),
       %perovskite
       fitdata.Material(k).dpc001=3.905; %initial distance between two planes (default x=0 i.e. SrTiO3)
       fitdata.Material(k).dpc111=3.905/sqrt(3); %initial distance between two planes (default x=0 i.e. SrTiO3)
       fitdata.Material(k).d=fitdata.Material(k).dpc001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xBSTpanel,'Visible','on');
 
    elseif strcmp(fitdata.Material(k).Type,'BaTiO3'),
       %perovskite
       fitdata.Material(k).dpc001=4.036; %initial distance between two planes
       fitdata.Material(k).dpc111=4.036/sqrt(3); %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
    
    elseif strcmp(fitdata.Material(k).Type,'BiFeO3'),
       %perovskite
       fitdata.Material(k).dpc001=4.924; %initial distance between two planes
       fitdata.Material(k).dpc111=4.924/sqrt(3); %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       
    elseif strcmp(fitdata.Material(k).Type,{'CaCuO2'}),
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       fitdata.Material(k).d=3.20546; %initial distance between two planes
    
    elseif strcmp(fitdata.Material(k).Type,{'Ca2RuO4'}),
       %pseudotetragonal
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       fitdata.Material(k).d=11.968145; %initial distance between two planes
       
    elseif strcmp(fitdata.Material(k).Type,'CaTiO3'),
       %perovskite
       fitdata.Material(k).dpc001=3.889; %initial distance between two planes
       fitdata.Material(k).dpc111=3.889/sqrt(3); %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       
    elseif strcmp(fitdata.Material(k).Type,'CaVO3'),
       %perovskite
       fitdata.Material(k).dpc001=3.830; %initial distance between two planes
       fitdata.Material(k).dpc111=3.830/sqrt(3); %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
          
    elseif strcmp(fitdata.Material(k).Type,'CoO'),
       fitdata.Material(k).dpc001=4.263; %initial distance between two planes
       fitdata.Material(k).dpc111=4.263/sqrt(3); %initial distance between two planes 
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       
    elseif strcmp(fitdata.Material(k).Type,'LaAlO3'),
       fitdata.Material(k).dpc001=3.787; %initial distance between two planes        
       fitdata.Material(k).dpc111=3.787/sqrt(3); %initial distance between two planes 
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
 
    elseif strcmp(fitdata.Material(k).Type,'LaCoO3'),
       fitdata.Material(k).dpc001=3.816; %initial distance between two planes
       fitdata.Material(k).dpc111=3.816/sqrt(3); %initial distance between two planes 
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
    
    elseif strcmp(fitdata.Material(k).Type,'La2CuO4'),
       fitdata.Material(k).d=5.398; %initial distance between two planes  
       set(plotdata.Material(k).OrientationPanel,'visible','off');
 
    elseif strcmp(fitdata.Material(k).Type,'LaFeO3'),
       fitdata.Material(k).dpc001=3.959; %initial distance between two planes 
       fitdata.Material(k).dpc111=3.959/sqrt(3); %initial distance between two planes 
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';   
       
    elseif strcmp(fitdata.Material(k).Type,'LaMnO3'),
       fitdata.Material(k).dpc001=3.945; %initial distance between two planes
       fitdata.Material(k).dpc111=3.945/sqrt(3); %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';         
  
    elseif strcmp(fitdata.Material(k).Type,'LaNiO3'),
       fitdata.Material(k).dpc001=3.857; %initial distance between two planes
       fitdata.Material(k).dpc111=3.857/sqrt(3); %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';  
   
    elseif strcmp(fitdata.Material(k).Type,'La2NiMnO6'),%doubleperovskite
       fitdata.Material(k).dpc001=3.871; %initial distance between two planes
       fitdata.Material(k).dpc111=3.871/sqrt(3)*2; %initial distance between two planes 
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
    
    elseif strcmp(fitdata.Material(k).Type,'LaTiO3'),
       fitdata.Material(k).dpc001=3.959; %initial distance between two planes https://materialsproject.org/materials/mp-8020/
       fitdata.Material(k).dpc111=3.959/sqrt(3); %initial distance between two planes https://materialsproject.org/materials/mp-8020/
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)'; 
 
    elseif strcmp(fitdata.Material(k).Type,'LaVO3'),
       fitdata.Material(k).dpc001=3.951; %initial distance between two planes https://materialsproject.org/materials/mp-19053/
       fitdata.Material(k).dpc111=3.951/sqrt(3); %initial distance between two planes https://materialsproject.org/materials/mp-19053/
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)'; 
       
    elseif strcmp(fitdata.Material(k).Type,'LSMO'),
       fitdata.Material(k).dpc001=3.873; %initial distance between two planes http://ematweb.cmi.ua.ac.be/emat/pdf/1214.pdf
       fitdata.Material(k).dpc111=3.873/sqrt(3); %initial distance between two planes http://ematweb.cmi.ua.ac.be/emat/pdf/1214.pdf
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)'; 
       
    elseif strcmp(fitdata.Material(k).Type,'Mg3N2'),
       fitdata.Material(k).dpc001=9.9528; %initial distance between two planes
       fitdata.Material(k).dpc111=9.9528/sqrt(3)/2; %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';  
    
    elseif strcmp(fitdata.Material(k).Type,'MgO'),
       fitdata.Material(k).dpc001=4.212; %initial distance between two planes
       fitdata.Material(k).dpc111=4.212/sqrt(3); %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';  
       
    elseif strcmp(fitdata.Material(k).Type,'MnTiO3'),
       fitdata.Material(k).dpc001=3.832; %initial distance between two planes https://materialsproject.org/materials/mp-19082/ V=112.583 for double unit cell
       fitdata.Material(k).dpc111=3.832/sqrt(3); %initial distance between two planes https://materialsproject.org/materials/mp-19082/ V=112.583 for double unit cell
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)'; 
 
    elseif strcmp(fitdata.Material(k).Type,'(Nd_x,La_{1-x})NiO3'),
       fitdata.Material(k).dpc001=3.857; %initial distance between 2 planes - default x=0: LaNiO3
       fitdata.Material(k).dpc111=3.857/sqrt(3); %initial distance between 2 planes - default x=0: LaNiO3
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xNLNpanel,'Visible','on');
      
   elseif strcmp(fitdata.Material(k).Type,'NdNiO2'),
       % tetragonal a=b=3.962, c=3.268 https://materialsproject.org/materials/mp-31063/
       fitdata.Material(k).d=3.268; %initial distance between 2 planes
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       
    elseif strcmp(fitdata.Material(k).Type,'NdNiO3'),
       fitdata.Material(k).dpc001=3.861; %initial distance between 2 planes - NdNiO3: pseudocubic a=3.861 cubicroot(230.143/4) https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
       fitdata.Material(k).dpc111=3.861/sqrt(3); %initial distance between 2 planes - NdNiO3: pseudocubic a=3.861 cubicroot(230.143/4) https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes       
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
 
    elseif strcmp(fitdata.Material(k).Type,'Nd2NiMnO6'),%doubleperovskite
       fitdata.Material(k).dpc001=3.851; %initial distance between 2 planes https://aip.scitation.org/doi/full/10.1063/1.4906989?ver=pdfcov a=5.4162, b=5.4863, c=7.6718 V=a*b*c 4 unit cells 
       fitdata.Material(k).dpc111=3.851/sqrt(3)*2; %initial distance between 2 planes https://aip.scitation.org/doi/full/10.1063/1.4906989?ver=pdfcov a=5.4162, b=5.4863, c=7.6718 V=a*b*c 4 unit cells 
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes https://aip.scitation.org/doi/full/10.1063/1.4906989?ver=pdfcov a=5.4162, b=5.4863, c=7.6718 V=a*b*c 4 unit cells 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
   
    elseif strcmp(fitdata.Material(k).Type,'PbNiO3'),
       fitdata.Material(k).dpc001=3.811; %initial distance between 2 planes - https://materialsproject.org/materials/mp-974108/ V=55.330
       fitdata.Material(k).dpc111=3.811/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-974108/ V=55.330
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       
    elseif strcmp(fitdata.Material(k).Type,'(Pb_x,Sr_{1-x})TiO3'),
       fitdata.Material(k).dpc001=3.905; %initial distance between 2 planes - default x=0: SrTiO3
       fitdata.Material(k).dpc111=3.905/sqrt(3); %initial distance between 2 planes - default x=0: SrTiO3
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xPSTpanel,'Visible','on');

    elseif strcmp(fitdata.Material(k).Type,'PbTiO3'),
       fitdata.Material(k).dpc001=4.152; %initial distance between 2 planes - known bulk value for PbTiO3 with a=b=3.904
       fitdata.Material(k).dpc111=4.152/sqrt(3); %initial distance between 2 planes - known bulk value for PbTiO3 with a=b=3.904
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).TBPolarization,'Visible','on');
       
    elseif strcmp(fitdata.Material(k).Type,'Pb(Zr_x,Ti_{1-x})O3'),
       fitdata.Material(k).dpc001=4.152; %initial distance between 2 planes - default x=0: PbTiO3
       fitdata.Material(k).dpc111=4.152/sqrt(3); %initial distance between 2 planes - default x=0: PbTiO3
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xPZTpanel,'Visible','on');

    elseif strcmp(fitdata.Material(k).Type,{'PrBa2Cu3O7'}),
       fitdata.Material(k).d=11.916; %initial distance between 2 planes
       set(plotdata.Material(k).OrientationPanel,'visible','off');   
   
   elseif strcmp(fitdata.Material(k).Type,'PrNiO2'),
       % using the lattice parameters of SrCuO2 tetragonal a=b=3.948, c=3.485 https://materialsproject.org/materials/mp-37514/
       fitdata.Material(k).d=3.485; %initial distance between 2 planes
       set(plotdata.Material(k).OrientationPanel,'visible','off');
        
    elseif strcmp(fitdata.Material(k).Type,'PrNiO3'),
       fitdata.Material(k).dpc001=3.872; %initial distance between 2 planes - https://materialsproject.org/materials/mp-22280/ V=232.262 4 unit cells
       fitdata.Material(k).dpc111=3.872/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-22280/ V=232.262 4 unit cells
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';   
       
    elseif strcmp(fitdata.Material(k).Type,'PrVO3'),
       fitdata.Material(k).dpc001=3.936; %initial distance between 2 planes - https://materialsproject.org/materials/mp-1069346/
       fitdata.Material(k).dpc111=3.936/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-1069346/
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';   
       
    elseif strcmp(fitdata.Material(k).Type,'SmNiO3'),
       fitdata.Material(k).dpc001=3.794; %initial distance between 2 planes - https://materialsproject.org/materials/mp-1099668/
       fitdata.Material(k).dpc111=3.794/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-1099668/
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';  
       
   elseif strcmp(fitdata.Material(k).Type,{'Sr3Al2O6'}),
       fitdata.Material(k).d=7.999957; %initial distance between 2 planes https://materialsproject.org/materials/mp-3393/  
       set(plotdata.Material(k).OrientationPanel,'visible','off');
    
   elseif strcmp(fitdata.Material(k).Type,{'SrCoO2.5'}),
       fitdata.Material(k).d=7.8194; %initial distance between 2 planes    
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       
   elseif strcmp(fitdata.Material(k).Type,'SrCoO3'),
       %https://materialsproject.org/materials/mp-505766/ alpha=beta=gamma=90
       %a=3.860; b=3.853; c=b; d=cubicroot(a*b*c)
       fitdata.Material(k).dpc001=3.855; %initial distance between 2 planes
       fitdata.Material(k).dpc111=3.855/sqrt(3); %initial distance between 2 planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)'; 
   
   elseif strcmp(fitdata.Material(k).Type,'SrCrO3'),
       fitdata.Material(k).dpc001=3.8185; %initial distance between 2 planes
       fitdata.Material(k).dpc111=3.8185/sqrt(3); %initial distance between 2 planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)'; 
       
   elseif strcmp(fitdata.Material(k).Type,'SrCuO2'),
       % tetragonal a=b=3.948, c=3.485 https://materialsproject.org/materials/mp-37514/
       fitdata.Material(k).d=3.485; %initial distance between 2 planes
       set(plotdata.Material(k).OrientationPanel,'visible','off');
  
   elseif strcmp(fitdata.Material(k).Type,'SrIrO3'),
       fitdata.Material(k).dpc001=3.998; %initial distance between 2 planes - https://materialsproject.org/materials/mp-1016848/
       fitdata.Material(k).dpc111=3.998/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-1016848/
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)'; 
       
    elseif strcmp(fitdata.Material(k).Type,'SrMoO3'),
       fitdata.Material(k).dpc001=4.082; %initial distance between 2 planes - https://materialsproject.org/materials/mp-18747/
       fitdata.Material(k).dpc111=4.082/sqrt(3); %initial distance between 2 planes - https://materialsproject.org/materials/mp-18747/
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)'; 
       
    elseif strcmp(fitdata.Material(k).Type,'SrRuO3'),
       fitdata.Material(k).dpc001=3.985; %initial distance between 2 planes - https://materialsproject.org/materials/mp-4346/
       fitdata.Material(k).dpc111=fitdata.Material(k).dpc001/sqrt(3); 
       fitdata.Material(k).d=fitdata.Material(k).dpc001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)'; 
       
    elseif strcmp(fitdata.Material(k).Type,'SrTiO3'),
       fitdata.Material(k).dpc001=3.905; %initial distance between 2 planes
       fitdata.Material(k).dpc111=fitdata.Material(k).dpc001/sqrt(3); 
       fitdata.Material(k).d=fitdata.Material(k).dpc001;      
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)'; 
       
     elseif strcmp(fitdata.Material(k).Type,'SrVO3'),
       fitdata.Material(k).dpc001=3.901; %initial distance between 2 planes - cubic - https://materialsproject.org/materials/mp-18717/
       fitdata.Material(k).dpc111=fitdata.Material(k).dpc001/sqrt(3); 
       fitdata.Material(k).d=fitdata.Material(k).dpc001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';  
       
     elseif strcmp(fitdata.Material(k).Type,'Tm3Fe5O12'),
       fitdata.Material(k).dpc001=12.2325/4; %initial distance between 2 planes - known bulk value for Tm3Fe5O12 with a=b=c=12.2325 (Landolt/Bornstein)
       fitdata.Material(k).dpc111=12.2325/4/sqrt(3);
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
  
     elseif strcmp(fitdata.Material(k).Type,{'YBa2Cu3O7'}),
       fitdata.Material(k).d=11.824; %initial distance between 2 planes - tetragonal - https://materialsproject.org/materials/mp-20674/
       set(plotdata.Material(k).OrientationPanel,'visible','off');
    
     elseif strcmp(fitdata.Material(k).Type,'YBiO3'),
       %perovskite
       %https://materialsproject.org/materials/mvc-13598/
       fitdata.Material(k).dpc001=4.408542; %initial distance between two planes
       fitdata.Material(k).dpc111=4.408542/sqrt(3); %initial distance between two planes
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between two planes
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       
   elseif strcmp(fitdata.Material(k).Type,'Y3Fe5O12'),
       fitdata.Material(k).dpc001=12.376/4; %initial distance between 2 planes - known bulk value for Y3Fe5O12 with a=b=c=12.376 (Landolt/Bornstein)
       fitdata.Material(k).dpc111=12.376/4/sqrt(3);
       fitdata.Material(k).d=fitdata.Material(k).dpc001; %initial distance between 2 planes 
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
   
   elseif strcmp(fitdata.Material(k).Type,'(Y_xTm_{3-x})Fe5O12'),
       fitdata.Material(k).dpc001=12.2325/4; %initial distance between two planes (default x=0 i.e. Tm_3Fe5O12)
       fitdata.Material(k).dpc111=12.2325/4/sqrt(3); %initial distance between two planes (default x=0 i.e. Tm_3Fe5O12)
       fitdata.Material(k).d=fitdata.Material(k).dpc001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';
       set(plotdata.Material(k).xYTmIGpanel,'Visible','on');
       
   elseif strcmp(fitdata.Material(k).Type,'YNiO3'),
       %https://materialsproject.org/materials/mvc-15448/
       %a=3.753 b=c=3.759 alpha=89.967 beta=gamma=90 V=53.025
       %d=cubicroot(53.025)=3.757
       fitdata.Material(k).d=3.757; %initial distance between 2 planes
       fitdata.Material(k).dpc111=fitdata.Material(k).dpc001/sqrt(3); 
       fitdata.Material(k).d=fitdata.Material(k).dpc001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosepcOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='pc(001)';    
       
   elseif strcmp(fitdata.Material(k).Type,'(Zn_x,Mg_{1-x})O'),
       % same structure as ZnO but with Mg substituting Zn for Mg<30%
       % ZnO: hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
       % substituting with Mg increases out-of-plane lattice parameter
       fitdata.Material(k).dh0001=5.20380; %initial distance between 2 planes
       fitdata.Material(k).dh11bar20=3.2494/2;
       fitdata.Material(k).dh10bar10=3.2494*sqrt(3)/2;
       fitdata.Material(k).d=fitdata.Material(k).dh0001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosehOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='h(0001)';
       set(plotdata.Material(k).xZMOpanel,'Visible','on');
            
   elseif strcmp(fitdata.Material(k).Type,'Zn3N2'),
       % Ia-3 a=b=c=9.76910 alpha=bet=gamma=90
       fitdata.Material(k).d=9.76910; %distance initiale entre 2 plans
       set(plotdata.Material(k).OrientationPanel,'visible','off');
       
   elseif strcmp(fitdata.Material(k).Type,'ZnO'),
       fitdata.Material(k).dh0001=5.20380; %initial distance between 2 planes
       fitdata.Material(k).dh11bar20=3.2494/2;
       fitdata.Material(k).dh10bar10=3.2494*sqrt(3)/2;
       fitdata.Material(k).d=fitdata.Material(k).dh0001;
       set(plotdata.Material(k).OrientationPanel,'visible','on');
       set(plotdata.Material(k).choosehOrientation,'visible','on','value',1);
       fitdata.Material(k).Orientation='h(0001)';
   end;
   fitdata.Material(k).d=round(fitdata.Material(k).d,4);
   set(fitdata.Material(k).dWrite,'String',fitdata.Material(k).d);
   set(fitdata.Material(k).dSlide,'Value',fitdata.Material(k).d);
   ProgcNSimu
   
   %%   
case {'choosepcOrientation(1)','choosepcOrientation(2)','choosepcOrientation(3)','choosepcOrientation(4)','choosepcOrientation(5)','choosepcOrientation(6)'}
   k=str2num(entry(21));
   str=get(plotdata.Material(k).choosepcOrientation,'String');
   entry=get(plotdata.Material(k).choosepcOrientation,'Value');
   MaterialpcOrientation=strsplit(str(entry,:));
   fitdata.Material(k).Orientation=MaterialpcOrientation(1);
   if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
       fitdata.Material(k).d=round(fitdata.Material(k).dpc001,4);
   elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
       fitdata.Material(k).d=round(fitdata.Material(k).dpc111,4);
   else
        msgbox('Error with orientation');
   end;
   set(fitdata.Material(k).dWrite,'String',fitdata.Material(k).d);
   set(fitdata.Material(k).dSlide,'Value',fitdata.Material(k).d);
   ProgcNSimu
%%
case {'choosehOrientation(1)','choosehOrientation(2)','choosehOrientation(3)','choosehOrientation(4)','choosehOrientation(5)','choosehOrientation(6)'}
   k=str2num(entry(20));
   str=get(plotdata.Material(k).choosehOrientation,'String');
   entry=get(plotdata.Material(k).choosehOrientation,'Value');
   MaterialhOrientation=strsplit(str(entry,:));
   fitdata.Material(k).Orientation=MaterialhOrientation(1);
   if strcmp(fitdata.Material(k).Orientation,'h(0001)'),
              fitdata.Material(k).d=round(fitdata.Material(k).dh0001,4);
   elseif strcmp(fitdata.Material(k).Orientation,'h(11-20)'),
              fitdata.Material(k).d=round(fitdata.Material(k).dh11bar20,4);
   elseif strcmp(fitdata.Material(k).Orientation,'h(10-10)'),
              fitdata.Material(k).d=round(fitdata.Material(k).dh10bar10,4);
   else
        msgbox('Error with orientation');
   end;
   set(fitdata.Material(k).dWrite,'String',fitdata.Material(k).d);
   set(fitdata.Material(k).dSlide,'Value',fitdata.Material(k).d);
   ProgcNSimu
%%
case {'chooseddistrib(1)','chooseddistrib(2)','chooseddistrib(3)','chooseddistrib(4)','chooseddistrib(5)','chooseddistrib(6)'}
   k=str2num(entry(16));
   
   if fitdata.Material(k).ddistrib.Value==1,
    fitdata.Material(k).dDistribution='constant';
    set(fitdata.Material(k).dconst,'visible','on');
    set(fitdata.Material(k).dexpdistrib,'visible','off');
    set(fitdata.Material(k).dJacobidistrib,'visible','off');
    set(fitdata.Material(k).dupload,'visible','off');
   elseif fitdata.Material(k).ddistrib.Value==2,
    fitdata.Material(k).dDistribution='exp';
    set(fitdata.Material(k).dconst,'visible','off');
    set(fitdata.Material(k).dexpdistrib,'visible','on');
    set(fitdata.Material(k).dJacobidistrib,'visible','off');    
    set(fitdata.Material(k).dupload,'visible','off');
   elseif fitdata.Material(k).ddistrib.Value==3,
    fitdata.Material(k).dDistribution='Jacobi';
    set(fitdata.Material(k).dconst,'visible','off');
    set(fitdata.Material(k).dexpdistrib,'visible','off');
    set(fitdata.Material(k).dJacobidistrib,'visible','on');
    set(fitdata.Material(k).dupload,'visible','off');
   elseif fitdata.Material(k).ddistrib.Value==4,
    fitdata.Material(k).dDistribution='upload';
    set(fitdata.Material(k).dconst,'visible','off');
    set(fitdata.Material(k).dexpdistrib,'visible','off');
    set(fitdata.Material(k).dJacobidistrib,'visible','off');
    set(fitdata.Material(k).dupload,'visible','on');
else
    warning('Error in toggle function for dDistribution')
end
ProgcNSimu;
%%
case 'chooseXMin'   
   value=str2num(get(plotdata.chooseXMin,'string'));
   
   if strcmp(fitdata.FitDisplay.haxis,'L')
      if value<plotdata.LMax,
          plotdata.LMin=value;
          axis([plotdata.LMin plotdata.LMax 10^-9 1]);
          set(plotdata.chooseXMin,'String',plotdata.LMin);
      end;        
       
   elseif strcmp(fitdata.FitDisplay.haxis,'2Theta')
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

   if strcmp(fitdata.FitDisplay.haxis,'L')
      if value>plotdata.LMin,
          plotdata.LMax=value;
          axis([plotdata.LMin plotdata.LMax 10^-9 1]);
          set(plotdata.chooseXMax,'String',plotdata.LMax);
      end;        
       
   elseif strcmp(fitdata.FitDisplay.haxis,'2Theta')
      if value>plotdata.TThetaMin,
          plotdata.TThetaMax=value;
          axis([plotdata.TThetaMin plotdata.TThetaMax 10^-9 1]);
          set(plotdata.chooseXMax,'String',plotdata.TThetaMax);
      end; 
   else warning('Error in chooseXMax')
   end;
%%
    case 'FullScale'
   if strcmp(fitdata.FitDisplay.haxis,'L')
      plotdata.LMin=0;
      set(plotdata.chooseXMin,'String',plotdata.LMin);
      plotdata.LMax=2*fitdata.Substrate.d/plotdata.lambda;
      set(plotdata.chooseXMax,'String',plotdata.LMax);
      axis([plotdata.LMin plotdata.LMax 10^-9 1]);
   elseif strcmp(fitdata.FitDisplay.haxis,'2Theta')
      plotdata.TThetaMin=0;
      set(plotdata.chooseXMin,'String',plotdata.TThetaMin);
      plotdata.TThetaMax=180;
      set(plotdata.chooseXMax,'String',plotdata.TThetaMax);
      axis([plotdata.TThetaMin plotdata.TThetaMax 10^-9 1]);
   else warning('Error in FullScale')
   end;
%%
case 'chooseScalingIntensity'
   fitdata.ScalingIntensity=str2num(get(plotdata.chooseScalingIntensity,'string'));
   PlotFitAndData
%%   
   case 'chooseRMSMinMax'   
   delete(plotdata.hVerticalLines);
   if strcmp(fitdata.FitDisplay.haxis,'L')
          plotdata.RMS.LMin=str2num(get(plotdata.chooseRMSLMin,'string'));
          plotdata.RMS.LMax=str2num(get(plotdata.chooseRMSLMax,'string'));
          plotdata.hVerticalLines = line([plotdata.RMS.LMin,plotdata.RMS.LMin],get(plotdata.hAxes,'Ylim'),'Color','green');
          plotdata.hVerticalLines = [plotdata.hVerticalLines line([plotdata.RMS.LMax,plotdata.RMS.LMax],get(plotdata.hAxes,'Ylim'),'Color','green')];   
   elseif strcmp(fitdata.FitDisplay.haxis,'2Theta')
          plotdata.RMS.TThetaMin=str2num(get(plotdata.chooseRMSTThetaMin,'string'));
          plotdata.RMS.TThetaMax=str2num(get(plotdata.chooseRMSTThetaMax,'string'));
          plotdata.hVerticalLines = line([plotdata.RMS.TThetaMin,plotdata.RMS.TThetaMin],get(plotdata.hAxes,'Ylim'),'Color','green');
          plotdata.hVerticalLines = [plotdata.hVerticalLines line([plotdata.RMS.TThetaMax,plotdata.RMS.TThetaMax],get(plotdata.hAxes,'Ylim'),'Color','green')];
   else warning('Error in fitdata.FitDisplay.haxis value');
      end; 
   RMS       
   
end
%===A utility to center the window on the screen============
function pos = centerfig(width,height)

% Find the screen size in pixels
screen_s = get(0,'ScreenSize');
pos = [screen_s(3)/2 - width/2, screen_s(4)/2 - height/2, width, height];

function toggleSubstrateDyScO3orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'(001)o')
    fitdata.Substrate.d=dsubstrateDyScO3o001;
    plotdata.Substrate.g=gsubstrateDyScO3o001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'(110)o')
    fitdata.Substrate.d=dsubstrateDyScO3pc001;
    plotdata.Substrate.g=gsubstrateDyScO3pc001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'(101)o')
    fitdata.Substrate.d=dsubstrateDyScO3pc111;
    plotdata.Substrate.g=gsubstrateDyScO3pc111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateGd3Ga5O12orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'cubic (001)')
    fitdata.Substrate.d=dsubstrateGd3Ga5O12c001;
    plotdata.Substrate.g=gsubstrateGd3Ga5O12c001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'cubic (111)')
    fitdata.Substrate.d=dsubstrateGd3Ga5O12c111;
    plotdata.Substrate.g=gsubstrateGd3Ga5O12c111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateGdScO3orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'pseudo-cubic (001)')
    fitdata.Substrate.d=dsubstrateGdScO3pc001;
    plotdata.Substrate.g=gsubstrateGdScO3pc001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'pseudo-cubic (111)')
    fitdata.Substrate.d=dsubstrateGdScO3pc111;
    plotdata.Substrate.g=gsubstrateGdScO3pc111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateKTaO3orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'cubic (001)')
    fitdata.Substrate.d=dsubstrateKTaO3c001;
    plotdata.Substrate.g=gsubstrateKTaO3c001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'cubic (111)')
    fitdata.Substrate.d=dsubstrateKTaO3c111;
    plotdata.Substrate.g=gsubstrateKTaO3c111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateLaAlO3orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'pseudo-cubic (001)')
    fitdata.Substrate.d=dsubstrateLaAlO3pc001;
    plotdata.Substrate.g=gsubstrateLaAlO3pc001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'pseudo-cubic (111)')
    fitdata.Substrate.d=dsubstrateLaAlO3pc111;
    plotdata.Substrate.g=gsubstrateLaAlO3pc111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateLaGaO3orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'(110)o')
    fitdata.Substrate.d=dsubstrateLaGaO3o110;
    plotdata.Substrate.g=gsubstrateLaGaO3o110;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'(101)o')
    fitdata.Substrate.d=dsubstrateLaGaO3o101;
    plotdata.Substrate.g=gsubstrateLaGaO3o101;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'(001)o')
    fitdata.Substrate.d=dsubstrateLaGaO3o001;
    plotdata.Substrate.g=gsubstrateLaGaO3o001;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateLSATorientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'pseudo-cubic (001)')
    fitdata.Substrate.d=dsubstrateLSATpc001;
    plotdata.Substrate.g=gsubstrateLSATpc001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'pseudo-cubic (111)')
    fitdata.Substrate.d=dsubstrateLSATpc111;
    plotdata.Substrate.g=gsubstrateLSATpc111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateMgOorientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'cubic (001)')
    fitdata.Substrate.d=dsubstrateMgOc001;
    plotdata.Substrate.g=gsubstrateMgOc001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'cubic (111)')
    fitdata.Substrate.d=dsubstrateMgOc111;
    plotdata.Substrate.g=gsubstrateMgOc111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateNb05STOorientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'cubic (001)')
    fitdata.Substrate.d=dsubstrateNb05STOc001;
    plotdata.Substrate.g=gsubstrateNb05STOc001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'cubic (111)')
    fitdata.Substrate.d=dsubstrateNb05STOc111;
    plotdata.Substrate.g=gsubstrateNb05STOc111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateNdAlO3orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'pseudo-cubic (001)')
    fitdata.Substrate.d=dsubstrateNdAlO3pc001;
    plotdata.Substrate.g=gsubstrateNdAlO3pc001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'pseudo-cubic (111)')
    fitdata.Substrate.d=dsubstrateNdAlO3pc111;
    plotdata.Substrate.g=gsubstrateNdAlO3pc111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateNdGaO3orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'(110)o')
    fitdata.Substrate.d=dsubstrateNdGaO3pc001;
    plotdata.Substrate.g=gsubstrateNdGaO3pc001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'(101)o')
    fitdata.Substrate.d=dsubstrateNdGaO3pc111;
    plotdata.Substrate.g=gsubstrateNdGaO3pc111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstratePMN71PT29orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'tetragonal (001)')
    fitdata.Substrate.d=dsubstratePMN71PT29t001;
    plotdata.Substrate.g=gsubstratePMN71PT29t001;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateSiorientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'cubic (001)')
    fitdata.Substrate.d=dsubstrateSic001;
    plotdata.Substrate.g=gsubstrateSic001;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateSrLaAlO4orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'tetragonal (001)')
    fitdata.Substrate.d=dsubstrateSrLaAlO4t001;
    plotdata.Substrate.g=gsubstrateSrLaAlO4t001;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateSrLaGaO4orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'tetragonal (001)')
    fitdata.Substrate.d=dsubstrateSrLaGaO4t001;
    plotdata.Substrate.g=gsubstrateSrLaGaO4t001;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateSrTiO3orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'cubic (001)')
    fitdata.Substrate.d=dsubstrateSrTiO3c001;
    plotdata.Substrate.g=gsubstrateSrTiO3c001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'cubic (111)')
    fitdata.Substrate.d=dsubstrateSrTiO3c111;
    plotdata.Substrate.g=gsubstrateSrTiO3c111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateTbScO3orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'orthorhombic (001)')
    fitdata.Substrate.d=dsubstrateTbScO3o001;
    plotdata.Substrate.g=gsubstrateTbScO3o001;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateTiO2orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'tetragonal (001)')
    fitdata.Substrate.d=dsubstrateTiO2t001;
    plotdata.Substrate.g=gsubstrateTiO2t001;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateYAlO3orientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'pseudo-cubic (001)')
    fitdata.Substrate.d=dsubstrateYAlO3pc001;
    plotdata.Substrate.g=gsubstrateYAlO3pc001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'pseudo-cubic (111)')
    fitdata.Substrate.d=dsubstrateYAlO3pc001;
    plotdata.Substrate.g=gsubstrateYAlO3pc111;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateYSZorientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'cubic (001)')
    fitdata.Substrate.d=dsubstrateYSZc001;
    plotdata.Substrate.g=gsubstrateYSZc001;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggleSubstrateZnOorientation(source,eventdata)
global fitdata plotdata;
load('Substrates.mat');
if strcmp(get(get(source,'SelectedObject'),'String'),'hexagonal (0001)')
    fitdata.Substrate.d=dsubstrateZnOh0001;
    plotdata.Substrate.g=gsubstrateZnOh0001;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'hexagonal (11-20)')
    fitdata.Substrate.d=dsubstrateZnOh11bar20;
    plotdata.Substrate.g=gsubstrateZnOh11bar20;
elseif strcmp(get(get(source,'SelectedObject'),'String'),'hexagonal (10-10)')
    fitdata.Substrate.d=dsubstrateZnOh10bar10;
    plotdata.Substrate.g=gsubstrateZnOh10bar10;
else
    warning('Error in Toggle function (substrate orientation)')
end
ProgcNSimu

function toggle(source,eventdata)
global fitdata plotdata;
if strcmp(get(get(source,'SelectedObject'),'String'),'L')
    fitdata.FitDisplay.haxis='L';
    plotdata.LMin=2*fitdata.Substrate.d*sin(plotdata.TThetaMin/2*pi/180)/plotdata.lambda;
    plotdata.LMax=2*fitdata.Substrate.d*sin(plotdata.TThetaMax/2*pi/180)/plotdata.lambda;
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
    fitdata.FitDisplay.haxis='2Theta';   
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
set(plotdata.chooseXMinText,'String',[fitdata.FitDisplay.haxis,' min:']);
set(plotdata.chooseXMaxText,'String',[fitdata.FitDisplay.haxis,' max:']);
PlotFitAndData
%%
function togglePolarization(source,eventdata)
global fitdata plotdata;
k=str2num(source.Parent.Tag(3));
if strcmp(get(get(source,'SelectedObject'),'String'),'Yes')
elseif strcmp(get(get(source,'SelectedObject'),'String'),'No')
    fitdata.Material(k).Polarization=0;
    set(plotdata.Material(k).PolarizationWrite,'String',fitdata.Material(k).Polarization);
    set(plotdata.Material(k).PolarizationSlide,'value',fitdata.Material(k).Polarization);    
else
    warning('Error in toggle function for Polarization')
end
ProgcNSimu;