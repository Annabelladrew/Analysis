% Created by Celine Lichtensteiger
% Loads the correct substrate data and intensity
%*********************************

global plotdata;

load('Substrates.mat');

if strcmp(plotdata.orientation,'(001)'),
    
   if strcmp(plotdata.Substrate.Type,'none'), 
       plotdata.Substrate.c=0;
       plotdata.Substrate.g=gsubstrateKTaO3001.*0;
   elseif strcmp(plotdata.Substrate.Type,'DyScO3'), 
       plotdata.Substrate.c=3.9403;
       plotdata.Substrate.g=gsubstrateDyScO3001;
   elseif strcmp(plotdata.Substrate.Type,'GdScO3'), 
       plotdata.Substrate.c=3.9636;
       plotdata.Substrate.g=gsubstrateGdScO3001;
   elseif strcmp(plotdata.Substrate.Type,'KTaO3'), 
       plotdata.Substrate.c=3.989;
       plotdata.Substrate.g=gsubstrateKTaO3001;
   elseif strcmp(plotdata.Substrate.Type,'LaAlO3'),
       plotdata.Substrate.c=3.789;
       plotdata.Substrate.g=gsubstrateLaAlO3001;
   elseif strcmp(plotdata.Substrate.Type,'LaSrAlO4'),
       plotdata.Substrate.c=12.6377;
       plotdata.Substrate.g=gsubstrateLaSrAlO4001;
   elseif strcmp(plotdata.Substrate.Type,'LSAT'),
      plotdata.Substrate.c=3.868;       
      plotdata.Substrate.g=gsubstrateLSAT001;
   elseif strcmp(plotdata.Substrate.Type,'NdAlO3'),
      plotdata.Substrate.c=3.74;       
      plotdata.Substrate.g=gsubstrateNdAlO3001;
   elseif strcmp(plotdata.Substrate.Type,'NdGaO3'),       
      plotdata.Substrate.c=3.864;       
      plotdata.Substrate.g=gsubstrateNdGaO3001;
   elseif strcmp(plotdata.Substrate.Type,'Si'),       
      plotdata.Substrate.c=5.4307;       
      plotdata.Substrate.g=gsubstrateSi001;
   elseif strcmp(plotdata.Substrate.Type,'SrTiO3'),       
      plotdata.Substrate.c=3.905;       
      plotdata.Substrate.g=gsubstrateSrTiO3001;
   elseif strcmp(plotdata.Substrate.Type,'TbScO3'),       
      plotdata.Substrate.c=7.917;       
      plotdata.Substrate.g=gsubstrateTbScO3001;
   elseif strcmp(plotdata.Substrate.Type,'YAlO3'),       
      plotdata.Substrate.c=3.71;       
      plotdata.Substrate.g=gsubstrateYAlO3001;
   else warning('Error in substrate type (beginning of ProgcNSimu)');
   end
   
   elseif strcmp(plotdata.orientation,'(111)'),
       
   if strcmp(plotdata.Substrate.Type,'none'); 
       plotdata.Substrate.d=0;
       plotdata.Substrate.g=gsubstrateKTaO3111.*0;
   elseif strcmp(plotdata.Substrate.Type,'DyScO3'), 
       plotdata.Substrate.d=3.9403/sqrt(3);
       plotdata.Substrate.g=gsubstrateDyScO3111;
   elseif strcmp(plotdata.Substrate.Type,'GdScO3'), 
       plotdata.Substrate.d=3.9636/sqrt(3);
       plotdata.Substrate.g=gsubstrateGdScO3111;
   elseif strcmp(plotdata.Substrate.Type,'KTaO3'), 
       plotdata.Substrate.d=3.989/sqrt(3);
       plotdata.Substrate.g=gsubstrateKTaO3111;
   elseif strcmp(plotdata.Substrate.Type,'LaAlO3'),
       plotdata.Substrate.d=3.789/sqrt(3);
       plotdata.Substrate.g=gsubstrateLaAlO3111;
   elseif strcmp(plotdata.Substrate.Type,'LaSrAlO4'),
      msgbox('LaSrAlO4 is not available in the (111) orientation. If needed, please contact Celine.Lichtensteiger@unige.ch.');
   elseif strcmp(plotdata.Substrate.Type,'LSAT'),
      plotdata.Substrate.d=3.868/sqrt(3);       
      plotdata.Substrate.g=gsubstrateLSAT111;
   elseif strcmp(plotdata.Substrate.Type,'NdAlO3'),
      plotdata.Substrate.d=3.74/sqrt(3);       
      plotdata.Substrate.g=gsubstrateNdAlO3111;
   elseif strcmp(plotdata.Substrate.Type,'NdGaO3'),       
      plotdata.Substrate.d=3.864/sqrt(3);       
      plotdata.Substrate.g=gsubstrateNdGaO3111;
   elseif strcmp(plotdata.Substrate.Type,'Si'),
      msgbox('Si is not available in the (111) orientation. If needed, please contact Celine.Lichtensteiger@unige.ch.');
   elseif strcmp(plotdata.Substrate.Type,'SrTiO3'),       
      plotdata.Substrate.d=3.905/sqrt(3);       
      plotdata.Substrate.g=gsubstrateSrTiO3111;
   elseif strcmp(plotdata.Substrate.Type,'TbScO3'),       
      msgbox('TbScO3 is not available in the (111) orientation. If needed, please contact Celine.Lichtensteiger@unige.ch.');
   elseif strcmp(plotdata.Substrate.Type,'YAlO3'),       
      plotdata.Substrate.d=3.71/sqrt(3);       
      plotdata.Substrate.g=gsubstrateYAlO3111;
   else warning('Error in substrate type (beginning of ProgcNSimu');
   end
   
   else warning('Error in plotdata.orientation'),
   end;