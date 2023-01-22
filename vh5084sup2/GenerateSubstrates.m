% Created by Celine Lichtensteiger
% Generate substrates among:
% DyScO3, GdScO3 , KTaO3, LaAlO3, LaSrAlO4, LSAT, NdAlO3, NdGaO3, Si, SrTiO3, TbScO3, YAlO3
%*********************************

tic;
global plotdata;

%% Part to be changed to generate substrate
plotdata.Substrate.Type='LaSrAlO4';
plotdata.orientation='(001)';   % '(001)' or '(111)'
%% Main part of the program generating substrates
fO=AtomicScatteringFactor('O',plotdata.QQ);
switch plotdata.Substrate.Type
    
    case {'DyScO3','GdScO3','KTaO3','LaAlO3','LSAT','NdAlO3','NdGaO3','SrTiO3','YAlO3'} %% Perovskite structure (001 and 111 oriented)
        if strcmp(plotdata.Substrate.Type,'DyScO3'),
            plotdata.Substrate.c=3.9403;
            fAsubstrate=AtomicScatteringFactor('Dy',plotdata.QQ);
            fBsubstrate=AtomicScatteringFactor('Sc',plotdata.QQ);
        elseif strcmp(plotdata.Substrate.Type,'GdScO3'),
            plotdata.Substrate.c=3.9636;
            fAsubstrate=AtomicScatteringFactor('Gd',plotdata.QQ);
            fBsubstrate=AtomicScatteringFactor('Sc',plotdata.QQ);
        elseif strcmp(plotdata.Substrate.Type,'KTaO3'), 
            plotdata.Substrate.c=3.989;
            fAsubstrate=AtomicScatteringFactor('K',plotdata.QQ);
            fBsubstrate=AtomicScatteringFactor('Ta',plotdata.QQ);
        elseif strcmp(plotdata.Substrate.Type,'LaAlO3'), 
            plotdata.Substrate.c=3.789;
            fAsubstrate=AtomicScatteringFactor('La',plotdata.QQ); 
            fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);
        elseif strcmp(plotdata.Substrate.Type,'LSAT'), 
            plotdata.Substrate.c=3.868;
            fAsubstrate=0.3*AtomicScatteringFactor('La',plotdata.QQ)+0.7*AtomicScatteringFactor('Sr',plotdata.QQ);
            fBsubstrate=0.65*AtomicScatteringFactor('Al',plotdata.QQ)+0.35*AtomicScatteringFactor('Ta',plotdata.QQ);
        elseif strcmp(plotdata.Substrate.Type,'NdAlO3'), 
            plotdata.Substrate.c=3.74;
            fAsubstrate=AtomicScatteringFactor('Nd',plotdata.QQ);
            fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);
        elseif strcmp(plotdata.Substrate.Type,'NdGaO3'), 
            plotdata.Substrate.c=3.864;
            fAsubstrate=AtomicScatteringFactor('Nd',plotdata.QQ);
            fBsubstrate=AtomicScatteringFactor('Ga',plotdata.QQ);
        elseif strcmp(plotdata.Substrate.Type,'SrTiO3'),
            plotdata.Substrate.c=3.905;
            fAsubstrate=AtomicScatteringFactor('Sr',plotdata.QQ);
            fBsubstrate=AtomicScatteringFactor('Ti',plotdata.QQ);
        elseif strcmp(plotdata.Substrate.Type,'YAlO3'), 
            plotdata.Substrate.c=3.71;
            fAsubstrate=AtomicScatteringFactor('Y',plotdata.QQ);
            fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);
        else warning('Error in perovskite choice in GenerateSubstrates.m');
        end;  
        plotdata.Substrate.d=plotdata.Substrate.c/sqrt(3);
        switch plotdata.orientation
            case '(001)'
                %001 Atomic positions in the unit cell, expressed using the substrate unit cell vectors:
                zperovskite(1)=0;   %A
                zperovskite(2)=0.5; %B
                zperovskite(3)=0;	%O1
                zperovskite(4)=0.5;	%O2
                zperovskite(5)=0.5;	%O3
                %Structure factor F 
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite,plotdata.Substrate.c);
            case '(111)'
                %111 Atomic positions in the unit cell, expressed using the substrate unit cell vectors:
                zperovskite(1)=0;   %A
                zperovskite(2)=0.5; %B
                zperovskite(3)=0;	%O1
                zperovskite(4)=0;	%O2
                zperovskite(5)=0;	%O3
                %Structure factor F
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite,plotdata.Substrate.d);
            otherwise
                warning('Error in plotdata.orientation in GenerateSubstrates.');
        end
    case 'LaSrGaO4'
        LaSrGaO4 001
        zLaSrGaO4(1)=0;         % 1 Ga - 2 O
        zLaSrGaO4(2)=0.1412;    % 0.5 La - 0.5 Sr
        zLaSrGaO4(3)=0.168;     % 1 O
        zLaSrGaO4(4)=0.332;     % 1 O
        zLaSrGaO4(5)=0.3588;    % 0.5 La - 0.5 Sr
        zLaSrGaO4(6)=0.5;       % Ga - 2 O
        zLaSrGaO4(7)=0.6412;    % 0.5 La - 0.5 Sr
        zLaSrGaO4(8)=0.668;     % 1 O
        zLaSrGaO4(9)=0.832;     % 1 O
        zLaSrGaO4(10)=0.8588;   % 0.5 La - 0.5 Sr
        fLa=AtomicScatteringFactor('La',plotdata.QQ);
        fSr=AtomicScatteringFactor('Sr',plotdata.QQ);
        fGa=AtomicScatteringFactor('Ga',plotdata.QQ);
        plotdata.Substrate.c=12.67;
        Fsubstrate=(fGa+2*fO).*exp(1i*(plotdata.Q.*zLaSrGaO4(1).*plotdata.Substrate.c))...
        +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zLaSrGaO4(2).*plotdata.Substrate.c))...
        +fO.*exp(1i*(plotdata.Q.*zLaSrGaO4(3).*plotdata.Substrate.c))...
        +fO.*exp(1i*(plotdata.Q.*zLaSrGaO4(4).*plotdata.Substrate.c))...
        +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zLaSrGaO4(5).*plotdata.Substrate.c))...
        +(fGa+2*fO).*exp(1i*(plotdata.Q.*zLaSrGaO4(6).*plotdata.Substrate.c))...
        +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zLaSrGaO4(7).*plotdata.Substrate.c))...
        +fO.*exp(1i*(plotdata.Q.*zLaSrGaO4(8).*plotdata.Substrate.c))...
        +fO.*exp(1i*(plotdata.Q.*zLaSrGaO4(9).*plotdata.Substrate.c))...
        +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zLaSrGaO4(10).*plotdata.Substrate.c));
    case 'LaSrAlO4'
        zLaSrAlO4(1)=0;                     % 1 Al - 2 O
        zLaSrAlO4(2)=0.141099973121255;     % 0.5 La - 0.5 Sr
        zLaSrAlO4(3)=0.162699969959049;     % 1 O
        zLaSrAlO4(4)=0.337300030040951;     % 1 O
        zLaSrAlO4(5)=0.358900026878745;     % 0.5 La - 0.5 Sr
        zLaSrAlO4(6)=0.5;                   % 1 Al - 2 O
        zLaSrAlO4(7)=0.641099973121255;     % 0.5 La - 0.5 Sr
        zLaSrAlO4(8)=0.662699969959049;     % 1 O
        zLaSrAlO4(9)=0.837300030040951;     % 1 O
        zLaSrAlO4(10)=0.858900026878745;    % 0.5 La - 0.5 Sr
        fLa=AtomicScatteringFactor('La',plotdata.QQ);
        fSr=AtomicScatteringFactor('Sr',plotdata.QQ);
        fAl=AtomicScatteringFactor('Al',plotdata.QQ);
        plotdata.Substrate.c=12.6377;
        Fsubstrate=(fAl+2*fO).*exp(1i*(plotdata.Q.*zLaSrAlO4(1).*plotdata.Substrate.c))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zLaSrAlO4(2).*plotdata.Substrate.c))...
            +fO.*exp(1i*(plotdata.Q.*zLaSrAlO4(3).*plotdata.Substrate.c))...
            +fO.*exp(1i*(plotdata.Q.*zLaSrAlO4(4).*plotdata.Substrate.c))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zLaSrAlO4(5).*plotdata.Substrate.c))...
            +(fAl+2*fO).*exp(1i*(plotdata.Q.*zLaSrAlO4(6).*plotdata.Substrate.c))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zLaSrAlO4(7).*plotdata.Substrate.c))...
            +fO.*exp(1i*(plotdata.Q.*zLaSrAlO4(8).*plotdata.Substrate.c))...
            +fO.*exp(1i*(plotdata.Q.*zLaSrAlO4(9).*plotdata.Substrate.c))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zLaSrAlO4(10).*plotdata.Substrate.c));
    case 'TbScO3'
        zTbScO3(1)=0;                     % 2 Sc
        zTbScO3(2)=0.0564000252620942;    % 2 O
        zTbScO3(3)=0.25;                  % 2 Tb - 2 O
        zTbScO3(4)=0.443599974737906;     % 2 O
        zTbScO3(5)=0.5;                   % 2 Sc
        zTbScO3(6)=0.556400025262094;     % 2 O
        zTbScO3(7)=0.75;                  % 2 Tb - 2 O
        zTbScO3(8)=0.943599974737906;     % 2 O
        fTb=AtomicScatteringFactor('Tb',plotdata.QQ);
        fSc=AtomicScatteringFactor('Sc',plotdata.QQ);
        plotdata.Substrate.c=7.917;
        Fsubstrate=2*fSc.*exp(1i*(plotdata.Q.*zTbScO3(1).*plotdata.Substrate.c))...
            +2*fO.*exp(1i*(plotdata.Q.*zTbScO3(2).*plotdata.Substrate.c))...
            +(2*fTb+2*fO).*exp(1i*(plotdata.Q.*zTbScO3(3).*plotdata.Substrate.c))...
            +2*fO.*exp(1i*(plotdata.Q.*zTbScO3(4).*plotdata.Substrate.c))...
            +2*fSc.*exp(1i*(plotdata.Q.*zTbScO3(5).*plotdata.Substrate.c))...
            +2*fO.*exp(1i*(plotdata.Q.*zTbScO3(6).*plotdata.Substrate.c))...
            +(2*fTb+2*fO).*exp(1i*(plotdata.Q.*zTbScO3(7).*plotdata.Substrate.c))...
            +2*fO.*exp(1i*(plotdata.Q.*zTbScO3(8).*plotdata.Substrate.c));
    case 'Si' %rotated in-plane so one unit-cell surface of silicon correspond to 2 unit-cell surface of traditional perovskite: we therefor consider only half of the atoms to rescale for the relative intensity
        zSi(1)=0;                     %2Si
        zSi(2)=0.25;                  %2Si
        zSi(3)=0.5;                   %2Si
        zSi(4)=0.75;                  %2Si
        plotdata.Substrate.c=5.4307;
        fSi=AtomicScatteringFactor('Si',plotdata.QQ); 
        Fsubstrate=fSi.*exp(1i*(plotdata.Q.*zSi(1).*plotdata.Substrate.c))...
            +fSi.*exp(1i*(plotdata.Q.*zSi(2).*plotdata.Substrate.c))...
            +fSi.*exp(1i*(plotdata.Q.*zSi(3).*plotdata.Substrate.c))...
            +fSi.*exp(1i*(plotdata.Q.*zSi(4).*plotdata.Substrate.c));
end;

%% Calculations

expzmsubstrate=zeros(size(plotdata.Q));

if strcmp(plotdata.orientation,'(001)'),
    NSubstrate=plotdata.Substrate.N001;
elseif strcmp(plotdata.orientation,'(111)'),
    NSubstrate=plotdata.Substrate.N111; %(111): using more layers when calculating the substrate in (111) to avoid finite size oscillations
else warning('Error in plotdata.orientation in GenerateSubstrates.m'),
end;

h = waitbar(0.2,'Please wait...');
for zmsubstrate=1:NSubstrate,
 expzmsubstrate=expzmsubstrate+exp(1i*plotdata.Q.*zmsubstrate.*plotdata.Substrate.c).*exp(-(plotdata.TotalThickness-zmsubstrate.*plotdata.Substrate.c)./plotdata.mu);
end;
close(h) 

%Scattering amplitude g:
plotdata.Substrate.g=Fsubstrate.*expzmsubstrate;

toc;
if strcmp(plotdata.orientation,'(001)'),
    name=['gsubstrate',num2str(plotdata.Substrate.Type),'001'];  % name used to store the calculated intensity
elseif strcmp(plotdata.orientation,'(111)'),
        name=['gsubstrate',num2str(plotdata.Substrate.Type),'111'];  % name used to store the calculated intensity
else warning('Error in plotdata.orientation in GenerateSubstrates.m');
end;
save('Substrates.mat',name,'-append');