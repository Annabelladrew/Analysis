% Created by Celine Lichtensteiger
% Generate substrates among:
% DyScO3 (o001,pc001,pc111)
% Gd3Ga5O12 (garnet c001, c111)
% GdScO3 (pc001, pc111)
% KTaO3 (c001, c111)
% LaAlO3 (pc001, pc111)
% LaGaO3 (pc001, pc111)
% LSAT (pc001, pc111)
% MgO (c001,c111)
% Nb:SrTiO3 (c001, c111) 0.5%wt (Nb05STO)
% NdAlO3 (pc001, pc111)
% NdGaO3 (pc001, pc111)
% PMN71-PT29 (t001)
% Si(c001)
% SrLaAlO4 (t001)
% SrLaGaO4 (t001)
% SrTiO3 (c001, c111)
% TbScO3 (o001)
% TiO2 (t001)
% YAlO3 (pc001, pc111)
% YSZ (c001) 
% ZnO (hexagonal h0001, h11bar20, h10bar10) 
% Needs to run InteractiveXRDFit.m first
%*********************************

tic;
global plotdata;

%% Part to be changed to generate substrate
Substrate.Type='YSZc001';  %Choose between DyScO3, Gd3Ga5O12, GdScO3, KTaO3, LaAlO3, LaGaO3, LSAT, Nb05STOc001, Nb05STOc111, NdAlO3, NdGaO3, MgO, PMN71PT29, Si, SrLaAlO4, SrLaGaO4, SrTiO3, TbScO3, TiO2, YAlO3, YSZc001, ZnO (h0001, h10bar10, h11bar20, pc001, pc111, o001, t001, ...)

%% Main part of the program generating substrates
fO=AtomicScatteringFactor('O',plotdata.QQ);

zperovskite001(1)=0;    %A
zperovskite001(2)=0.5;  %B
zperovskite001(3)=0;    %O1
zperovskite001(4)=0.5;	%O2
zperovskite001(5)=0.5;	%O3

zperovskite111(1)=0;    %A
zperovskite111(2)=0.5;  %B
zperovskite111(3)=0;	%O1
zperovskite111(4)=0;	%O2
zperovskite111(5)=0;	%O3
            
switch Substrate.Type
    
    case 'DyScO3o001'
        % 4xDyScO3
        % rotated in-plane so that one unit-cell [ab] of DyScO3 corresponds to 2 unit-cell surface of traditional perovskite
        % orthorhombic a+b-b- with a,b in-plane, c out-of-plane
        % Pearson's Crystal Data
        a=5.443; b=5.717; c=7.901;  %alpha=beta=gamma=90 - b out-of-plane, a,c in-plane
        d=c;                        %out-of-plane lattice parameter
        V=a*b*c;                    %volume of the unit cell
        zDyScO3o001(1)=0;                     % 2 Sc
        zDyScO3o001(2)=0.0683869690058407;    % 2 O
        zDyScO3o001(3)=0.250000031295978;     % 2 Dy - 2 O
        zDyScO3o001(4)=0.431612968402204;     % 2 O
        zDyScO3o001(5)=0.500000062591955;     % 2 Sc
        zDyScO3o001(6)=0.568387031597796;     % 2 O
        zDyScO3o001(7)=0.750000093887933;     % 2 Dy - 2 O
        zDyScO3o001(8)=0.93161315617807;      % 2 O
        fDy=AtomicScatteringFactor('Dy',plotdata.QQ);
        fSc=AtomicScatteringFactor('Sc',plotdata.QQ);
        Fsubstrate=2*fSc.*exp(1i*(plotdata.Q.*zDyScO3o001(1).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zDyScO3o001(2).*d))...
        +2*(fDy+fO).*exp(1i*(plotdata.Q.*zDyScO3o001(3).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zDyScO3o001(4).*d))...
        +2*fSc.*exp(1i*(plotdata.Q.*zDyScO3o001(5).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zDyScO3o001(6).*d))...
        +2*(fDy+fO).*exp(1i*(plotdata.Q.*zDyScO3o001(7).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zDyScO3o001(8).*d));
        
    case 'DyScO3pc001' %orthorhombic (Pearson?s Crystal Data)
        a=5.443; b=5.717; c=7.901;
        h=1; k=1; l=0;
        d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2)); %(110)o
        V=a*b*c/4; %volume of the unit cell
        orthorhombic.a=a; orthorhombic.b=b; orthorhombic.c=c;
        orthorhombic.h=h; orthorhombic.k=k; orthorhombic.l=l;
        name='DyScO3o110';  % name used to store the values
        assignin('base',name,orthorhombic);
        save('Substrates.mat',name,'-append');        
        fAsubstrate=AtomicScatteringFactor('Dy',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Sc',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
        
    case 'DyScO3pc111' %orthorhombic (Pearson?s Crystal Data)
        a=5.443; b=5.717; c=7.901;
        h=1; k=0; l=1;
        d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2)); %(101)o
        d=d/2;
        V=a*b*c/4; %volume of the unit cell
        orthorhombic.a=a; orthorhombic.b=b; orthorhombic.c=c;
        orthorhombic.h=h; orthorhombic.k=k; orthorhombic.l=l;
        name='DyScO3o101';  % name used to store the values
        assignin('base',name,orthorhombic);
        save('Substrates.mat',name,'-append');  
        fAsubstrate=AtomicScatteringFactor('Dy',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Sc',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 

    case 'Gd3Ga5O12c001'
        % Gadolinium Gallium garnet (GGG)
        % cubic with c= 12.51341A from Materials Projet mp-2921
        % from Landolt-Bornstein: 1.2376nm
        % from Saint Gobain and Zanjani et al, AIP Advances 9, 035024 (2019): 1.2376nm (this is the value we use here)
        % one formula unit contains 20 atoms
        % conventional unit cell contains 160 atoms = 8 formula units
        % "new" unit cell contains 40 atoms = conventional unit cell/4
        a=12.383; b=a; c=a;         %a=12.376
        d=c/4;                      %out-of-plane lattice parameter of the "new" unit cell
        V=a*b*c/4;                  %volume of the "new" unit cell
        zGd3Ga5O12(1)=0*4;          % 4 Ga + 2 Gd
        zGd3Ga5O12(2)=0.027942*4;   % 4 O
        zGd3Ga5O12(3)=0.054232*4;   % 4 O
        zGd3Ga5O12(4)=0.100164*4;   % 4 O
        zGd3Ga5O12(5)=0.125000*4;   % 2 Ga + 2 Gd
        zGd3Ga5O12(6)=0.149836*4;   % 4 O
        zGd3Ga5O12(7)=0.195768*4;   % 4 O
        zGd3Ga5O12(8)=0.222058*4;   % 4 O
        zGd3Ga5O12(9)=0.250000*4;   % 4 Ga + 2 Gd
        fGd=AtomicScatteringFactor('Gd',plotdata.QQ);
        fGa=AtomicScatteringFactor('Ga',plotdata.QQ);
        Fsubstrate=(4*fGa+2*fGd).*exp(1i*(plotdata.Q.*zGd3Ga5O12(1).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(2).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(3).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(4).*d))...
        +(2*fGa+2*fGd).*exp(1i*(plotdata.Q.*zGd3Ga5O12(5).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(6).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(7).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(8).*d))...
        +(4*fGa+2*fGd).*exp(1i*(plotdata.Q.*zGd3Ga5O12(9).*d));
    
    case 'Gd3Ga5O12c111'
        % Gadolinium Gallium garnet (GGG)
        % cubic with c= 12.51341A from Materials Projet mp-2921
        % from Landolt-Bornstein: 1.2376nm
        % from Saint Gobain and Zanjani et al, AIP Advances 9, 035024 (2019): 1.2376nm (this is the value we use here)
        % one formula unit contains 20 atoms
        % conventional 111 unit cell contains 160 * 6 atoms = 8*6 formula units
        % "new" unit cell contains 40 atoms = conventional 111 unit cell/(4*6)
        a=12.383; b=a; c=a;         % a=12.376 %conventional 001 unit cell = 160 atoms
        d=c/4/sqrt(3);              %out-of-plane lattice parameter of the "new" unit cell
        V=a*b*c/4;                  %volume of the "new" unit cell
        zGd3Ga5O12(1)=0*24;          % 2 Ga
        zGd3Ga5O12(2)=0.0029983310*24;% 3 O
        zGd3Ga5O12(3)=0.0112770123*24;% 3 O
        zGd3Ga5O12(4)=0.0123123419*24;% 3 O
        zGd3Ga5O12(5)=0.0205910099*24;% 3 O
        zGd3Ga5O12(6)=0.0208333233*24;% 6 Ga + 6 Gd
        zGd3Ga5O12(7)=0.0210756768*24;% 3 O
        zGd3Ga5O12(8)=0.0293543314*24;% 3 O
        zGd3Ga5O12(9)=0.0303896610*24;% 3 0
        zGd3Ga5O12(10)=0.0386683290*24;% 3 0
        zGd3Ga5O12(11)=0.0416666600*24;% 2 Ga
        fGd=AtomicScatteringFactor('Gd',plotdata.QQ);
        fGa=AtomicScatteringFactor('Ga',plotdata.QQ);
        Fsubstrate=2*fGa.*exp(1i*(plotdata.Q.*zGd3Ga5O12(1).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(2).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(3).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(4).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(5).*d))...
        +(6*fGa+6*fGd).*exp(1i*(plotdata.Q.*zGd3Ga5O12(6).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(7).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(8).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(9).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zGd3Ga5O12(10).*d))...
        +2*fGa.*exp(1i*(plotdata.Q.*zGd3Ga5O12(11).*d));
    
    case 'GdScO3pc001'
        d=3.9636;
        V=d^3; %volume of the unit cell
        fAsubstrate=AtomicScatteringFactor('Gd',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Sc',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F

    case 'GdScO3pc111'
        d=3.9636/sqrt(3);
        V=3.9636^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Gd',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Sc',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
    
    case 'KTaO3c001'
        d=3.989;
        V=d^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('K',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ta',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
    
    case 'KTaO3c111'
        d=3.989/sqrt(3);
        V=3.989^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('K',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ta',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        
    case 'LaAlO3pc001'
        d=3.789;
        V=d^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('La',plotdata.QQ); 
        fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F
        
    case 'LaAlO3pc111'
        d=3.789/sqrt(3);
        V=3.789^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('La',plotdata.QQ); 
        fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
     
   case 'LaGaO3o110' %orthorhombic Ref:Marti-Journal of Physics: Condensed Matter-1994 4xLaGaO3
        a=5.5269; b=5.4943; c=7.7774;
        h=1; k=1; l=0;
        d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2)); %(110)o equivalent to (001)pc
        V=a*b*c/4; %unit cell volume
        orthorhombic.a=a; orthorhombic.b=b; orthorhombic.c=c;
        orthorhombic.h=h; orthorhombic.k=k; orthorhombic.l=l;
        name='LaGaO3o110';  % name used to store the values
        assignin('base',name,orthorhombic);
        save('Substrates.mat',name,'-append');        
        fAsubstrate=AtomicScatteringFactor('La',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ga',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
        
    case 'LaGaO3o101' %orthorhombic Ref:Marti-Journal of Physics: Condensed Matter-1994 4xLaGaO3
        a=5.5269; b=5.4943; c=7.7774;
        h=1; k=0; l=1;
        d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2)); %(101)o equivalent to (111)pc
        d=d/2;
        V=a*b*c/4; %unit cell volume
        orthorhombic.a=a; orthorhombic.b=b; orthorhombic.c=c;
        orthorhombic.h=h; orthorhombic.k=k; orthorhombic.l=l;
        name='LaGaO3o101';  % name used to store the values
        assignin('base',name,orthorhombic);
        save('Substrates.mat',name,'-append');  
        fAsubstrate=AtomicScatteringFactor('La',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ga',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F
    case 'LaGaO3o001' %orthorhombic Ref:Marti-Journal of Physics: Condensed Matter-1994 4xLaGaO3
        a=5.5269; b=5.4943; c=7.7774;
        h=0; k=0; l=1;
        d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2)); %(001)o equivalent to (100)pc
        d=d/2;
        V=a*b*c/4; %unit cell volume
        orthorhombic.a=a; orthorhombic.b=b; orthorhombic.c=c;
        orthorhombic.h=h; orthorhombic.k=k; orthorhombic.l=l;
        name='LaGaO3o001';  % name used to store the values
        assignin('base',name,orthorhombic);
        save('Substrates.mat',name,'-append');  
        fAsubstrate=AtomicScatteringFactor('La',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ga',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F
        
    case 'LSATpc001'
        d=3.868;
        V=d^3; %unit cell volume
        fAsubstrate=0.3*AtomicScatteringFactor('La',plotdata.QQ)+0.7*AtomicScatteringFactor('Sr',plotdata.QQ);
        fBsubstrate=0.65*AtomicScatteringFactor('Al',plotdata.QQ)+0.35*AtomicScatteringFactor('Ta',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F
        
    case 'LSATpc111'
        d=3.868/sqrt(3);
        V=3.868^3; %unit cell volume
        fAsubstrate=0.3*AtomicScatteringFactor('La',plotdata.QQ)+0.7*AtomicScatteringFactor('Sr',plotdata.QQ);
        fBsubstrate=0.65*AtomicScatteringFactor('Al',plotdata.QQ)+0.35*AtomicScatteringFactor('Ta',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F
        
    case 'MgOc001'
        % cubic a=4.212 A
        % (Crystec Datasheets http://www.crystec.de/daten/mgo.pdf)
        d=4.212;
        V=d^3; %unit cell volume
        zMgOc001(1)=0;                      % Mg-O
        zMgOc001(2)=0.5;                    % 2Mg-2O
        zMgOc001(3)=1;                      % MgO
        fMg=AtomicScatteringFactor('Mg',plotdata.QQ);
        Fsubstrate=(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(1).*d))...
        +2*(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(2).*d))...
        +(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(3).*d));  %Structure factor F 
    
    case 'MgOc111'
        d=4.212/sqrt(3);
        V=4.212^3; %unit cell volume
        zMgOc111(1)=0;                    % 2Mg
        zMgOc111(2)=0.5;                  % 4O
        zMgOc111(3)=1;                    % 2Mg
        fMg=AtomicScatteringFactor('Mg',plotdata.QQ);    
        Fsubstrate=2*fMg.*exp(1i*(plotdata.Q.*zMgOc111(1).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zMgOc111(2).*d))...
        +2*fMg.*exp(1i*(plotdata.Q.*zMgOc111(3).*d));  %Structure factor F 
    
            
    case 'Nb05STOc001' %0.5%wt Nb:STO corresponds to Sr (Ti 1-x, Nb x) O3 with x = 0.99%, i.e 0.0099Nb and 0.9901Ti
        d=3.9068;
        V=d^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Sr',plotdata.QQ);
        fBsubstrate=0.0099*AtomicScatteringFactor('Nb',plotdata.QQ)+(1-0.0099)*AtomicScatteringFactor('Ti',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F
        
    case 'Nb05STOc111' %0.5%wt Nb:STO corresponds to Sr (Ti 1-x, Nb x) O3 with x = 0.99%, i.e 0.0099Nb and 0.9901Ti
        d=3.9068/sqrt(3);
        V=3.9068^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Sr',plotdata.QQ);
        fBsubstrate=0.0099*AtomicScatteringFactor('Nb',plotdata.QQ)+(1-0.0099)*AtomicScatteringFactor('Ti',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F
        
    case 'NdAlO3pc001'
        d=3.74;
        V=d^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Nd',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
        
    case 'NdAlO3pc111' 
        d=3.74/sqrt(3);
        V=3.74^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Nd',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F
        
    case 'NdGaO3pc001' %orthorhombic (10.1103/PhysRevB.80.064408) 4xNdGaO3
        a=5.426; b=5.496; c=7.706;
        h=1; k=1; l=0;
        d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2)); %(110)o
        V=a*b*c/4; %unit cell volume
        orthorhombic.a=a; orthorhombic.b=b; orthorhombic.c=c;
        orthorhombic.h=h; orthorhombic.k=k; orthorhombic.l=l;
        name='NdGaO3o110';  % name used to store the values
        assignin('base',name,orthorhombic);
        save('Substrates.mat',name,'-append');        
        fAsubstrate=AtomicScatteringFactor('Nd',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ga',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
        
    case 'NdGaO3pc111' %orthorhombic (10.1103/PhysRevB.80.064408) 4xNdGaO3
        a=5.426; b=5.496; c=7.706;
        h=1; k=0; l=1;
        d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2)); %(101)o
        d=d/2;
        V=a*b*c/4; %unit cell volume
        orthorhombic.a=a; orthorhombic.b=b; orthorhombic.c=c;
        orthorhombic.h=h; orthorhombic.k=k; orthorhombic.l=l;
        name='NdGaO3o101';  % name used to store the values
        assignin('base',name,orthorhombic);
        save('Substrates.mat',name,'-append');  
        fAsubstrate=AtomicScatteringFactor('Nd',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ga',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F
        
    case 'PMN71PT29t001'
        %PMN-PT tetragonal 001 a=b=4.0166 c=4.02870 (PMN71PT29_P4mm_155863.cif)
        a=4.0166; b=a; c=4.02870;
        d=c;
        V=a*b*c; %unit cell volume
        zPMN71PT29(1)=0;         % 1 Pb
        zPMN71PT29(2)=0.52100;   % 0.29 Ti - 0.473 Nb - 0.237 Mg
        zPMN71PT29(3)=0.03800;   % 1 O
        zPMN71PT29(4)=0.54950;   % 2 O
        fPb=AtomicScatteringFactor('Pb',plotdata.QQ);
        fTiNbMg=0.29*AtomicScatteringFactor('Ti',plotdata.QQ)+0.473*AtomicScatteringFactor('Nb',plotdata.QQ)+0.237*AtomicScatteringFactor('Mg',plotdata.QQ);     
        Fsubstrate=(fPb).*exp(1i*(plotdata.Q.*zPMN71PT29(1).*d))...
            +(fTiNbMg).*exp(1i*(plotdata.Q.*zPMN71PT29(2).*d))...
            +fO.*exp(1i*(plotdata.Q.*zPMN71PT29(3).*d))...
            +2*fO.*exp(1i*(plotdata.Q.*zPMN71PT29(4).*d));
        
    case 'Sic001' %8Si
        %rotated in-plane so one unit-cell surface of silicon correspond to 2 unit-cell surface of traditional perovskite
        d=5.4307;
        V=d^3; %unit cell volume
        fSi=AtomicScatteringFactor('Si',plotdata.QQ);
        zSi(1)=0;                     %2Si
        zSi(2)=0.25;                  %2Si
        zSi(3)=0.5;                   %2Si
        zSi(4)=0.75;                  %2Si
        Fsubstrate=2*fSi.*exp(1i*(plotdata.Q.*zSi(1).*d))...
            +2*fSi.*exp(1i*(plotdata.Q.*zSi(2).*d))...
            +2*fSi.*exp(1i*(plotdata.Q.*zSi(3).*d))...
            +2*fSi.*exp(1i*(plotdata.Q.*zSi(4).*d));
        
    case 'SrLaAlO4t001'
        %1531895 CIF a=3.7544 b=3.7544 c=12.6494 alpha=beta=gamma=90
        a=3.7544; b=a; c=12.6494;
        d=12.6377;
        V=a*b*c;
        zSrLaAlO4(1)=0;                     % 1 Al - 2 O
        zSrLaAlO4(2)=0.141099973121255;     % 0.5 La - 0.5 Sr
        zSrLaAlO4(3)=0.162699969959049;     % 1 O
        zSrLaAlO4(4)=0.337300030040951;     % 1 O
        zSrLaAlO4(5)=0.358900026878745;     % 0.5 La - 0.5 Sr
        zSrLaAlO4(6)=0.5;                   % 1 Al - 2 O
        zSrLaAlO4(7)=0.641099973121255;     % 0.5 La - 0.5 Sr
        zSrLaAlO4(8)=0.662699969959049;     % 1 O
        zSrLaAlO4(9)=0.837300030040951;     % 1 O
        zSrLaAlO4(10)=0.858900026878745;    % 0.5 La - 0.5 Sr
        fLa=AtomicScatteringFactor('La',plotdata.QQ);
        fSr=AtomicScatteringFactor('Sr',plotdata.QQ);
        fAl=AtomicScatteringFactor('Al',plotdata.QQ);
        Fsubstrate=(fAl+2*fO).*exp(1i*(plotdata.Q.*zSrLaAlO4(1).*d))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zSrLaAlO4(2).*d))...
            +fO.*exp(1i*(plotdata.Q.*zSrLaAlO4(3).*d))...
            +fO.*exp(1i*(plotdata.Q.*zSrLaAlO4(4).*d))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zSrLaAlO4(5).*d))...
            +(fAl+2*fO).*exp(1i*(plotdata.Q.*zSrLaAlO4(6).*d))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zSrLaAlO4(7).*d))...
            +fO.*exp(1i*(plotdata.Q.*zSrLaAlO4(8).*d))...
            +fO.*exp(1i*(plotdata.Q.*zSrLaAlO4(9).*d))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zSrLaAlO4(10).*d));
        
    case 'SrLaGaO4t001'
        %SrLaGaO4 tetragonal 001 a=b=3.84 c=12.68 (Crystec datasheets)
        a=3.84; b=a; c=12.68;
        d=c;
        V=a*b*c; %unit cell volume
        zSrLaGaO4(1)=0;         % 1 Ga - 2 O
        zSrLaGaO4(2)=0.1412;    % 0.5 La - 0.5 Sr
        zSrLaGaO4(3)=0.168;     % 1 O
        zSrLaGaO4(4)=0.332;     % 1 O
        zSrLaGaO4(5)=0.3588;    % 0.5 La - 0.5 Sr
        zSrLaGaO4(6)=0.5;       % Ga - 2 O
        zSrLaGaO4(7)=0.6412;    % 0.5 La - 0.5 Sr
        zSrLaGaO4(8)=0.668;     % 1 O
        zSrLaGaO4(9)=0.832;     % 1 O
        zSrLaGaO4(10)=0.8588;   % 0.5 La - 0.5 Sr
        fLa=AtomicScatteringFactor('La',plotdata.QQ);
        fSr=AtomicScatteringFactor('Sr',plotdata.QQ);
        fGa=AtomicScatteringFactor('Ga',plotdata.QQ);
        Fsubstrate=(fGa+2*fO).*exp(1i*(plotdata.Q.*zSrLaGaO4(1).*d))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zSrLaGaO4(2).*d))...
            +fO.*exp(1i*(plotdata.Q.*zSrLaGaO4(3).*d))...
            +fO.*exp(1i*(plotdata.Q.*zSrLaGaO4(4).*d))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zSrLaGaO4(5).*d))...
            +(fGa+2*fO).*exp(1i*(plotdata.Q.*zSrLaGaO4(6).*d))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zSrLaGaO4(7).*d))...
            +fO.*exp(1i*(plotdata.Q.*zSrLaGaO4(8).*d))...
            +fO.*exp(1i*(plotdata.Q.*zSrLaGaO4(9).*d))...
            +(0.5*fLa+0.5*fSr).*exp(1i*(plotdata.Q.*zSrLaGaO4(10).*d));
        
    case 'SrTiO3c001'
        d=3.905;
        V=d^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Sr',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ti',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F
        
    case 'SrTiO3c111'
        d=3.905/sqrt(3);
        V=3.905^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Sr',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ti',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F
        
    case 'TbScO3o001'
        %rotated in-plane so that one unit-cell [ab] of TbScO3 corresponds to 2 unit-cell surface of traditional perovskite
        %orthorhombic a=5.4654A b=5.7292A c=7.917A with a,b in-plane, c out-of-plane
        a=5.4654; b=5.7292; c=7.917;
        d=c;
        V=a*b*c; %unit cell volume
        zTbScO3o(1)=0;                     % 2 Sc
        zTbScO3o(2)=0.0564000252620942;    % 2 O
        zTbScO3o(3)=0.25;                  % 2 Tb - 2 O
        zTbScO3o(4)=0.443599974737906;     % 2 O
        zTbScO3o(5)=0.5;                   % 2 Sc
        zTbScO3o(6)=0.556400025262094;     % 2 O
        zTbScO3o(7)=0.75;                  % 2 Tb - 2 O
        zTbScO3o(8)=0.943599974737906;     % 2 O
        fTb=AtomicScatteringFactor('Tb',plotdata.QQ);
        fSc=AtomicScatteringFactor('Sc',plotdata.QQ);
        Fsubstrate=2*fSc.*exp(1i*(plotdata.Q.*zTbScO3o(1).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zTbScO3o(2).*d))...
        +2*(fTb+fO).*exp(1i*(plotdata.Q.*zTbScO3o(3).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zTbScO3o(4).*d))...
        +2*fSc.*exp(1i*(plotdata.Q.*zTbScO3o(5).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zTbScO3o(6).*d))...
        +2*(fTb+fO).*exp(1i*(plotdata.Q.*zTbScO3o(7).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zTbScO3o(8).*d));
    
    case 'TiO2t001'
        %TiO2-rutile tetragonal 001 a=b=4.606777A c=2.991757A (Materialsproject.org mp-2657)
        a=4.606777; b=a; c=2.991757;
        d=c;
        V=a*b*c; %unit cell volume
        zTiO2(1)=0;     % 0.5 Ti - 1 O
        zTiO2(2)=0.5;   % 1 Ti - 2 O
        zTiO2(3)=1;     % 0.5 Ti - 1 O
        fTi=AtomicScatteringFactor('Ti',plotdata.QQ);
        Fsubstrate=(0.5*fTi+fO).*exp(1i*(plotdata.Q.*zTiO2(1).*d))...
            +(fTi+2*fO).*exp(1i*(plotdata.Q.*zTiO2(2).*d))...
            +(0.5*fTi+fO).*exp(1i*(plotdata.Q.*zTiO2(3).*d));
        
    case 'YAlO3pc001' 
        d=3.71;
        V=d^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Y',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);
        Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            
    case 'YSZc001' %yttria stabilized zirconia, ZrO2 + 9.5mol% correspond to nominal formula Zr_0.826Y_0.174O_1.913
        d=5.1453; %as measured
        V=d^3; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Zr',plotdata.QQ)*0.826+AtomicScatteringFactor('Y',plotdata.QQ)*0.174;
        fBsubstrate=AtomicScatteringFactor('O',plotdata.QQ)*1.913;        
        zYSZc001(1)=0;                     % 2 A (2xZr_0.826Y_0.174)
        zYSZc001(2)=0.25;                  % 2 B (2xO_1.913)
        zYSZc001(3)=0.5;                   % 2 A (2xZr_0.826Y_0.174)
        zYSZc001(4)=0.75;                  % 2 B (2xO_1.913)
        Fsubstrate=fAsubstrate.*exp(1i*(plotdata.Q.*zYSZc001(1).*d))...
        +fBsubstrate.*exp(1i*(plotdata.Q.*zYSZc001(2).*d))...
        +fAsubstrate.*exp(1i*(plotdata.Q.*zYSZc001(3).*d))...
        +fBsubstrate.*exp(1i*(plotdata.Q.*zYSZc001(4).*d));

    case 'ZnOh0001' %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
                    % (0001): z along c-axis
                    % inplane surface=sqrt(3)/2*a^2
                    % V(unit cell)=sqrt(3)*a^2*c/2
        a=3.2494; b=a; c=5.20380;
        d=c;
        V=sqrt(3)*a^2*c/2; %unit cell volume
        zZnOh0001(1)=0;     %Zn
        zZnOh0001(2)=0.5;   %Zn
        zZnOh0001(3)=0.3821;%O
        zZnOh0001(4)=0.8821;%O
        fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
        Fsubstrate=fZn.*exp(1i*(plotdata.Q.*zZnOh0001(1).*d))...
            +fZn.*exp(1i*(plotdata.Q.*zZnOh0001(2).*d))...
            +fO.*exp(1i*(plotdata.Q.*zZnOh0001(3).*d))...
            +fO.*exp(1i*(plotdata.Q.*zZnOh0001(4).*d)); 
        
    case 'ZnOh11bar20' %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
                    % (11-20): z along bissectrice a,b
                    % inplane surface S=sqrt(3)*a*c
                    % distance between planes = V(unit cell)/S = a/2
        a=3.2494; b=a; c=5.20380;
        d=a/2;
        V=sqrt(3)*a^2*c/2; %unit cell volume
        zZnOh11bar20(1)=0;     %2 Zn + 2 O
        fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
        Fsubstrate=(2*fZn+2*fO).*exp(1i*(plotdata.Q.*zZnOh11bar20(1).*d));
        
    case 'ZnOh10bar10' %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
                    % (10-10): z perpendicular to a
                    % inplane surface S=ac
                    % distance between planes = V(unit cell)/S = sqrt(3)a/2
        a=3.2494; b=a; c=5.20380;
        d=a*sqrt(3)/2;
        V=sqrt(3)*a^2*c/2; %unit cell volume
        zZnOh10bar10(1)=1/3;   % Zn + O
        zZnOh10bar10(2)=2/3;     % Zn + O
        fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
        Fsubstrate=(fZn+fO).*exp(1i*(plotdata.Q.*zZnOh10bar10(1).*d))+(fZn+fO).*exp(1i*(plotdata.Q.*zZnOh10bar10(2).*d));
        
    otherwise
        warning('Error in Substrate.Type in GenerateSubstrates.')
end;  

%% Calculations

expzmsubstrate=zeros(size(plotdata.Q));

NSubstrate=round(plotdata.Substrate.thickness/d);
h = waitbar(0.2,'Please wait...');
for zmsubstrate=1:NSubstrate,
   expzmsubstrate=expzmsubstrate+exp(1i*plotdata.Q.*zmsubstrate.*d).*exp(-(plotdata.TotalThickness-zmsubstrate.*d)./plotdata.mu);
end;
close(h) 

Fsubstrate=Fsubstrate*d/V; %renormalized by the in-plane surface of the unit cell

%Scattering amplitude g:
Substrate.g=Fsubstrate.*expzmsubstrate;

toc;

gname=['gsubstrate',num2str(Substrate.Type)];  % name used to store the calculated intensity
dname=['dsubstrate',num2str(Substrate.Type)];  % name used to store the out-of-plane lattice parameter
Vname=['Vsubstrate',num2str(Substrate.Type)];  % name used to store the unit cell volume

assignin('base',gname,Substrate.g);
assignin('base',dname,d);
assignin('base',Vname,V);

save('Substrates.mat',gname,'-append');
save('Substrates.mat',dname,'-append');
save('Substrates.mat',Vname,'-append');