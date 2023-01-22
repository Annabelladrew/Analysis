% Created by Celine Lichtensteiger
% Calculate new substrate intensities (user defined lattice parameter and penetration depth) among:
% Al2O3 (h0001))
% DyScO3 (o001,o110,o101)
% beta-Ga2O3 (m010) monoclinic (010) i.e. out-of-plane b
% Gd3Ga5O12 (garnet c001, c111)
% GdScO3 (o001,o110,o101)
% KTaO3 (c001, c111)
% LaAlO3 (pc001, pc111)
% LaGaO3 (pc001, pc111)
% LSAT (pc001, pc111)
% MgO (c001,c111)
% Nb:SrTiO3 (c001, c111) 0.5%wt (Nb05STO)
% NdAlO3 (pc001, pc111)
% NdGaO3 (pc001, pc111)
% NdScO3 (o001,pc001,pc111)
% PMN71-PT29 (t001)
% Si(c001,c111)
% SrLaAlO4 (t001)
% SrLaGaO4 (t001)
% SrTiO3 (c001, c111)
% TbScO3 (o001)
% TiO2 (t001)
% YAlO3 (pc001, pc111)
% YSZ (c001, c111) 
% ZnO (hexagonal h0001, h11bar20, h10bar10) 
%*********************************

%tic;
global plotdata;
load('Substrates.mat');
SubstrateTable=struct2table(SubstrateData);
condition=find(strcmp(SubstrateTable.Name,char(plotdata.Substrate.Type)));
a=SubstrateData(condition).a; b=SubstrateData(condition).b; c=SubstrateData(condition).c;
alpha=SubstrateData(condition).alpha*pi/180; beta=SubstrateData(condition).beta*pi/180; gamma=SubstrateData(condition).gamma*pi/180;
z=[]; f=[];
fA=zeros(size(plotdata.QQ));
fB=zeros(size(plotdata.QQ));
Fsubstrate=fA.*exp(1i*(plotdata.Q.*0));

%% Part to be changed to generate substrate
%Choose between Al2O3, Gd3Ga5O12, GdScO3, KTaO3, LaAlO3, LaGaO3, LSAT, MgO, Nb05STOc001, Nb05STOc111, NdAlO3, NdGaO3, PMN71PT29, Si, SrLaAlO4, SrLaGaO4, SrTiO3, TbScO3, TiO2, YAlO3, YSZc001, ZnO (h0001, h10bar10, h11bar20, pc001, pc111, o001, t001, ...)
%% Main part of the program generating substrates

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

fO=AtomicScatteringFactor('O',plotdata.QQ);

switch char(plotdata.Substrate.Type)
    
    case 'Al2O3'
        %cf saphir with orientation c out-of-plane
        %alpha= beta=90.0000  %gamma=120.0000*pi/180;
        %30 atoms/unit cell: 12 Al + 18 O
        %a=4.76060; b=a; c=12.99400;
        switch plotdata.Substrate.Orientation
            case '(0001)h'
                d=c;%out-of-plane lattice parameter
                V=a*b*c*sin(gamma);%V=a*b*d*sin(gamma);
                zAl2O3h0001(1)=0.0188367708173003;  % 1 Al
                zAl2O3h0001(2)=0.0833333846390642;  % 3 O
                zAl2O3h0001(3)=0.147829998460828;   % 1 Al
                zAl2O3h0001(4)=0.185503232261044;   % 1 Al
                zAl2O3h0001(5)=0.25;                % 3 O
                zAl2O3h0001(6)=0.31449669078036;    % 1 Al
                zAl2O3h0001(7)=0.352170001539172;   % 1 Al
                zAl2O3h0001(8)=0.416666692319532;   % 3 O
                zAl2O3h0001(9)=0.481163306141296;   % 1 Al
                zAl2O3h0001(10)=0.518836693858704;  % 1 Al
                zAl2O3h0001(11)=0.583333384639064;  % 3 O
                zAl2O3h0001(12)=0.647830075419424;  % 1 Al
                zAl2O3h0001(13)=0.685503386178236;  % 1 Al
                zAl2O3h0001(14)=0.750000076958596;  % 3 O
                zAl2O3h0001(15)=0.81449669078036;   % 1 Al
                zAl2O3h0001(16)=0.852170001539172;  % 1 Al
                zAl2O3h0001(17)=0.916666692319532;  % 3 O
                zAl2O3h0001(18)=0.981163383099892;  % 1 Al
                fAl=AtomicScatteringFactor('Al',plotdata.QQ);
                Fsubstrate=fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(1).*d))...
                +3*fO.*exp(1i*(plotdata.Q.*zAl2O3h0001(2).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(3).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(4).*d))...
                +3*fO.*exp(1i*(plotdata.Q.*zAl2O3h0001(5).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(6).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(7).*d))...
                +3*fO.*exp(1i*(plotdata.Q.*zAl2O3h0001(8).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(9).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(10).*d))...
                +3*fO.*exp(1i*(plotdata.Q.*zAl2O3h0001(11).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(12).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(13).*d))...
                +3*fO.*exp(1i*(plotdata.Q.*zAl2O3h0001(14).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(15).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(16).*d))...
                +3*fO.*exp(1i*(plotdata.Q.*zAl2O3h0001(17).*d))...
                +fAl.*exp(1i*(plotdata.Q.*zAl2O3h0001(18).*d));
        end
                
    case 'DyScO3'%orthorhombic
        % 4xDyScO3
        % rotated in-plane so that one unit-cell [ab] of DyScO3 corresponds to 2 unit-cell surface of traditional perovskite
        % orthorhombic a+b-b- with a,b in-plane, c out-of-plane
        % Pearson's Crystal Data
        %a=5.443; b=5.717; c=7.901;  %alpha=beta=gamma=90 - b out-of-plane, a,c in-plane
        %d=c;                        %out-of-plane lattice parameter
        %V=a*b*d;                    %volume of the unit cell
        fAsubstrate=AtomicScatteringFactor('Dy',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Sc',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(001)o'
                d=c;%out-of-plane lattice parameter
                V=a*b*c;
                zDyScO3o001(1)=0;                     % 2 Sc
                zDyScO3o001(2)=0.0683869690058407;    % 2 O
                zDyScO3o001(3)=0.250000031295978;     % 2 Dy - 2 O
                zDyScO3o001(4)=0.431612968402204;     % 2 O
                zDyScO3o001(5)=0.500000062591955;     % 2 Sc
                zDyScO3o001(6)=0.568387031597796;     % 2 O
                zDyScO3o001(7)=0.750000093887933;     % 2 Dy - 2 O
                zDyScO3o001(8)=0.93161315617807;      % 2 O
                Fsubstrate=2*fSc.*exp(1i*(plotdata.Q.*zDyScO3o001(1).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zDyScO3o001(2).*d))...
                +2*(fDy+fO).*exp(1i*(plotdata.Q.*zDyScO3o001(3).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zDyScO3o001(4).*d))...
                +2*fSc.*exp(1i*(plotdata.Q.*zDyScO3o001(5).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zDyScO3o001(6).*d))...
                +2*(fDy+fO).*exp(1i*(plotdata.Q.*zDyScO3o001(7).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zDyScO3o001(8).*d));
            case '(110)o'%(001)pc
                h=1; k=1; l=0;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                V=a*b*c/4; %volume of the unit cell
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(101)o'%(111)pc
                h=1; k=0; l=1;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                d=d/2;
                V=a*b*c/4; %volume of the unit cell
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end
        
    case 'beta-Ga2O3'%monoclinic
        % 8xGa + 11xO +8xC
        % growth on (010) planes: a,c in-plane, b out-of-plane
        % a=12.214; b=3.0371; c=5.7981; alpha=gamma=90 beta=103.83°
        % d=b;                          %out-of-plane lattice parameter
        % S=a*c*sin(beta)               %in-plane surface
        % V=a*b*c*sin(beta);            %volume of the unit cell
        V=a*b*c*sin(beta*pi/180);
        fGa=AtomicScatteringFactor('Ga',plotdata.QQ);
        fC=AtomicScatteringFactor('C',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(010)m'
                d=b;%out-of-plane lattice parameter
                zbetaGa2O3(1)=0;       % 3 O + 2 Ga
                zbetaGa2O3(2)=0.141618;% 2 C
                zbetaGa2O3(3)=0.358382;% 2 C
                zbetaGa2O3(4)=0.5;     % 5 O + 4 Ga
                zbetaGa2O3(5)=0.641618;% 2 C
                zbetaGa2O3(6)=0.858382;% 2 C
                zbetaGa2O3(7)=1;       % 3 O + 2 Ga
                Fsubstrate=(3*fO+2*fGa).*exp(1i*(plotdata.Q.*zbetaGa2O3(1).*d))...
                +2*fC.*exp(1i*(plotdata.Q.*zbetaGa2O3(2).*d))...
                +2*fC.*exp(1i*(plotdata.Q.*zbetaGa2O3(3).*d))...
                +(5*fO+4*fGa).*exp(1i*(plotdata.Q.*zbetaGa2O3(4).*d))...
                +2*fC.*exp(1i*(plotdata.Q.*zbetaGa2O3(5).*d))...
                +2*fC.*exp(1i*(plotdata.Q.*zbetaGa2O3(6).*d))...
                +(3*fO+2*fGa).*exp(1i*(plotdata.Q.*zbetaGa2O3(7).*d));
        end
        
    case 'Gd3Ga5O12'%Cubic
        % Gadolinium Gallium garnet (GGG)
        % cubic with c= 12.51341A from Materials Projet mp-2921
        % from Landolt-Bornstein: 1.2376nm from Saint Gobain and Zanjani et al, AIP Advances 9, 035024 (2019)
        % c=12.383 (default)
        % : 1.2376nm (this is the value we use here)
        % one formula unit contains 20 atoms
        % conventional unit cell contains 160 atoms = 8 formula units
        % "new" unit cell contains 40 atoms = conventional unit cell/4
        switch plotdata.Substrate.Orientation
            case '(001)c'
                d=c/4;                    %out-of-plane lattice parameter of the "new" unit cell
                V=a*b*c;                  %volume of the "new" unit cell
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
            case '(111)c'
                d=c/4/sqrt(3);              %out-of-plane lattice parameter of the "new" unit cell
                V=a*b*c*sqrt(3);            %volume of the "new" unit cell
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
        end
        
     case 'GdScO3'%orthorhombic
        % 4xGdScO3
        % rotated in-plane so that one unit-cell [ab] of GdScO3 corresponds to 2 unit-cell surface of traditional perovskite
        % orthorhombic with a,b in-plane, c out-of-plane
        %a=5.480; b=5.746; c=7.932;  %alpha=beta=gamma=90 - c out-of-plane, a,b in-plane
        %d=c;                        %out-of-plane lattice parameter
        %V=a*b*d;                    %volume of the unit cell
        fAsubstrate=AtomicScatteringFactor('Gd',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Sc',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(001)o'
                d=c;%out-of-plane lattice parameter
                V=a*b*c;
                zGdScO3o001(1)=0;                   % 2 Gd
                zGdScO3o001(2)=0.065105;            % 2 O
                zGdScO3o001(3)=0.250000;            % 2 Gd - 2 O
                zGdScO3o001(4)=0.434895;            % 2 O
                zGdScO3o001(5)=0.500000;            % 2 Sc
                zGdScO3o001(6)=0.565105;            % 2 O
                zGdScO3o001(7)=0.750000;            % 2 Gd - 2 O
                zGdScO3o001(8)=0.934895;            % 2 O
                Fsubstrate=2*fBsubstrate.*exp(1i*(plotdata.Q.*zGdScO3o001(1).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zGdScO3o001(2).*d))...
                +2*(fAsubstrate+fO).*exp(1i*(plotdata.Q.*zGdScO3o001(3).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zGdScO3o001(4).*d))...
                +2*fBsubstrate.*exp(1i*(plotdata.Q.*zGdScO3o001(5).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zGdScO3o001(6).*d))...
                +2*(fAsubstrate+fO).*exp(1i*(plotdata.Q.*zGdScO3o001(7).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zGdScO3o001(8).*d));
            case '(110)o'%(001)pc
                h=1; k=1; l=0;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                V=a*b*c/4; %volume of the unit cell
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(101)o'%(111)pc
                h=1; k=0; l=1;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                d=d/2;
                V=a*b*c/4; %volume of the unit cell
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end    
    
    case 'KTaO3' %cubic 
        %d=3.989;
        fAsubstrate=AtomicScatteringFactor('K',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ta',plotdata.QQ);
        V=a*b*c; %unit cell volume
        switch plotdata.Substrate.Orientation
            case '(001)c'
                d=c;
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(111)c'                
                d=c/sqrt(3);
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end

    case 'LaAlO3' %PseudoCubic
        fAsubstrate=AtomicScatteringFactor('La',plotdata.QQ); 
        fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);
        V=a*b*c; %unit cell volume
        switch plotdata.Substrate.Orientation
            case '(001)pc'
                d=c;
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(111)pc'                
                d=c/sqrt(3);
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end
        
    case 'LaGaO3'%orthorhombic
        %Ref:Marti-Journal of Physics: Condensed Matter-1994 4xLaGaO3
        %a=5.5269; b=5.4943; c=7.7774;
        V=a*b*c/4; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('La',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ga',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(110)o' %equivalent to (001)pc
                h=1; k=1; l=0;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(101)o'
                h=1; k=0; l=1;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2)); %equivalent to (111)pc
                d=d/2;
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F
            case '(001)o' %CHECK
                h=0; k=0; l=1;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2)); %equivalent to (100)pc CHECK - should be zperosvkite001
                d=d/2;
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F %CHECK
        end;
        
    case 'LSAT'%PseudoCubic
        %d=3.868;        
        V=a*b*c; %unit cell volume
        fAsubstrate=0.3*AtomicScatteringFactor('La',plotdata.QQ)+0.7*AtomicScatteringFactor('Sr',plotdata.QQ);
        fBsubstrate=0.65*AtomicScatteringFactor('Al',plotdata.QQ)+0.35*AtomicScatteringFactor('Ta',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(001)pc'
                d=c;
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(111)pc'                
                d=c/sqrt(3);
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end   
        
    case 'MgO'%Cubic
        % cubic a=4.212 A
        % (Crystec Datasheets http://www.crystec.de/daten/mgo.pdf)
        V=a*b*c; %unit cell volume
        fMg=AtomicScatteringFactor('Mg',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(001)c'
                d=c;
                zMgOc001(1)=0;                      % Mg-O
                zMgOc001(2)=0.5;                    % 2Mg-2O
                zMgOc001(3)=1;                      % MgO
                Fsubstrate=(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(1).*d))...
                +2*(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(2).*d))...
                +(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(3).*d));  %Structure factor F             
            case '(111)c'                
                d=c/sqrt(3);
                zMgOc111(1)=0;                    % 2Mg
                zMgOc111(2)=0.5;                  % 4O
                zMgOc111(3)=1;                    % 2Mg
                Fsubstrate=2*fMg.*exp(1i*(plotdata.Q.*zMgOc111(1).*d))...
                +4*fO.*exp(1i*(plotdata.Q.*zMgOc111(2).*d))...
                +2*fMg.*exp(1i*(plotdata.Q.*zMgOc111(3).*d));  %Structure factor F         
        end
            
    case 'Nb:SrTiO30.5%wt'%Cubic
        %0.5%wt Nb:STO corresponds to Sr (Ti 1-x, Nb x) O3 with x = 0.99%, i.e 0.0099Nb and 0.9901Ti
        %d=3.9068;
        V=a*b*c; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Sr',plotdata.QQ);
        fBsubstrate=0.0099*AtomicScatteringFactor('Nb',plotdata.QQ)+(1-0.0099)*AtomicScatteringFactor('Ti',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(001)c'
                d=c;
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(111)c'                
                d=c/sqrt(3);
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end
        
    case 'NdAlO3'%PseudoCubic
        %d=3.74;
        V=a*b*c; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Nd',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);      
        switch plotdata.Substrate.Orientation
            case '(001)pc'
                d=c;
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(111)pc'                
                d=c/sqrt(3);
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end   
        
    case 'NdGaO3'%Orthorhombic
        %(10.1103/PhysRevB.80.064408) 4xNdGaO3
        %a=5.426; b=5.496; c=7.706;        
        V=a*b*c/4; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Nd',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ga',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(110)o'%equivalent to (001)pc
                h=1; k=1; l=0;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(101)o'%equivalent to (111)pc
                h=1; k=0; l=1;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2)); 
                d=d/2;
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F
        end

    case 'NdScO3'%orthorhombic
        % 4xNdScO3
        % rotated in-plane so that one unit-cell [ab] of NdScO3 corresponds to 2 unit-cell surface of traditional perovskite
        % orthorhombic a+b-b- with a,b in-plane, c out-of-plane
        %a=5.430; b=5.500; c=7.710;  %alpha=beta=gamma=90 - b out-of-plane, a,c in-plane
        %d=c;                        %out-of-plane lattice parameter
        %V=a*b*d;                    %volume of the unit cell
        fAsubstrate=AtomicScatteringFactor('Nd',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Sc',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(001)o'
                d=c;%out-of-plane lattice parameter
                V=a*b*c;
                zNdScO3o001(1)=0;                   % 2 Sc
                zNdScO3o001(2)=0.058368;            % 2 O
                zNdScO3o001(3)=0.250000;            % 2 Nd - 2 O
                zNdScO3o001(4)=0.441632;            % 2 O
                zNdScO3o001(5)=0.500000;            % 2 Sc
                zNdScO3o001(6)=0.558368;            % 2 O
                zNdScO3o001(7)=0.750000;            % 2 Nd - 2 O
                zNdScO3o001(8)=0.941632;            % 2 O
                Fsubstrate=2*fBsubstrate.*exp(1i*(plotdata.Q.*zNdScO3o001(1).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zNdScO3o001(2).*d))...
                +2*(fAsubstrate+fO).*exp(1i*(plotdata.Q.*zNdScO3o001(3).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zNdScO3o001(4).*d))...
                +2*fBsubstrate.*exp(1i*(plotdata.Q.*zNdScO3o001(5).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zNdScO3o001(6).*d))...
                +2*(fAsubstrate+fO).*exp(1i*(plotdata.Q.*zNdScO3o001(7).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zNdScO3o001(8).*d));
            case '(110)o'%(001)pc
                h=1; k=1; l=0;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                V=a*b*c/4; %volume of the unit cell
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(101)o'%(111)pc
                h=1; k=0; l=1;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                d=d/2;
                V=a*b*c/4; %volume of the unit cell
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end
        
    case 'PMN-PT'%Tetragonal
        %PMN71PT29t001 tetragonal 001 a=b=4.0166 c=4.02870 (PMN71PT29_P4mm_155863.cif)
        %a=4.0166; b=a; c=4.02870;
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
        
    case 'Si'%Cubic 
        V=a*b*c; %unit cell volume
        %8Si
        fSi=AtomicScatteringFactor('Si',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(001)c'
                h=0; k=0; l=1;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));%=5.4307;
                %rotated in-plane so one unit-cell surface of silicon correspond to 2 unit-cell surface of traditional perovskite
                zSi(1)=0;                     %2Si
                zSi(2)=0.25;                  %2Si
                zSi(3)=0.5;                   %2Si
                zSi(4)=0.75;                  %2Si
                Fsubstrate=2*fSi.*exp(1i*(plotdata.Q.*zSi(1).*d))...
                    +2*fSi.*exp(1i*(plotdata.Q.*zSi(2).*d))...
                    +2*fSi.*exp(1i*(plotdata.Q.*zSi(3).*d))...
                    +2*fSi.*exp(1i*(plotdata.Q.*zSi(4).*d));            
            case '(111)c'
                h=1; k=1; l=1;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                z(1)=0; f(1,:)=2*fSi; 
                z(2)=3/4; f(2,:)=4*fSi;
                z(3)=1;f(3,:)=2*fSi;
                for j=1:3,
                    Fsubstrate=Fsubstrate+f(j,:).*exp(1i*(plotdata.Q.*z(j).*d));
                end;
        end
        
             case 'SmScO3'%orthorhombic
        % 4xSmScO3
        % rotated in-plane so that one unit-cell [ab] of SmScO3 corresponds to 2 unit-cell surface of traditional perovskite
        % orthorhombic with a,b in-plane, c out-of-plane
        %a=5.527; b=5.758; c=7.965;  %alpha=beta=gamma=90 - c out-of-plane, a,b in-plane
        %d=c;                        %out-of-plane lattice parameter
        %V=a*b*d;                    %volume of the unit cell
        fAsubstrate=AtomicScatteringFactor('Sm',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Sc',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(001)o'
                d=c;%out-of-plane lattice parameter
                V=a*b*c;
                zSmScO3o001(1)=0;                   % 2 Sc
                zSmScO3o001(2)=0.062238;            % 2 O
                zSmScO3o001(3)=0.250000;            % 2 Sm - 2 O
                zSmScO3o001(4)=0.437762;            % 2 O
                zSmScO3o001(5)=0.500000;            % 2 Sc
                zSmScO3o001(6)=0.562238;            % 2 O
                zSmScO3o001(7)=0.750000;            % 2 Sm - 2 O
                zSmScO3o001(8)=0.937762;            % 2 O
                Fsubstrate=2*fBsubstrate.*exp(1i*(plotdata.Q.*zSmScO3o001(1).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zSmScO3o001(2).*d))...
                +2*(fAsubstrate+fO).*exp(1i*(plotdata.Q.*zSmScO3o001(3).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zSmScO3o001(4).*d))...
                +2*fBsubstrate.*exp(1i*(plotdata.Q.*zSmScO3o001(5).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zSmScO3o001(6).*d))...
                +2*(fAsubstrate+fO).*exp(1i*(plotdata.Q.*zSmScO3o001(7).*d))...
                +2*fO.*exp(1i*(plotdata.Q.*zSmScO3o001(8).*d));
            case '(110)o'%(001)pc
                h=1; k=1; l=0;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                V=a*b*c/4; %volume of the unit cell
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(101)o'%(111)pc
                h=1; k=0; l=1;
                d=1/(sqrt(h^2/a^2+k^2/b^2+l^2/c^2));
                d=d/2;
                V=a*b*c/4; %volume of the unit cell
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end    
        
    case 'SrLaAlO4'%Tetragnoal
        %1531895 CIF a=3.7544 b=3.7544 c=12.6494 alpha=beta=gamma=90
        %a=3.7544; b=a; c=12.6494;
        %d=12.6377;
        d=c;
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
        
    case 'SrLaGaO4'%Tetragonal
        %SrLaGaO4 tetragonal 001 a=b=3.84 c=12.68 (Crystec datasheets)
        %a=3.84; b=a; c=12.68;
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
        
    case 'SrTiO3'%Cubic
        %d=3.905;
        V=a*b*c; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Sr',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Ti',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(001)c'
                d=c;
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(111)c'                
                d=c/sqrt(3);
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end        
  
    case 'TbScO3'%Orthorhombic
        %(001)o
        %rotated in-plane so that one unit-cell [ab] of TbScO3 corresponds to 2 unit-cell surface of traditional perovskite
        %orthorhombic a=5.4654A b=5.7292A c=7.917A with a,b in-plane, c out-of-plane
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
    
    case 'TiO2'%Tetragonal
        %(001)t
        %TiO2-rutile tetragonal 001 a=b=4.606777A c=2.991757A (Materialsproject.org mp-2657)
        d=c;
        V=a*b*c; %unit cell volume
        zTiO2(1)=0;     % 0.5 Ti - 1 O
        zTiO2(2)=0.5;   % 1 Ti - 2 O
        zTiO2(3)=1;     % 0.5 Ti - 1 O
        fTi=AtomicScatteringFactor('Ti',plotdata.QQ);
        Fsubstrate=(0.5*fTi+fO).*exp(1i*(plotdata.Q.*zTiO2(1).*d))...
            +(fTi+2*fO).*exp(1i*(plotdata.Q.*zTiO2(2).*d))...
            +(0.5*fTi+fO).*exp(1i*(plotdata.Q.*zTiO2(3).*d));
        
    case 'YAlO3'%PseudoCubic 
        %d=3.71;
        d=c;
        V=a*b*c; %unit cell volume
        fAsubstrate=AtomicScatteringFactor('Y',plotdata.QQ);
        fBsubstrate=AtomicScatteringFactor('Al',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(001)pc'
                d=c;
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite001,d); %Structure factor F 
            case '(111)pc'                
                d=c/sqrt(3);
                Fsubstrate=StructureFactor(fAsubstrate,fBsubstrate,fO,plotdata.Q,zperovskite111,d); %Structure factor F 
        end 
        
    case 'YSZ'%Cubic fluorite AB2
        %yttria stabilized zirconia, ZrO2 + 9.5mol% correspond to nominal formula Zr_0.826Y_0.174O_1.913        
        %a=b=c=5.1453 alpha=beta=gamma=90 as measured
        fAsubstrate=AtomicScatteringFactor('Zr',plotdata.QQ)*0.826+AtomicScatteringFactor('Y',plotdata.QQ)*0.174;
        fBsubstrate=AtomicScatteringFactor('O',plotdata.QQ)*1.913/2;
        switch plotdata.Substrate.Orientation        
            case '(001)c'   %fluorite 4xAO2=4x3atoms=12atoms
                d=c;        
                V=a*b*c; %unit cell volume
                Fsubstrate=StructureFactorFluorite001(fAsubstrate,fBsubstrate,plotdata.Q,d); %Structure Factor Fluorite 001
            case '(111)c' %fluorite 3*4xAO2=12x3atoms=36atoms
                d=c/sqrt(3)*6/2;
                V=a*b*c*6/2; %unit cell volume
                Fsubstrate=StructureFactorFluorite111(fAsubstrate,fBsubstrate,plotdata.Q,d); %Structure Factor Fluorite 111
        end
    case 'ZnO'%Hexagonal
        %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
        fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
        switch plotdata.Substrate.Orientation
            case '(0001)h'
                %z along c-axis
                % inplane surface=sqrt(3)/2*a^2
                d=c;
                V=sqrt(3)*a*b*c/2; %unit cell volume
                zZnOh0001(1)=0;     %Zn
                zZnOh0001(2)=0.5;   %Zn
                zZnOh0001(3)=0.3821;%O
                zZnOh0001(4)=0.8821;%O
                Fsubstrate=fZn.*exp(1i*(plotdata.Q.*zZnOh0001(1).*d))...
                +fZn.*exp(1i*(plotdata.Q.*zZnOh0001(2).*d))...
                +fO.*exp(1i*(plotdata.Q.*zZnOh0001(3).*d))...
                +fO.*exp(1i*(plotdata.Q.*zZnOh0001(4).*d));
            case '(11-20)h'
                %z along bissectrice a,b
                % inplane surface S=sqrt(3)*a*c
                % distance between planes = V(unit cell)/S = a/2
                d=a/2;
                V=sqrt(3)*a*b*c/2; %unit cell volume
                zZnOh11bar20(1)=0;     %2 Zn + 2 O
                Fsubstrate=(2*fZn+2*fO).*exp(1i*(plotdata.Q.*zZnOh11bar20(1).*d));
            case '(10-10)h'
                %z perpendicular to a
                % inplane surface S=ac
                % distance between planes = V(unit cell)/S = sqrt(3)a/2
                d=a*sqrt(3)/2;
                V=sqrt(3)*a*b*c/2; %unit cell volume
                zZnOh10bar10(1)=1/3;   % Zn + O
                zZnOh10bar10(2)=2/3;     % Zn + O
                Fsubstrate=(fZn+fO).*exp(1i*(plotdata.Q.*zZnOh10bar10(1).*d))+(fZn+fO).*exp(1i*(plotdata.Q.*zZnOh10bar10(2).*d));
        end
    case 'none'                
        d=0;
        V=1; %unit cell volume
        Fsubstrate=0.*plotdata.Q;
end;
 
        

%% Calculations

expzmsubstrate=zeros(size(plotdata.Q));
h = waitbar(0.2,'Please wait...');
%expzmsubstrate=1./(1-exp(1i*plotdata.Q.*d+plotdata.epsilon));%epsilon:absorption coefficient
for zmsubstrate=-plotdata.Substrate.N:0,
   expzmsubstrate=expzmsubstrate+exp(1i*plotdata.Q.*zmsubstrate.*d).*exp(plotdata.epsilon*zmsubstrate);
end;
close(h) 

Fsubstrate=Fsubstrate*d/V; %renormalized by the in-plane surface of the unit cell

%Scattering amplitude g:
plotdata.Substrate.g=Fsubstrate.*expzmsubstrate;

%toc;

SubstrateData(condition).d=d; %store the out-of-plane lattice parameter
plotdata.Substrate.d=SubstrateData(condition).d;
SubstrateData(condition).V=V; %store the unit cell volume
SubstrateData(condition).g=plotdata.Substrate.g; %store the calculated intensity
save('Substrates.mat','SubstrateData','-append'); %adds new variables to an existing file. If a variable already exists in a MAT-file, then save overwrites it with the value in the workspace.
