% Created by Celine Lichtensteiger
% Calculates the form factors of the different materials depending on their
% crystaollographic structures
% "z" are the out-of-plane atomic positions in the unit cells, expressed 
% in relative unit
%*********************************
function[FLayer]=FLayer(k,d)
global fitdata plotdata;

%% Perovskite
%Glazer-Acta Crystallogr B-1978.pdf

zperovskite001(1)=0;                                                   %A
zperovskite001(2)=0.5+fitdata.Material(k).Polarization*0.162/4.156;    %B
zperovskite001(3)=0+fitdata.Material(k).Polarization*0.473/4.156;      %O1
zperovskite001(4)=0.5+fitdata.Material(k).Polarization*0.486/4.156;    %O2
zperovskite001(5)=0.5+fitdata.Material(k).Polarization*0.486/4.156;    %O3

zperovskite111(1)=0;   %A
zperovskite111(2)=0.5; %B
zperovskite111(3)=0;	%O1
zperovskite111(4)=0;	%O2
zperovskite111(5)=0;	%O3

%% Double perovskite
zdbperovskite111(1)=0;   %A
zdbperovskite111(2)=0.25;%B-1
zdbperovskite111(3)=0;   %O1
zdbperovskite111(4)=0;	  %O2
zdbperovskite111(5)=0;	  %O3
zdbperovskite111(6)=0.5; %A
zdbperovskite111(7)=0.75;%B-2
zdbperovskite111(8)=0.5; %O1
zdbperovskite111(9)=0.5; %O2
zdbperovskite111(10)=0.5;%O3

%% AO:
zAO(1)=0; %A
zAO(2)=0; %O

%% BO2:
zBO2(1)=0; %B
zBO2(2)=0; %O
zBO2(3)=0; %O

fO=AtomicScatteringFactor('O',plotdata.QQ);
fA=zeros(size(plotdata.QQ));
fB=zeros(size(plotdata.QQ));

%% none
if strcmp(fitdata.Material(k).Type,'none'),
    fO=zeros(size(plotdata.QQ));
    FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
end;

%% AlO2
if strcmp(fitdata.Material(k).Type,'AlO2'), %%BO2-type monolayer
       %plotdata.Vunitcell=54.052; % for Al2O4 https://materialsproject.org/materials/mp-1182858/
       %fitdata.Material(k).Sunitcell=plotdata.Material(k).Vunitcell/d;
       plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
       fB=AtomicScatteringFactor('Al',plotdata.QQ);
       FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
       %FLayer=Flayer/plotdata.Material(k).Sunitcell; %renormalized by the surface size to take into account the relative atomic density
       %if strained, the Sunitcell should be the same as the one of the
       %substrate
end;

%% BaBiO3
if strcmp(fitdata.Material(k).Type,'BaBiO3'),
    plotdata.Material(k).V=4.43647^3; %unit cell volume https://materialsproject.org/materials/mp-1227955/
    fA=AtomicScatteringFactor('Ba',plotdata.QQ); 
    fB=AtomicScatteringFactor('Bi',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with BaBiO3 orientation');
    end;
end;

%% BaO
if strcmp(fitdata.Material(k).Type,'BaO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('Ba',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% BaSnO3
if strcmp(fitdata.Material(k).Type,'BaSnO3'),
    plotdata.Material(k).V=4.188634^3; %unit cell volume https://materialsproject.org/materials/mp-3163/
    fA=AtomicScatteringFactor('Ba',plotdata.QQ); 
    fB=AtomicScatteringFactor('Sn',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with BaSnO3 orientation');
    end;
end;

%% (Ba_x,Sr_{1-x})TiO3
if strcmp(fitdata.Material(k).Type,'(Ba_x,Sr_{1-x})TiO3'),
    %BaTiO3: cubic a=4.036 https://materialsproject.org/materials/mp-2998/
    %SrTiO3: cubic a=3.945 https://materialsproject.org/materials/mp-5229/    
    plotdata.Material(k).V=fitdata.Material(k).xBST*4.036^3+(1-fitdata.Material(k).xBST)*3.945^3;%volume of the unit cell
    fA=fitdata.Material(k).xBST*AtomicScatteringFactor('Ba',plotdata.QQ)+(1-fitdata.Material(k).xBST)*AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
    %Sunitcell=Vunitcell/fitdata.Material(k).d; %in-plane area of the unit cell
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
        %FLayer=Flayer/Sunitcell; %renormalized by the surface size to take into account the relative atomic density
        %if strained, the Sunitcell should be the same as the one of the
        %substrate
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
        %FLayer=Flayer/Sunitcell; %renormalized by the surface size to take into account the relative atomic density
        %if strained, the Sunitcell should be the same as the one of the
        %substrate
    else
        msgbox('Error with (Ba_x,Sr_{1-x})TiO3 orientation');
    end;
end;

%% BaTiO3
if strcmp(fitdata.Material(k).Type,'BaTiO3'),
    plotdata.Material(k).V=4.036^3; %unit cell volume https://materialsproject.org/materials/mp-2998/
    fA=AtomicScatteringFactor('Ba',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with BaTiO3 orientation');
    end;
end;

%% BiFeO3
if strcmp(fitdata.Material(k).Type,'BiFeO3'),
    plotdata.Material(k).V=3.753^2*4.924; %unit cell volume https://materialsproject.org/materials/mp-1069079/
    fA=AtomicScatteringFactor('Bi',plotdata.QQ); 
    fB=AtomicScatteringFactor('Fe',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with BiFeO3 orientation');
    end;
end;

%% CaCuO2 (CIF CaCuO2 tetragonal a=b=3.87328 c=3.20546):
if strcmp(fitdata.Material(k).Type,'CaCuO2'),
    plotdata.Material(k).V=3.87328^2*3.20546;
    zCaCuO2(1)=0.5;                     %Ca
    zCaCuO2(2)=0;                       %Cu
    zCaCuO2(3)=0;                       %O
    zCaCuO2(4)=0;                       %O
    fCa=AtomicScatteringFactor('Ca',plotdata.QQ); 
    fCu=AtomicScatteringFactor('Cu',plotdata.QQ);
    FLayer=fCa.*exp(1i*(plotdata.Q.*zCaCuO2(1).*d))...
        +fCu.*exp(1i*(plotdata.Q.*zCaCuO2(2).*d))...
        +fO.*exp(1i*(plotdata.Q.*zCaCuO2(3).*d))...
        +fO.*exp(1i*(plotdata.Q.*zCaCuO2(4).*d));
end;

%% Ca2RuO4 (pseudo-tetragonal a=5.490294 b=5.525868 c=11.968145, Materialsproject.org mp-21466):
if strcmp(fitdata.Material(k).Type,'Ca2RuO4'),
    plotdata.Material(k).V=5.490294*5.525868*11.968145;
    zCa2RuO4(1)=0;                     %Ru
    zCa2RuO4(2)=0.0279150194119473;    %2O
    zCa2RuO4(3)=0.149891984096115;     %2Ca
    zCa2RuO4(4)=0.167109021489964;     %2O
    zCa2RuO4(5)=0.332891020287605;     %2O
    zCa2RuO4(6)=0.350108057681454;     %2Ca
    zCa2RuO4(7)=0.472085022365621;     %2O
    zCa2RuO4(8)=0.500000041777569;     %2Ru
    zCa2RuO4(9)=0.527914977634379;     %2O
    zCa2RuO4(10)=0.649892025873684;    %2Ca
    zCa2RuO4(11)=0.667108979712395;    %2O
    zCa2RuO4(12)=0.832891062065174;    %2O
    zCa2RuO4(13)=0.850108015903885;    %2Ca
    zCa2RuO4(14)=0.972084980588053;    %2O
    zCa2RuO4(15)=1;                    %Ru
    fCa=AtomicScatteringFactor('Ca',plotdata.QQ); 
    fRu=AtomicScatteringFactor('Ru',plotdata.QQ);
    FLayer=fRu.*exp(1i*(plotdata.Q.*zCa2RuO4(1).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCa2RuO4(2).*d))...
        +2*fCa.*exp(1i*(plotdata.Q.*zCa2RuO4(3).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCa2RuO4(4).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCa2RuO4(5).*d))...
        +2*fCa.*exp(1i*(plotdata.Q.*zCa2RuO4(6).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCa2RuO4(7).*d))...
        +2*fRu.*exp(1i*(plotdata.Q.*zCa2RuO4(8).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCa2RuO4(9).*d))...
        +2*fCa.*exp(1i*(plotdata.Q.*zCa2RuO4(10).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCa2RuO4(11).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCa2RuO4(12).*d))...
        +2*fCa.*exp(1i*(plotdata.Q.*zCa2RuO4(13).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCa2RuO4(14).*d))...
        +fRu.*exp(1i*(plotdata.Q.*zCa2RuO4(15).*d));
end;

%% CaTiO3
if strcmp(fitdata.Material(k).Type,'CaTiO3'),
    plotdata.Material(k).V=3.889^3; %unit cell volume https://materialsproject.org/materials/mp-5827/
    fA=AtomicScatteringFactor('Ca',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with CaTiO3 orientation');
    end;
end;

%% CaVO3
if strcmp(fitdata.Material(k).Type,'CaVO3'),
    plotdata.Material(k).V=3.830^3; %unit cell volume https://materialsproject.org/materials/mp-1016853/
    fA=AtomicScatteringFactor('Ca',plotdata.QQ); 
    fB=AtomicScatteringFactor('V',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with CaVO3 orientation');
    end;
end;

%% CoO (CIF CoO cubic a=b=d=4.263 all angles=90?) - idem MgO
if strcmp(fitdata.Material(k).Type,'CoO'),
    plotdata.Material(k).V=4.263^3; 
    fCo=AtomicScatteringFactor('Co',plotdata.QQ); 
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        zCoOc001(1)=0;                      % Co-O
        zCoOc001(2)=0.5;                    % 2Co-2O
        zCoOc001(3)=1;                      % CoO
        FLayer=(fCo+fO).*exp(1i*(plotdata.Q.*zCoOc001(1).*d))...
        +2*(fCo+fO).*exp(1i*(plotdata.Q.*zCoOc001(2).*d))...
        +(fCo+fO).*exp(1i*(plotdata.Q.*zCoOc001(3).*d));  %Structure factor F 
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        zCoOc111(1)=0;                    % 2Co
        zCoOc111(2)=0.5;                  % 4O
        zCoOc111(3)=1;                    % 2Co
        FLayer=2*fCo.*exp(1i*(plotdata.Q.*zCoOc111(1).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zCoOc111(2).*d))...
        +2*fCo.*exp(1i*(plotdata.Q.*zCoOc111(3).*d));  %Structure factor F 
    else
        msgbox('Error with CoO orientation');
    end;
end;


%% LaAlO3
if strcmp(fitdata.Material(k).Type,'LaAlO3'),
    plotdata.Material(k).V=3.787^3; %Wikipedia... (room temperature pseudocubic)
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Al',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaAlO3 orientation');
    end;
end;

%% LaCoO3
if strcmp(fitdata.Material(k).Type,'LaCoO3'),
    plotdata.Material(k).V=3.816^3; %https://materialsproject.org/materials/mp-573180/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Co',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaCoO3 orientation');
    end;
end;

%% La2CuO4 (CIF orthorhombic Cmca a=5.35 b=13.148 c=5.398):
if strcmp(fitdata.Material(k).Type,'La2CuO4'),
    plotdata.Material(k).V=5.35*13.148*5.398;
    zLa2CuO4(1)=0;                      %Cu
    zLa2CuO4(2)=0.00696995740797079;    %O
    zLa2CuO4(3)=0.138890021296015;      %La
    zLa2CuO4(4)=0.18306997261941;       %O
    zLa2CuO4(5)=0.31693002738059;       %O
    zLa2CuO4(6)=0.361109978703985;      %La
    zLa2CuO4(7)=0.493030042592029;      %O
    zLa2CuO4(8)=0.5;                    %Cu
    zLa2CuO4(9)=0.506969957407971;      %O
    zLa2CuO4(10)=0.638890021296015;     %La
    zLa2CuO4(11)=0.683070048676605;     %O
    zLa2CuO4(12)=0.816929951323395;     %O
    zLa2CuO4(13)=0.861109978703985;     %La
    zLa2CuO4(14)=0.993029966534834;     %O
    fLa=AtomicScatteringFactor('La',plotdata.QQ); 
    fCu=AtomicScatteringFactor('Cu',plotdata.QQ);
    FLayer=fCu.*exp(1i*(plotdata.Q.*zLa2CuO4(1).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(2).*d))...
        +fLa.*exp(1i*(plotdata.Q.*zLa2CuO4(3).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(4).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(5).*d))...
        +fLa.*exp(1i*(plotdata.Q.*zLa2CuO4(6).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(7).*d))...
        +fCu.*exp(1i*(plotdata.Q.*zLa2CuO4(8).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(9).*d))...
        +fLa.*exp(1i*(plotdata.Q.*zLa2CuO4(10).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(11).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(12).*d))...
        +fLa.*exp(1i*(plotdata.Q.*zLa2CuO4(13).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(14).*d));
end;

%% LaFeO3
if strcmp(fitdata.Material(k).Type,'LaFeO3'),
    plotdata.Material(k).V=3.959^3; %https://materialsproject.org/materials/mp-552676/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Fe',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaFeO3 orientation');
    end;
end;

%% LaMnO3
if strcmp(fitdata.Material(k).Type,'LaMnO3'),
    plotdata.Material(k).V=3.945^3; %https://materialsproject.org/materials/mp-19025/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Mn',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaMnO3 orientation');
    end;
end;

%% LaNiO3
if strcmp(fitdata.Material(k).Type,'LaNiO3'),
    plotdata.Material(k).V=3.857^3; %https://materialsproject.org/materials/mp-1075921/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaNiO3 orientation');
    end;
end;

%% La2NiMnO6 (double perovskites)
if strcmp(fitdata.Material(k).Type,'La2NiMnO6'),
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        plotdata.Material(k).V=3.871^3; %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802844/
        fA=AtomicScatteringFactor('La',plotdata.QQ); 
        fB=0.5*AtomicScatteringFactor('Ni',plotdata.QQ)+0.5*AtomicScatteringFactor('Mn',plotdata.QQ);
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        plotdata.Material(k).V=3.871^3*2; %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802844/
        fLa=AtomicScatteringFactor('La',plotdata.QQ);
        fNi=AtomicScatteringFactor('Ni',plotdata.QQ); 
        fMn=AtomicScatteringFactor('Mn',plotdata.QQ);
        FLayer=fLa.*exp(1i*(plotdata.Q.*zdbperovskite111(1).*d))...
        +fNi.*exp(1i*(plotdata.Q.*zdbperovskite111(2).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(3).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(4).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(5).*d))...
        +fLa.*exp(1i*(plotdata.Q.*zdbperovskite111(6).*d))...
        +fMn.*exp(1i*(plotdata.Q.*zdbperovskite111(7).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(8).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(9).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(10).*d));
    else
        msgbox('Error with La2NiMnO6 orientation');
    end;
end;

%% LaO
if strcmp(fitdata.Material(k).Type,'LaO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('La',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% LaTiO3
if strcmp(fitdata.Material(k).Type,'LaTiO3'),
    plotdata.Material(k).V=3.959^3; %https://materialsproject.org/materials/mp-8020/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaTiO3 orientation');
    end;
end;

%% LaVO3
if strcmp(fitdata.Material(k).Type,'LaVO3'),
    plotdata.Material(k).V=3.951^3; %https://materialsproject.org/materials/mp-19053/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('V',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaVO3 orientation');
    end;
end;

%% LSMO
if strcmp(fitdata.Material(k).Type,'LSMO'),
    plotdata.Material(k).V=3.875^3;%http://ematweb.cmi.ua.ac.be/emat/pdf/1214.pdf
    fA=0.67*AtomicScatteringFactor('La',plotdata.QQ)+0.33*AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Mn',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LSMO orientation');
    end;
end;

%% Mg3N2
% Ia-3 a=b=c=9.9528 alpha=bet=gamma=90 80 atoms
if strcmp(fitdata.Material(k).Type,'Mg3N2'),
    fMg=AtomicScatteringFactor('Mg',plotdata.QQ);
    fN=AtomicScatteringFactor('N',plotdata.QQ);
    a=9.9528; b=a; c=a;
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        plotdata.Material(k).V=a^3;
        zMg3N2(1)=0;        %2N
        zMg3N2(2)=0.303461; %2N
        zMg3N2(3)=1.104761; %4Mg
        zMg3N2(4)=1.171445; %4Mg
        zMg3N2(5)=1.512825; %4Mg
        zMg3N2(6)=2.4882;   %8N
        zMg3N2(7)=3.463574; %4Mg
        zMg3N2(8)=3.804955; %4Mg
        zMg3N2(9)=3.871639; %4Mg
        zMg3N2(10)=4.67294; %2N
        zMg3N2(11)=4.9764;  %4N
        zMg3N2(12)=5.27986; %2N
        zMg3N2(13)=6.081161;%4Mg
        zMg3N2(14)=6.147844;%4Mg
        zMg3N2(15)=6.489225;%4Mg
        zMg3N2(16)=7.4646;  %8N
        zMg3N2(17)=8.439974;%4Mg
        zMg3N2(18)=8.781356;%4Mg
        zMg3N2(19)=8.848039;%4Mg
        zMg3N2(20)=9.649339;%2N
        zMg3N2(21)=9.9528;  %2N
        zMg3N2=zMg3N2/c; %to go from absolute to relative units
        FLayer=2*fN.*exp(1i*(plotdata.Q.*zMg3N2(1).*d))...
        +2*fN.*exp(1i*(plotdata.Q.*zMg3N2(2).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(3).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(4).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(5).*d))...
        +8*fN.*exp(1i*(plotdata.Q.*zMg3N2(6).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(7).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(8).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(9).*d))...
        +2*fN.*exp(1i*(plotdata.Q.*zMg3N2(10).*d))...
        +4*fN.*exp(1i*(plotdata.Q.*zMg3N2(11).*d))...
        +2*fN.*exp(1i*(plotdata.Q.*zMg3N2(12).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(13).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(14).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(15).*d))...
        +8*fN.*exp(1i*(plotdata.Q.*zMg3N2(16).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(17).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(18).*d))...
        +4*fMg.*exp(1i*(plotdata.Q.*zMg3N2(19).*d))...
        +2*fN.*exp(1i*(plotdata.Q.*zMg3N2(20).*d))...
        +2*fN.*exp(1i*(plotdata.Q.*zMg3N2(21).*d));
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        %need to use hexagonal unit cell (6*V) -> result in 80*6=480 atoms divided in 12 sub-layers (12Mg/16N/12Mg)
        % -> new unit cell is 2*a*sqrt(3)/12=a/(2*sqrt(3))
        plotdata.Material(k).V=a^3/2;
        zMg3N2(1)=0.15339984728;    %3Mg
        zMg3N2(2)=0.23859979101;    %3Mg
        zMg3N2(3)=0.29059983120;    %3Mg
        zMg3N2(4)=0.31739992766;    %3Mg
        zMg3N2(5)=0.43901997428;    %6N
        zMg3N2(6)=0.50000000000;    %4N
        zMg3N2(7)=0.56097982477;    %6N
        zMg3N2(8)=0.68259987139;    %3Mg
        zMg3N2(9)=0.70939996785;    %3Mg
        zMg3N2(10)=0.76140000804;   %3Mg
        zMg3N2(11)=0.84659995177;   %3Mg       
        FLayer=3*fMg.*exp(1i*(plotdata.Q.*zMg3N2(1).*d))...
        +3*fMg.*exp(1i*(plotdata.Q.*zMg3N2(2).*d))...
        +3*fMg.*exp(1i*(plotdata.Q.*zMg3N2(3).*d))...
        +3*fMg.*exp(1i*(plotdata.Q.*zMg3N2(4).*d))...
        +6*fN.*exp(1i*(plotdata.Q.*zMg3N2(5).*d))...
        +4*fN.*exp(1i*(plotdata.Q.*zMg3N2(6).*d))...
        +6*fN.*exp(1i*(plotdata.Q.*zMg3N2(7).*d))...
        +3*fMg.*exp(1i*(plotdata.Q.*zMg3N2(8).*d))...
        +3*fMg.*exp(1i*(plotdata.Q.*zMg3N2(9).*d))...
        +3*fMg.*exp(1i*(plotdata.Q.*zMg3N2(10).*d))...
        +3*fMg.*exp(1i*(plotdata.Q.*zMg3N2(11).*d));             
    else
        msgbox('Error with Mg3N2 orientation');
    end;
end;

%% MgO
if strcmp(fitdata.Material(k).Type,'MgO'),
    % cubic a=4.212 A
    % (Crystec Datasheets http://www.crystec.de/daten/mgo.pdf)
    plotdata.Material(k).V=4.212^3; 
    fMg=AtomicScatteringFactor('Mg',plotdata.QQ); 
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        zMgOc001(1)=0;                      % Mg-O
        zMgOc001(2)=0.5;                    % 2Mg-2O
        zMgOc001(3)=1;                      % MgO
        FLayer=(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(1).*d))...
        +2*(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(2).*d))...
        +(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(3).*d));  %Structure factor F 
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        zMgOc111(1)=0;                    % 2Mg
        zMgOc111(2)=0.5;                  % 4O
        zMgOc111(3)=1;                    % 2Mg
        FLayer=2*fMg.*exp(1i*(plotdata.Q.*zMgOc111(1).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zMgOc111(2).*d))...
        +2*fMg.*exp(1i*(plotdata.Q.*zMgOc111(3).*d));  %Structure factor F 
    else
        msgbox('Error with MgO orientation');
    end;
end;

%% MnO
if strcmp(fitdata.Material(k).Type,'MnO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('Mn',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% MnO2
if strcmp(fitdata.Material(k).Type,'MnO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('Mn',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% MnTiO3
if strcmp(fitdata.Material(k).Type,'MnTiO3'),
    plotdata.Material(k).V=3.832^3; %https://materialsproject.org/materials/mp-19082/ V=112.583 for double unit cell
    fA=AtomicScatteringFactor('Mn',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with MnTiO3 orientation');
    end;
end;

%% (Nd_x,La_{1-x})NiO3
if strcmp(fitdata.Material(k).Type,'(Nd_x,La_{1-x})NiO3'),
    %NdNiO3: pseudocubic a=3.861 cubicroot(230.143/4) https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
    %LaNiO3: cubic a=3.857; %https://materialsproject.org/materials/mp-1075921/
    plotdata.Material(k).V=fitdata.Material(k).xNLN*3.861^3+(1-fitdata.Material(k).xNLN)*3.857^3;%volume of the unit cell    
    fA=fitdata.Material(k).xNLN*AtomicScatteringFactor('Nd',plotdata.QQ)+(1-fitdata.Material(k).xNLN)*AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with (Nd_x,La_{1-x})NiO3 orientation');
    end;
end;

%% NdNiO2 (cf CaCuO2) 
if strcmp(fitdata.Material(k).Type,'NdNiO2'), 
    % tetragonal a=b=3.962, c=3.268 https://materialsproject.org/materials/mp-31063/
    a=3.962; b=a; c=3.268;
    plotdata.Material(k).V=a*b*c;
    zNdNiO2(1)=0.5;                     %Nd
    zNdNiO2(2)=0;                       %Ni
    zNdNiO2(3)=0;                       %O
    zNdNiO2(4)=0;                       %O
    fNd=AtomicScatteringFactor('Nd',plotdata.QQ); 
    fNi=AtomicScatteringFactor('Ni',plotdata.QQ);
    FLayer=fNd.*exp(1i*(plotdata.Q.*zNdNiO2(1).*d))...
        +fNi.*exp(1i*(plotdata.Q.*zNdNiO2(2).*d))...
        +fO.*exp(1i*(plotdata.Q.*zNdNiO2(3).*d))...
        +fO.*exp(1i*(plotdata.Q.*zNdNiO2(4).*d));
end;

%% NdNiO3
if strcmp(fitdata.Material(k).Type,'NdNiO3'),
    plotdata.Material(k).V=3.861^3; %https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
    fA=AtomicScatteringFactor('Nd',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with NdNiO3 orientation');
    end;
end;

%% Nd2NiMnO6 (double perovskite)
if strcmp(fitdata.Material(k).Type,'Nd2NiMnO6'),%double perovskite
    %https://aip.scitation.org/doi/full/10.1063/1.4906989?ver=pdfcov
    %a=5.4162, b=5.4863, c=7.6718 V=a*b*c 4 unit cells
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        plotdata.Material(k).V=3.851^3;
        fA=AtomicScatteringFactor('Nd',plotdata.QQ); 
        fB=0.5*AtomicScatteringFactor('Ni',plotdata.QQ)+0.5*AtomicScatteringFactor('Mn',plotdata.QQ);
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        plotdata.Material(k).V=3.851^3*2;
        fNd=AtomicScatteringFactor('Nd',plotdata.QQ);
        fNi=AtomicScatteringFactor('Ni',plotdata.QQ); 
        fMn=AtomicScatteringFactor('Mn',plotdata.QQ);
        fO=AtomicScatteringFactor('O',plotdata.QQ);
        FLayer=fNd.*exp(1i*(plotdata.Q.*zdbperovskite111(1).*d))...
        +fNi.*exp(1i*(plotdata.Q.*zdbperovskite111(2).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(3).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(4).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(5).*d))...
        +fNd.*exp(1i*(plotdata.Q.*zdbperovskite111(6).*d))...
        +fMn.*exp(1i*(plotdata.Q.*zdbperovskite111(7).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(8).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(9).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(10).*d));
    else
        msgbox('Error with Nd2NiMnO6 orientation');
    end;
end;

% NdO
if strcmp(fitdata.Material(k).Type,'NdO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('Nd',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% NiO2
if strcmp(fitdata.Material(k).Type,'NiO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% PbO
if strcmp(fitdata.Material(k).Type,'PbO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('Pb',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% PbNiO3
if strcmp(fitdata.Material(k).Type,'PbNiO3'),
    plotdata.Material(k).V=55.330; %https://materialsproject.org/materials/mp-974108/ V=55.330
    fA=AtomicScatteringFactor('Pb',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with PbNiO3 orientation');
    end;
end;

%% (Pb_x,Sr_{1-x})TiO3
if strcmp(fitdata.Material(k).Type,'(Pb_x,Sr_{1-x})TiO3'),
    %PbTiO3: known bulk values a=b=3.904, c=4.152
    %SrTiO3: known bulk values a=b=c=3.905
    plotdata.Material(k).V=fitdata.Material(k).xPST*(3.904^2*4.152)+(1-fitdata.Material(k).xPST)*3.905^3;%volume of the unit cell   
    fA=fitdata.Material(k).xPST*AtomicScatteringFactor('Pb',plotdata.QQ)+(1-fitdata.Material(k).xPST)*AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with (Pb_x,Sr_{1-x})TiO3 orientation');
    end;
end;

%% PbTiO3
if strcmp(fitdata.Material(k).Type,'PbTiO3'),
    plotdata.Material(k).V=3.904^2*4.152;
    fA=AtomicScatteringFactor('Pb',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with PbTiO3 orientation');
    end;
end;

%% Pb(Zr_x,Ti_{1-x})O3
if strcmp(fitdata.Material(k).Type,'Pb(Zr_x,Ti_{1-x})O3'),
    %PTO: 3.904^2*4.152 PZO:4.209^3 https://materialsproject.org/materials/mp-1068577/
    plotdata.Material(k).V=fitdata.Material(k).xPZT*(4.209^3)+(1-fitdata.Material(k).xPZT)*(3.904^2*4.152); %volume of the unit cell   
    fA=AtomicScatteringFactor('Pb',plotdata.QQ); 
    fB=fitdata.Material(k).xPZT*AtomicScatteringFactor('Zr',plotdata.QQ)+(1-fitdata.Material(k).xPZT)*AtomicScatteringFactor('Ti',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with Pb(Zr_x,Ti_{1-x})O3 orientation');
    end;
end;

%% PrBa2Cu3O7
if strcmp(fitdata.Material(k).Type,'PrBa2Cu3O7'),
    % https://materialsproject.org/materials/mp-20936/
    % orthorhombic alpha=beta=gamma=90
    a=3.956; b=3.904; c=11.916;
    % c out-of-plane
    plotdata.Material(k).V=a*b*c; %unit cell volume
    zPBCO(1)=0.5;       %Pb
    zPBCO(2)=0.1789;	%Ba1
    zPBCO(3)=0.8211;	%Ba2
    zPBCO(4)=0;         %Cu1
    zPBCO(5)=0.3459;	%Cu2
    zPBCO(6)=0.6541;	%Cu3
    zPBCO(7)=0;         %O1
    zPBCO(8)=0.1573;	%O2
    zPBCO(9)=0.3702;	%O3
    zPBCO(10)=0.3702;	%O4
    zPBCO(11)=0.6298;	%O5
    zPBCO(12)=0.6298;	%O6
    zPBCO(13)=0.8427;	%O7
    fPr=AtomicScatteringFactor('Pr',plotdata.QQ);
    fBa=AtomicScatteringFactor('Ba',plotdata.QQ);
    fCu=AtomicScatteringFactor('Cu',plotdata.QQ);
    FLayer=fPr.*exp(1i*(plotdata.Q.*zPBCO(1).*d))...
        +fBa.*exp(1i*(plotdata.Q.*zPBCO(2).*d))...
        +fBa.*exp(1i*(plotdata.Q.*zPBCO(3).*d))...
        +fCu.*exp(1i*(plotdata.Q.*zPBCO(4).*d))...
        +fCu.*exp(1i*(plotdata.Q.*zPBCO(5).*d))...
        +fCu.*exp(1i*(plotdata.Q.*zPBCO(6).*d))...
        +fO.*exp(1i*(plotdata.Q.*zPBCO(7).*d))...
        +fO.*exp(1i*(plotdata.Q.*zPBCO(8).*d))...
        +fO.*exp(1i*(plotdata.Q.*zPBCO(9).*d))...
        +fO.*exp(1i*(plotdata.Q.*zPBCO(10).*d))...
        +fO.*exp(1i*(plotdata.Q.*zPBCO(11).*d))...
        +fO.*exp(1i*(plotdata.Q.*zPBCO(12).*d))...
        +fO.*exp(1i*(plotdata.Q.*zPBCO(13).*d));
end;

%% PrNiO2 (cf CaCuO2) 
if strcmp(fitdata.Material(k).Type,'PrNiO2'), 
    % using the lattice parameters of SrCuO2 tetragonal a=b=3.948, c=3.485 https://materialsproject.org/materials/mp-37514/
    a=3.948; b=a; c=3.485;
    plotdata.Material(k).V=a*b*c;
    zPrNiO2(1)=0.5;                     %Nd
    zPrNiO2(2)=0;                       %Ni
    zPrNiO2(3)=0;                       %O
    zPrNiO2(4)=0;                       %O
    fPr=AtomicScatteringFactor('Pr',plotdata.QQ); 
    fNi=AtomicScatteringFactor('Ni',plotdata.QQ);
    FLayer=fPr.*exp(1i*(plotdata.Q.*zPrNiO2(1).*d))...
        +fNi.*exp(1i*(plotdata.Q.*zPrNiO2(2).*d))...
        +fO.*exp(1i*(plotdata.Q.*zPrNiO2(3).*d))...
        +fO.*exp(1i*(plotdata.Q.*zPrNiO2(4).*d));
end;

%% PrNiO3
if strcmp(fitdata.Material(k).Type,'PrNiO3'),
    plotdata.Material(k).V=3.872^3; %volume of the unit cell - https://materialsproject.org/materials/mp-22280/ V=232.262 4 unit cells
    fA=AtomicScatteringFactor('Pr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with PrNiO3 orientation');
    end;
end;

%% PrVO3
if strcmp(fitdata.Material(k).Type,'PrVO3'),
    plotdata.Material(k).V=3.936^3; %volume of the unit cell - https://materialsproject.org/materials/mp-1069346/
    fA=AtomicScatteringFactor('Pr',plotdata.QQ); 
    fB=AtomicScatteringFactor('V',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with PrVO3 orientation');
    end;
end;

%% RuO2
if strcmp(fitdata.Material(k).Type,'RuO2'), %%BO2-type monolayer
       plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
       fB=AtomicScatteringFactor('Ru',plotdata.QQ);
       FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% SmNiO3
if strcmp(fitdata.Material(k).Type,'SmNiO3'),
    plotdata.Material(k).V=3.794^3; %https://materialsproject.org/materials/mp-1099668/
    fA=AtomicScatteringFactor('Sm',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SmNiO3 orientation');
    end;
end;

%% Sr3Al2O6 https://materialsproject.org/materials/mp-3393/  
%Sr72Al48O144 in cubic unit cell a=b=c=15.999914A
if strcmp(fitdata.Material(k).Type,'Sr3Al2O6'),
    %In the (001) orientation, in-plane 16 unit cells of STO will fit in 1 unit cell of Sr72Al48O144. This is taken into account by divinding by 16 the StructureFactor.
    %Out-of-plane, I just consider half of the Sr72Al48O144 unit cell, i.e. with a bulk lattice parameter of 7.999957A (that the user can then modify). 
    a=15.999914;
    plotdata.Material(k).V=a^3/2; 
    zSAO(1)=0;                      %2 Sr
    zSAO(2)=0.00283001521133176*2;  %2 O
    zSAO(3)=0.0041260221773692*2;   %2 Al
    zSAO(4)=0.0076369785487597*2;   %2 O
    zSAO(5)=0.0159040229841235*2;   %2 Al
    zSAO(6)=0.0168720281871515*2;   %2 O
    zSAO(7)=0.0175839695138361*2;   %2 O
    zSAO(8)=0.0199369821612791*2;   %2 O
    zSAO(9)=0.0200889829782835*2;   %2 Al
    zSAO(10)=0.0260960152660821*2;  %2 O
    zSAO(11)=0.094133005964907*2;   %2 O
    zSAO(12)=0.100253976365123*2;   %2 O
    zSAO(13)=0.115953998252741*2;   %2 Sr
    zSAO(14)=0.119861019252978*2;   %2 Sr
    zSAO(15)=0.124367980977898*2;   %2 Sr
    zSAO(16)=0.124571982074404*2;   %2 Sr
    zSAO(17)=0.124742982993534*2;   %2 Sr
    zSAO(18)=0.125220985562797*2;   %2 O
    zSAO(19)=0.129366007842292*2;   %2 O
    zSAO(20)=0.131721020500485*2;   %2 Sr
    zSAO(21)=0.134670036351445*2;   %2 Sr
    zSAO(22)=0.153380011917564*2;   %2 O
    zSAO(23)=0.153955015008206*2;   %2 O
    zSAO(24)=0.220413997225235*2;   %2 O
    zSAO(25)=0.222892010544557*2;   %2 O
    zSAO(26)=0.230195987303432*2;   %2 O
    zSAO(27)=0.233955007508165*2;   %2 Al
    zSAO(28)=0.234410009953803*2;   %2 O
    zSAO(29)=0.236916023423626*2;   %2 O
    zSAO(30)=0.238600032475175*2;   %2 Al
    zSAO(31)=0.244254000365252*2;   %2 O
    zSAO(32)=0.247728019038102*2;   %2 Sr
    zSAO(33)=0.248376022521121*2;   %2 Al
    zSAO(34)=0.251623977478879*2;   %2 Al
    zSAO(35)=0.252272043462234*2;   %2 Sr
    zSAO(36)=0.255745999634748*2;   %2 O
    zSAO(37)=0.261400030025161*2;   %2 Al
    zSAO(38)=0.263083976576374*2;   %2 O
    zSAO(39)=0.265589990046196*2;   %2 O
    zSAO(40)=0.266044992491835*2;   %2 Al
    zSAO(41)=0.269804012696568*2;   %2 O
    zSAO(42)=0.277107989455443*2;   %2 O
    zSAO(43)=0.279586002774765*2;   %2 O
    zSAO(44)=0.346044984991794*2;   %2 O
    zSAO(45)=0.346619988082436*2;   %2 O
    zSAO(46)=0.365330026148891*2;   %2 Sr
    zSAO(47)=0.368279041999851*2;   %2 Sr
    zSAO(48)=0.370633992157708*2;   %2 O
    zSAO(49)=0.374779014437203*2;   %2 O
    zSAO(50)=0.37525695450613*2;    %2 Sr
    zSAO(51)=0.37542795542526*2;    %2 Sr
    zSAO(52)=0.375631956521766*2;   %2 Sr
    zSAO(53)=0.380138980747022*2;   %2 Sr
    zSAO(54)=0.384046001747259*2;   %2 Sr
    zSAO(55)=0.399744961129166*2;   %2 O
    zSAO(56)=0.405866994035093*2;   %2 O
    zSAO(57)=0.473904047234254*2;   %2 O
    zSAO(58)=0.479911017021717*2;   %2 Al
    zSAO(59)=0.480063017838721*2;   %2 O
    zSAO(60)=0.482416030486164*2;   %2 O
    zSAO(61)=0.483128034313184*2;   %2 O
    zSAO(62)=0.484095977015876*2;   %2 Al
    zSAO(63)=0.49236302145124*2;    %2 O
    zSAO(64)=0.495873977822631*2;   %2 Al
    zSAO(65)=0.497170047289004*2;   %2 O
    zSAO(66)=0.5*2;                 %2 Sr
    fSr=AtomicScatteringFactor('Sr',plotdata.QQ);
    fAl=AtomicScatteringFactor('Al',plotdata.QQ); 
    fO=AtomicScatteringFactor('O',plotdata.QQ);
    FLayer=+2*fSr.*exp(1i*(plotdata.Q.*zSAO(1).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(2).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(3).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(4).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(5).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(6).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(7).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(8).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(9).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(10).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(11).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(12).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(13).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(14).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(15).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(16).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(17).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(18).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(19).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(20).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(21).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(22).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(23).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(24).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(25).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(26).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(27).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(28).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(29).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(30).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(31).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(32).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(33).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(34).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(35).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(36).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(37).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(38).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(39).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(40).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(41).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(42).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(43).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(44).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(45).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(46).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(47).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(48).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(49).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(50).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(51).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(52).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(53).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(54).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(55).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(56).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(57).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(58).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(59).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(60).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(61).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(62).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(63).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zSAO(64).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSAO(65).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSAO(66).*d));
end;

%% SrCoO2.5 SrCoO2_5_Pnma_Brownmillerite_162239.cif orthorhombic with a=5.45790, b=15.6388, c=5.5643
%Consider only up to y=15.6388/2 (symmetric) i.e. 36/2 atomes = 18 atomes
if strcmp(fitdata.Material(k).Type,'SrCoO2.5'),
    %In the (001) orientation, in-plane 2 unit cells of STO will fit in 1 unit cell of SrCoO2.5. This is taken into account by divinding by 2 the StructureFactor.
    %Out-of-plane, I just consider half of the SrCoO2.5 unit cell, i.e. with a bulk lattice parameter of 7.8194A (that the user can then modify). 
    a=5.45790; b=15.6388/2; c=5.5643;
    plotdata.Material(k).V=a*b*c;
    zSCO(1)=0;          %2 Co
    zSCO(2)=0.006000;   %2 O
    zSCO(3)=0.110200;   %2 Sr
    zSCO(4)=0.140900;   %2 O
    zSCO(5)=0.140900;   %2 O + 2 Co
    zSCO(6)=0.359100;   %2 O
    zSCO(7)=0.389800;   %2 Sr
    zSCO(8)=0.494000;   %2 O
    fSr=AtomicScatteringFactor('Sr',plotdata.QQ);
    fCo=AtomicScatteringFactor('Co',plotdata.QQ); 
    fO=AtomicScatteringFactor('O',plotdata.QQ);
    FLayer=2*fCo.*exp(1i*(plotdata.Q.*zSCO(1).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSCO(2).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSCO(3).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSCO(4).*d))...
        +2*(fO+fCo).*exp(1i*(plotdata.Q.*zSCO(5).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSCO(6).*d))...
        +2*fSr.*exp(1i*(plotdata.Q.*zSCO(7).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zSCO(8).*d));
end;

%% SrCoO3
if strcmp(fitdata.Material(k).Type,'SrCoO3'),
    %https://materialsproject.org/materials/mp-505766/ alpha=beta=gamma=90
    a=3.860; b=3.853; c=b; 
    plotdata.Material(k).V=a*b*c;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Co',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrCoO3 orientation');
    end;
end;

%% SrCrO3
if strcmp(fitdata.Material(k).Type,'SrCrO3'),
    plotdata.Material(k).V=3.8185^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Cr',plotdata.QQ);    
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrCrO3 orientation');
    end;
end;

%% SrCuO2 (cf CaCuO2) 
if strcmp(fitdata.Material(k).Type,'SrCuO2'), 
    % tetragonal a=b=3.948, c=3.485 https://materialsproject.org/materials/mp-37514/
    a=3.948; b=a; c=3.485;
    plotdata.Material(k).V=a*b*c;
    zSrCuO2(1)=0.5;                     %Sr
    zSrCuO2(2)=0;                       %Cu
    zSrCuO2(3)=0;                       %O
    zSrCuO2(4)=0;                       %O
    fSr=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fCu=AtomicScatteringFactor('Cu',plotdata.QQ);
    FLayer=fSr.*exp(1i*(plotdata.Q.*zSrCuO2(1).*d))...
        +fCu.*exp(1i*(plotdata.Q.*zSrCuO2(2).*d))...
        +fO.*exp(1i*(plotdata.Q.*zSrCuO2(3).*d))...
        +fO.*exp(1i*(plotdata.Q.*zSrCuO2(4).*d));
end;

%% SrIrO3
if strcmp(fitdata.Material(k).Type,'SrIrO3'),
    %https://materialsproject.org/materials/mp-1016848/
    %cubic
    a=3.998;
    plotdata.Material(k).V=a^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ir',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrIrO3 orientation');
    end;
end;

%% SrMoO3
if strcmp(fitdata.Material(k).Type,'SrMoO3'),
    %https://materialsproject.org/materials/mp-18747/
    %cubic
    a=4.082;
    plotdata.Material(k).V=a^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Mo',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrMoO3 orientation');
    end;
end;

%% SrO
if strcmp(fitdata.Material(k).Type,'SrO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('Sr',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% SrO2
if strcmp(fitdata.Material(k).Type,'SrO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('Sr',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% SrRuO3
if strcmp(fitdata.Material(k).Type,'SrRuO3'),
    %https://materialsproject.org/materials/mp-4346/
    %cubic
    a=3.985;
    plotdata.Material(k).V=a^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ru',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrRuO3 orientation');
    end;
end;

%% SrTiO3
if strcmp(fitdata.Material(k).Type,'SrTiO3'),
    %cubic
    a=3.905;
    plotdata.Material(k).V=a^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrTiO3 orientation');
    end;
end;

%% SrVO3
if strcmp(fitdata.Material(k).Type,'SrVO3'),
    %https://materialsproject.org/materials/mp-18717/
    %cubic
    a=3.901;
    plotdata.Material(k).V=a^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('V',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrVO3 orientation');
    end;
end;

%% TiO2
if strcmp(fitdata.Material(k).Type,'TiO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% Tm3Fe5O12
if strcmp(fitdata.Material(k).Type,'Tm3Fe5O12'),
    % Thulium Iron garnet (TmIG)
    a=12.325; b=a; c=a; %Landolt-Bornstein
    plotdata.Material(k).V=a*b*c/4; %volume of the "new" unit cell (see below)
    fTm=AtomicScatteringFactor('Tm',plotdata.QQ); 
    fFe=AtomicScatteringFactor('Fe',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        % one formula unit contains 20 atoms
        % conventional unit cell contains 160 atoms = 8 formula units
        % "new" unit cell contains 40 atoms = conventional unit cell/4
        %d=c/4;                          %out-of-plane lattice parameter of the "new" unit cell
        zTmIG001(1)=0*4;              % 4 Fe + 2 Tm
        zTmIG001(2)=0.024741*4;       % 4 O
        zTmIG001(3)=0.058764*4;       % 4 O
        zTmIG001(4)=0.098178*4;       % 4 O
        zTmIG001(5)=0.125000*4;       % 2 Fe + 2 Tm
        zTmIG001(6)=0.151822*4;       % 4 O
        zTmIG001(7)=0.191236*4;       % 4 O
        zTmIG001(8)=0.225259*4;       % 4 O
        zTmIG001(9)=0.250000*4;       % 4 Fe + 2 Tm
        FLayer=(4*fFe+2*fTm).*exp(1i*(plotdata.Q.*zTmIG001(1).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zTmIG001(2).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zTmIG001(3).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zTmIG001(4).*d))...
        +(2*fFe+2*fTm).*exp(1i*(plotdata.Q.*zTmIG001(5).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zTmIG001(6).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zTmIG001(7).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zTmIG001(8).*d))...
        +(4*fFe+2*fTm).*exp(1i*(plotdata.Q.*zTmIG001(9).*d));
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        % conventional 111 unit cell contains 160 * 6 atoms = 8*6 formula units
        % "new" unit cell contains 40 atoms = conventional 111 unit cell/(4*6)
        %conventional 001 unit cell = 160 atoms
        %d=c/4/sqrt(3);             %out-of-plane lattice parameter of the "new" unit cell
        zTmIG111(1)=0*24;         % 2 Fe
        zTmIG111(2)=0.0024455*24; % 3 O
        zTmIG111(3)=0.0106925*24; % 3 O
        zTmIG111(4)=0.0113862*24; % 3 O
        zTmIG111(5)=0.0196332*24; % 3 O
        zTmIG111(6)=0.0208333*24; % 6 Tm + 6 Fe
        zTmIG111(7)=0.0220335*24; % 3 O
        zTmIG111(8)=0.0302805*24; % 3 O
        zTmIG111(9)=0.0309742*24; % 3 0
        zTmIG111(10)=0.0392212*24;% 3 0
        zTmIG111(11)=0.0416667*24;% 2 Fe
        FLayer=2*fFe.*exp(1i*(plotdata.Q.*zTmIG111(1).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zTmIG111(2).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zTmIG111(3).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zTmIG111(4).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zTmIG111(5).*d))...
        +(6*fTm+6*fFe).*exp(1i*(plotdata.Q.*zTmIG111(6).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zTmIG111(7).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zTmIG111(8).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zTmIG111(9).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zTmIG111(10).*d))...
        +2*fFe.*exp(1i*(plotdata.Q.*zTmIG111(11).*d));
    else
        msgbox('Error with Tm3Fe5O12 orientation');
    end;
end;

%% VO2
if strcmp(fitdata.Material(k).Type,'VO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('V',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% YBCO:
if strcmp(fitdata.Material(k).Type,'YBa2Cu3O7'),
    %https://materialsproject.org/materials/mp-20674/
    % tetragonal 
    a=3.845; b=3.926; c=11.824;
    plotdata.Material(k).V=a*b*c; %unit cell volume
    zYBCO(1)=0.5;       %Y
    zYBCO(2)=0.1843;	%B1
    zYBCO(3)=0.815700;	%B2
    zYBCO(4)=0;         %C1
    zYBCO(5)=0.355600;	%C2
    zYBCO(6)=0.644400;	%C3
    zYBCO(7)=0;         %O1
    zYBCO(8)=0.3779;	%O2
    zYBCO(9)=0.622100;	%O3
    zYBCO(10)=0.377900;	%O4
    zYBCO(11)=0.622100;	%O5
    zYBCO(12)=0.159000;	%O6
    zYBCO(13)=0.841000;	%O7
    fY=AtomicScatteringFactor('Y',plotdata.QQ); 
    fBa=AtomicScatteringFactor('Ba',plotdata.QQ);
    fCu=AtomicScatteringFactor('Cu',plotdata.QQ);
    FLayer=fY.*exp(1i*(plotdata.Q.*zYBCO(1).*d))...
        +fBa.*exp(1i*(plotdata.Q.*zYBCO(2).*d))...
        +fBa.*exp(1i*(plotdata.Q.*zYBCO(3).*d))...
        +fCu.*exp(1i*(plotdata.Q.*zYBCO(4).*d))...
        +fCu.*exp(1i*(plotdata.Q.*zYBCO(5).*d))...
        +fCu.*exp(1i*(plotdata.Q.*zYBCO(6).*d))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(7).*d))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(8).*d))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(9).*d))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(10).*d))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(11).*d))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(12).*d))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(13).*d));
end;


%% YBiO3
if strcmp(fitdata.Material(k).Type,'YBiO3'),
    plotdata.Material(k).V=4.408542^3; %unit cell volume https://materialsproject.org/materials/mp-13598/
    fA=AtomicScatteringFactor('Y',plotdata.QQ); 
    fB=AtomicScatteringFactor('Bi',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with YBiO3 orientation');
    end;
end;

%% Y3Fe5O12
if strcmp(fitdata.Material(k).Type,'Y3Fe5O12'),
    % Ytrium Iron garnet (YIG)
    a=12.376; b=a; c=a; %Landolt-Bornstein
    plotdata.Material(k).V=a*b*c/4; %volume of the "new" unit cell (see below)
    fY=AtomicScatteringFactor('Y',plotdata.QQ); 
    fFe=AtomicScatteringFactor('Fe',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        % one formula unit contains 20 atoms
        % conventional unit cell contains 160 atoms = 8 formula units
        % "new" unit cell contains 40 atoms = conventional unit cell/4
        % d=c/4;                       %out-of-plane lattice parameter of the "new" unit cell
        zYIG001(1)=0*4;              % 4 Fe + 2 Y
        zYIG001(2)=0.026003*4;       % 4 O
        zYIG001(3)=0.057201*4;       % 4 O
        zYIG001(4)=0.098942*4;       % 4 O
        zYIG001(5)=0.125000*4;       % 2 Fe + 2 Y
        zYIG001(6)=0.151058*4;       % 4 O
        zYIG001(7)=0.192799*4;       % 4 O
        zYIG001(8)=0.223997*4;       % 4 O
        zYIG001(9)=0.250000*4;       % 4 Fe + 2 Y
        FLayer=(4*fFe+2*fY).*exp(1i*(plotdata.Q.*zYIG001(1).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYIG001(2).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYIG001(3).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYIG001(4).*d))...
        +(2*fFe+2*fY).*exp(1i*(plotdata.Q.*zYIG001(5).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYIG001(6).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYIG001(7).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYIG001(8).*d))...
        +(4*fFe+2*fY).*exp(1i*(plotdata.Q.*zYIG001(9).*d));
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        % conventional 111 unit cell contains 160 * 6 atoms = 8*6 formula units
        % "new" unit cell contains 40 atoms = conventional 111 unit cell/(4*6)
        % conventional 001 unit cell = 160 atoms
        % d=c/4/sqrt(3);           %out-of-plane lattice parameter of the "new" unit cell
        zYIG111(1)=0*24;         % 2 Fe
        zYIG111(2)=0.0026230*24; % 3 O
        zYIG111(3)=0.0112907*24; % 3 O
        zYIG111(4)=0.0113090*24; % 3 O
        zYIG111(5)=0.0199767*24; % 3 O
        zYIG111(6)=0.0208333*24; % 6 Y + 6 Fe
        zYIG111(7)=0.0216900*24; % 3 O
        zYIG111(8)=0.0303577*24; % 3 O
        zYIG111(9)=0.0303760*24; % 3 0
        zYIG111(10)=0.0390437*24;% 3 0
        zYIG111(11)=0.0416667*24;% 2 Fe
        FLayer=2*fFe.*exp(1i*(plotdata.Q.*zYIG111(1).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYIG111(2).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYIG111(3).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYIG111(4).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYIG111(5).*d))...
        +(6*fY+6*fFe).*exp(1i*(plotdata.Q.*zYIG111(6).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYIG111(7).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYIG111(8).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYIG111(9).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYIG111(10).*d))...
        +2*fFe.*exp(1i*(plotdata.Q.*zYIG111(11).*d));
    else
        msgbox('Error with Y3Fe5O12 orientation');
    end;
end;

%% (Y_x,Tm_{3-x})Fe5O12
if strcmp(fitdata.Material(k).Type,'(Y_xTm_{3-x})Fe5O12'),
    %YFe5O12: pseudo-cubic a=12.376 Landolt/Bornstein
    %TmFe5O12: pseudo- cubic a=12.325 Landolt/Bornstein
    aTmIG=12.325; aYIG=12.376; %Landolt-Bornstein
    VTmIG=aTmIG^3/4; VYIG=aYIG^3/4;
    plotdata.Material(k).V=(fitdata.Material(k).xYTmIG*VYIG+(3-fitdata.Material(k).xYTmIG)*VTmIG)/3; %volume of the "new" unit cell
    fYTm=(1/3)*(fitdata.Material(k).xYTmIG*AtomicScatteringFactor('Y',plotdata.QQ)+(3-fitdata.Material(k).xYTmIG)*AtomicScatteringFactor('Tm',plotdata.QQ)); 
    fFe=AtomicScatteringFactor('Fe',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        zTmIG001(1)=0*4;              % 4 Fe + 2 Tm
        zTmIG001(2)=0.024741*4;       % 4 O
        zTmIG001(3)=0.058764*4;       % 4 O
        zTmIG001(4)=0.098178*4;       % 4 O
        zTmIG001(5)=0.125000*4;       % 2 Fe + 2 Tm
        zTmIG001(6)=0.151822*4;       % 4 O
        zTmIG001(7)=0.191236*4;       % 4 O
        zTmIG001(8)=0.225259*4;       % 4 O
        zTmIG001(9)=0.250000*4;       % 4 Fe + 2 Tm
        zYIG001(1)=0*4;              % 4 Fe + 2 Y
        zYIG001(2)=0.026003*4;       % 4 O
        zYIG001(3)=0.057201*4;       % 4 O
        zYIG001(4)=0.098942*4;       % 4 O
        zYIG001(5)=0.125000*4;       % 2 Fe + 2 Y
        zYIG001(6)=0.151058*4;       % 4 O
        zYIG001(7)=0.192799*4;       % 4 O
        zYIG001(8)=0.223997*4;       % 4 O
        zYIG001(9)=0.250000*4;       % 4 Fe + 2 Y
        zYTmIG001=(fitdata.Material(k).xYTmIG*zYIG001+(3-fitdata.Material(k).xYTmIG)*zTmIG001)/3;
        FLayer=(4*fFe+2*fYTm).*exp(1i*(plotdata.Q.*zYTmIG001(1).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(2).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(3).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(4).*d))...
        +(2*fFe+2*fYTm).*exp(1i*(plotdata.Q.*zYTmIG001(5).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(6).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(7).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(8).*d))...
        +(4*fFe+2*fYTm).*exp(1i*(plotdata.Q.*zYTmIG001(9).*d));
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        zTmIG111(1)=0*24;         % 2 Fe
        zTmIG111(2)=0.0024455*24; % 3 O
        zTmIG111(3)=0.0106925*24; % 3 O
        zTmIG111(4)=0.0113862*24; % 3 O
        zTmIG111(5)=0.0196332*24; % 3 O
        zTmIG111(6)=0.0208333*24; % 6 Tm + 6 Fe
        zTmIG111(7)=0.0220335*24; % 3 O
        zTmIG111(8)=0.0302805*24; % 3 O
        zTmIG111(9)=0.0309742*24; % 3 0
        zTmIG111(10)=0.0392212*24;% 3 0
        zTmIG111(11)=0.0416667*24;% 2 Fe
        zYIG111(1)=0*24;         % 2 Fe
        zYIG111(2)=0.0026230*24; % 3 O
        zYIG111(3)=0.0112907*24; % 3 O
        zYIG111(4)=0.0113090*24; % 3 O
        zYIG111(5)=0.0199767*24; % 3 O
        zYIG111(6)=0.0208333*24; % 6 Y + 6 Fe
        zYIG111(7)=0.0216900*24; % 3 O
        zYIG111(8)=0.0303577*24; % 3 O
        zYIG111(9)=0.0303760*24; % 3 0
        zYIG111(10)=0.0390437*24;% 3 0
        zYIG111(11)=0.0416667*24;% 2 Fe
        zYTmIG111=(fitdata.Material(k).xYTmIG*zYIG111+(3-fitdata.Material(k).xYTmIG)*zTmIG111)/3;
        FLayer=2*fFe.*exp(1i*(plotdata.Q.*zYTmIG111(1).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYTmIG111(2).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYTmIG111(3).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYTmIG111(4).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYTmIG111(5).*d))...
        +(6*fYTm+6*fFe).*exp(1i*(plotdata.Q.*zYTmIG111(6).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYTmIG111(7).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYTmIG111(8).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYTmIG111(9).*d))...
        +3*fO.*exp(1i*(plotdata.Q.*zYTmIG111(10).*d))...
        +2*fFe.*exp(1i*(plotdata.Q.*zYTmIG111(11).*d));    
    else
        msgbox('Error with (Y_x,Tm_{3-x})Fe5O12 orientation');
    end;
end;

%% YNiO3
if strcmp(fitdata.Material(k).Type,'YNiO3'),
    %https://materialsproject.org/materials/mvc-15448/
    %a=3.753 b=c=3.759 alpha=89.967 beta=gamma=90 V=53.025
    plotdata.Material(k).V=53.025; %unit cell volume
    fA=AtomicScatteringFactor('Y',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
    if strcmp(fitdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(fitdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with YNiO3 orientation');
    end;
end;

%% ZnMgO
if strcmp(fitdata.Material(k).Type,'(Zn_x,Mg_{1-x})O'),
    % same structure as ZnO but with Mg substituting Zn for Mg<30%
    % ZnO: hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
    % substituting with Mg increases out-of-plane lattice parameter
    a=3.2494; b=a; c=5.20380;
    plotdata.Material(k).V=sqrt(3)*a^2*c/2; %unit cell volume
    fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
    fMg=AtomicScatteringFactor('Mg',plotdata.QQ);
    x=fitdata.Material(k).xZMO;
    if strcmp(fitdata.Material(k).Orientation,'h(0001)'),
        % (0001): z along c-axis
        % inplane surface=sqrt(3)/2*a^2
        % plotdata.Material(k).V(unit cell)=sqrt(3)*a^2*c/2
        zZnMgOh0001(1)=0;     %Zn
        zZnMgOh0001(2)=0.5;   %Zn
        zZnMgOh0001(3)=0.3821;%O
        zZnMgOh0001(4)=0.8821;%O
        FLayer=(x*fZn+(1-x)*fMg).*exp(1i*(plotdata.Q.*zZnMgOh0001(1).*d))...
            +(x*fZn+(1-x)*fMg).*exp(1i*(plotdata.Q.*zZnMgOh0001(2).*d))...
            +fO.*exp(1i*(plotdata.Q.*zZnMgOh0001(3).*d))...
            +fO.*exp(1i*(plotdata.Q.*zZnMgOh0001(4).*d));
    elseif strcmp(fitdata.Material(k).Orientation,'h(11-20)'), 
        %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
        % (11-20): z along bisectrice a,b
        % inplane surface S=sqrt(3)*a*c
        % distance between planes = V(unit cell)/S = a/2
        zZnMgOh11bar20(1)=0;     %2 Zn + 2 O
        %a=3.2494;
        %d=a/2;
        FLayer=(2*(x*fZn+(1-x)*fMg)+2*fO).*exp(1i*(plotdata.Q.*zZnMgOh11bar20(1).*d));
    elseif strcmp(fitdata.Material(k).Orientation,'h(10-10)'),
        %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
        % (10-10): z perpendicular to a
        % inplane surface S=ac
        % distance between planes = V(unit cell)/S = sqrt(3)a/2
        %a=3.2494;
        zZnMgOh10bar10(1)=1/3;   % Zn + O
        zZnMgOh10bar10(2)=2/3;     % Zn + O
        %d=a*sqrt(3)/2;
        FLayer=((x*fZn+(1-x)*fMg)+fO).*exp(1i*(plotdata.Q.*zZnMgOh10bar10(1).*d))+((x*fZn+(1-x)*fMg)+fO).*exp(1i*(plotdata.Q.*zZnMgOh10bar10(2).*d));  
    else
        msgbox('Error with ZnMgO orientation');
    end;
end;

%% Zn3N2
% Ia-3 a=b=c=9.76910 alpha=bet=gamma=90
if strcmp(fitdata.Material(k).Type,'Zn3N2'),
    a=9.76910; b=a; c=a;
    plotdata.Material(k).V=a^3; %unit cell volume
    zZn3N2(1)=0;        %2N
    zZn3N2(2)=0.211013; %2N
    zZn3N2(3)=1.001333; %4Zn
    zZn3N2(4)=1.212345; %4Zn
    zZn3N2(5)=1.463411; %4Zn
    zZn3N2(6)=2.442275; %8N
    zZn3N2(7)=3.421139; %4Zn
    zZn3N2(8)=3.672205; %4Zn
    zZn3N2(9)=3.883217; %4Zn
    zZn3N2(10)=4.673537;%2N
    zZn3N2(11)=4.88455; %4N
    zZn3N2(12)=5.095563;%2N
    zZn3N2(13)=5.885882;%4Zn
    zZn3N2(14)=6.096895;%4Zn
    zZn3N2(15)=6.347961;%4Zn
    zZn3N2(16)=7.326825;%8N
    zZn3N2(17)=8.305689;%4Zn
    zZn3N2(18)=8.556755;%4Zn
    zZn3N2(19)=8.767768;%4Zn
    zZn3N2(20)=9.558087;%2N
    zZn3N2(21)=9.7691;  %2N
    zZn3N2=zZn3N2/c; %to go from absolute to relative units
    fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
    fN=AtomicScatteringFactor('N',plotdata.QQ);
    FLayer=2*fN.*exp(1i*(plotdata.Q.*zZn3N2(1).*d))...
        +2*fN.*exp(1i*(plotdata.Q.*zZn3N2(2).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(3).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(4).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(5).*d))...
        +8*fN.*exp(1i*(plotdata.Q.*zZn3N2(6).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(7).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(8).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(9).*d))...
        +2*fN.*exp(1i*(plotdata.Q.*zZn3N2(10).*d))...
        +4*fN.*exp(1i*(plotdata.Q.*zZn3N2(11).*d))...
        +2*fN.*exp(1i*(plotdata.Q.*zZn3N2(12).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(13).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(14).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(15).*d))...
        +8*fN.*exp(1i*(plotdata.Q.*zZn3N2(16).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(17).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(18).*d))...
        +4*fZn.*exp(1i*(plotdata.Q.*zZn3N2(19).*d))...
        +2*fN.*exp(1i*(plotdata.Q.*zZn3N2(20).*d))...
        +2*fN.*exp(1i*(plotdata.Q.*zZn3N2(21).*d));
end;

%% ZnO
if strcmp(fitdata.Material(k).Type,'ZnO'), 
    %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
    a=3.2494; b=a; c=5.20380;
    plotdata.Material(k).V=sqrt(3)*a^2*c/2; %unit cell volume
    if strcmp(fitdata.Material(k).Orientation,'h(0001)'),
        % (0001): z along c-axis
        % inplane surface=sqrt(3)/2*a^2
        % V(unit cell)=sqrt(3)*a^2*c/2
        zZnOh0001(1)=0;     %Zn
        zZnOh0001(2)=0.5;   %Zn
        zZnOh0001(3)=0.3821;%O
        zZnOh0001(4)=0.8821;%O
        fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
        FLayer=fZn.*exp(1i*(plotdata.Q.*zZnOh0001(1).*d))...
            +fZn.*exp(1i*(plotdata.Q.*zZnOh0001(2).*d))...
            +fO.*exp(1i*(plotdata.Q.*zZnOh0001(3).*d))...
            +fO.*exp(1i*(plotdata.Q.*zZnOh0001(4).*d));
    elseif strcmp(fitdata.Material(k).Orientation,'h(11-20)'), 
        %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
        % (11-20): z along bisectrice (a,b)
        % inplane surface S=sqrt(3)*a*c
        % distance between planes = V(unit cell)/S = a/2
        zZnOh11bar20(1)=0;     %2 Zn + 2 O
        %a=3.2494;
        %d=a/2;
        fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
        FLayer=(2*fZn+2*fO).*exp(1i*(plotdata.Q.*zZnOh11bar20(1).*d));
    elseif strcmp(fitdata.Material(k).Orientation,'h(10-10)'),
        %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
        % (10-10): z perpendicular to a
        % inplane surface S=ac
        % distance between planes = V(unit cell)/S = sqrt(3)a/2
        % a=3.2494;
        zZnOh10bar10(1)=1/3;   % Zn + O
        zZnOh10bar10(2)=2/3;     % Zn + O
        %d=a*sqrt(3)/2;
        fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
        FLayer=(fZn+fO).*exp(1i*(plotdata.Q.*zZnOh10bar10(1).*d))+(fZn+fO).*exp(1i*(plotdata.Q.*zZnOh10bar10(2).*d));  
    else
        msgbox('Error with ZnO orientation');
    end;
end;

%% ZrO2
if strcmp(fitdata.Material(k).Type,'ZrO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('Zr',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;
FLayer=FLayer*d/plotdata.Material(k).V; %renormalized by the in-plane surface of the unit cell