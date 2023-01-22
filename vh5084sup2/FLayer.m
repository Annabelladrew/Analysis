% Created by Celine Lichtensteiger
% Calculates the form factors of the different materials depending on their
% crystaollographic structures
% "z" are the out-of-plane atomic positions in the unit cells, expressed 
% in relative unit
%*********************************
function[FLayer]=FLayer(k,c)
global plotdata;

%% Perovskite (zperovskite) and double perovskite (zdbpervoskite):
if strcmp(plotdata.orientation,'(001)'), %Glazer-Acta Crystallogr B-1978.pdf

zperovskite(1)=0;                                                   %A
zperovskite(2)=0.5+plotdata.Material(k).Polarization*0.162/4.156;   %B
zperovskite(3)=0+plotdata.Material(k).Polarization*0.473/4.156;     %O1
zperovskite(4)=0.5+plotdata.Material(k).Polarization*0.486/4.156;	%O2
zperovskite(5)=0.5+plotdata.Material(k).Polarization*0.486/4.156;	%O3

elseif strcmp(plotdata.orientation,'(111)'),

zperovskite(1)=0;   %A
zperovskite(2)=0.5; %B
zperovskite(3)=0;	%O1
zperovskite(4)=0;	%O2
zperovskite(5)=0;	%O3

zdbperovskite(1)=0;   %A
zdbperovskite(2)=0.25;%B-1
zdbperovskite(3)=0;   %O1
zdbperovskite(4)=0;	  %O2
zdbperovskite(5)=0;	  %O3
zdbperovskite(6)=0.5; %A
zdbperovskite(7)=0.75;%B-2
zdbperovskite(8)=0.5; %O1
zdbperovskite(9)=0.5; %O2
zdbperovskite(10)=0.5;%O3

else
    warning('Error in FLayer (orientation)');
end

%% Perovskite

fA=zeros(size(plotdata.QQ));
fB=zeros(size(plotdata.QQ));

if strcmp(plotdata.Material(k).Type,'BaTiO3'),
    fA=AtomicScatteringFactor('Ba',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'BiFeO3'),
    fA=AtomicScatteringFactor('Bi',plotdata.QQ); 
    fB=AtomicScatteringFactor('Fe',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'LaAlO3'),
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Al',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'LaFeO3'),
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Fe',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'LaMnO3'),
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Mn',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'LaNiO3'),
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'LSMO'),
    fA=0.67*AtomicScatteringFactor('La',plotdata.QQ)+0.33*AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Mn',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'MnTiO3'),
    fA=AtomicScatteringFactor('Mn',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'NdNiO3'),
    fA=AtomicScatteringFactor('Nd',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'PbTiO3'),
    fA=AtomicScatteringFactor('Pb',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'(Pb_x,Sr_{1-x})TiO3'),
    fA=plotdata.Material(k).xPST*AtomicScatteringFactor('Pb',plotdata.QQ)+(1-plotdata.Material(k).xPST)*AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'Pb(Zr_x,Ti_{1-x})O3'),
    fA=AtomicScatteringFactor('Pb',plotdata.QQ); 
    fB=plotdata.Material(k).xPZT*AtomicScatteringFactor('Zr',plotdata.QQ)+(1-plotdata.Material(k).xPZT)*AtomicScatteringFactor('Ti',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'SmNiO3'),
    fA=AtomicScatteringFactor('Sm',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'SrRuO3'),
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ru',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'SrTiO3'),
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
end;
if strcmp(plotdata.Material(k).Type,'SrVO3'),
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('V',plotdata.QQ);
end;

fO=AtomicScatteringFactor('O',plotdata.QQ);
FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite,c);

%% Different structures

if strcmp(plotdata.Material(k).Type,'none'),
    fA=zeros(size(plotdata.QQ));
    fB=zeros(size(plotdata.QQ));
    fO=zeros(size(plotdata.QQ));
    FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite,c);
end;

%% PrBa2Cu3O7 %MODIFY TO USE THE CORRECT ATOMIC POSITIONS
if strcmp(plotdata.Material(k).Type,'PrBa2Cu3O7'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate PrBa2Cu3O7-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
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
    fPr=AtomicScatteringFactor('Pr',plotdata.QQ);
    fBa=AtomicScatteringFactor('Ba',plotdata.QQ);
    fCu=AtomicScatteringFactor('Cu',plotdata.QQ);
    FLayer=fPr.*exp(1i*(plotdata.Q.*zYBCO(1).*c))...
        +fBa.*exp(1i*(plotdata.Q.*zYBCO(2).*c))...
        +fBa.*exp(1i*(plotdata.Q.*zYBCO(3).*c))...
        +fCu.*exp(1i*(plotdata.Q.*zYBCO(4).*c))...
        +fCu.*exp(1i*(plotdata.Q.*zYBCO(5).*c))...
        +fCu.*exp(1i*(plotdata.Q.*zYBCO(6).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(7).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(8).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(9).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(10).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(11).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(12).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(13).*c));
     end;
end;
%% YBCO:
if strcmp(plotdata.Material(k).Type,'YBa2Cu3O7'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate YBa2Cu3O7-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
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
    FLayer=fY.*exp(1i*(plotdata.Q.*zYBCO(1).*c))...
        +fBa.*exp(1i*(plotdata.Q.*zYBCO(2).*c))...
        +fBa.*exp(1i*(plotdata.Q.*zYBCO(3).*c))...
        +fCu.*exp(1i*(plotdata.Q.*zYBCO(4).*c))...
        +fCu.*exp(1i*(plotdata.Q.*zYBCO(5).*c))...
        +fCu.*exp(1i*(plotdata.Q.*zYBCO(6).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(7).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(8).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(9).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(10).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(11).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(12).*c))...
        +fO.*exp(1i*(plotdata.Q.*zYBCO(13).*c));
     end;
end;

%% AO:
zAO(1)=0; %A
zAO(2)=0; %O

if strcmp(plotdata.Material(k).Type,'BaO'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate BaO-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
         fA=AtomicScatteringFactor('Ba',plotdata.QQ);
         FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*c))+fO.*exp(1i*(plotdata.Q.*zAO(2).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'LaO'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate LaO-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
         fA=AtomicScatteringFactor('La',plotdata.QQ);
         FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*c))+fO.*exp(1i*(plotdata.Q.*zAO(2).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'MnO'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate MnO-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
         fA=AtomicScatteringFactor('Mn',plotdata.QQ);
         FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*c))+fO.*exp(1i*(plotdata.Q.*zAO(2).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'NdO'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate NdO-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
         fA=AtomicScatteringFactor('Nd',plotdata.QQ);
         FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*c))+fO.*exp(1i*(plotdata.Q.*zAO(2).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'PbO'),
         if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate PbO-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
    fA=AtomicScatteringFactor('Pb',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*c))+fO.*exp(1i*(plotdata.Q.*zAO(2).*c));
         end;
end;

if strcmp(plotdata.Material(k).Type,'SrO'),
         if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate SrO-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
    fA=AtomicScatteringFactor('Sr',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*c))+fO.*exp(1i*(plotdata.Q.*zAO(2).*c));
         end;
end;

%% BO2:
zBO2(1)=0; %B
zBO2(2)=0; %O
zBO2(3)=0; %O
if strcmp(plotdata.Material(k).Type,'AlO2'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate AlO2-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
         fB=AtomicScatteringFactor('Al',plotdata.QQ);
         FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'MnO2'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate MnO2-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
         fB=AtomicScatteringFactor('Mn',plotdata.QQ);
         FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'NiO2'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate NiO2-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
         fB=AtomicScatteringFactor('Ni',plotdata.QQ);
         FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'RuO2'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate NiO2-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
         fB=AtomicScatteringFactor('Ru',plotdata.QQ);
         FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'SrO2'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate SrO2-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
         fB=AtomicScatteringFactor('Sr',plotdata.QQ);
         FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'TiO2'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate TiO2-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
         fB=AtomicScatteringFactor('Ti',plotdata.QQ);
         FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'VO2'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate VO2-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else fB=AtomicScatteringFactor('V',plotdata.QQ);
         FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*c));
     end;
end;

if strcmp(plotdata.Material(k).Type,'ZrO2'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate ZrO2-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else fB=AtomicScatteringFactor('Zr',plotdata.QQ);
         FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*c))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*c));
     end;
     end;
%% La2CuO4 (CIF orthorhombic Cmca a=5.35 b=13.148 c=5.398):

if strcmp(plotdata.Material(k).Type,'La2CuO4'),
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
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate La2CuO4-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
    fLa=AtomicScatteringFactor('La',plotdata.QQ); 
    fCu=AtomicScatteringFactor('Cu',plotdata.QQ);
    FLayer=fCu.*exp(1i*(plotdata.Q.*zLa2CuO4(1).*c))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(2).*c))...
        +fLa.*exp(1i*(plotdata.Q.*zLa2CuO4(3).*c))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(4).*c))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(5).*c))...
        +fLa.*exp(1i*(plotdata.Q.*zLa2CuO4(6).*c))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(7).*c))...
        +fCu.*exp(1i*(plotdata.Q.*zLa2CuO4(8).*c))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(9).*c))...
        +fLa.*exp(1i*(plotdata.Q.*zLa2CuO4(10).*c))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(11).*c))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(12).*c))...
        +fLa.*exp(1i*(plotdata.Q.*zLa2CuO4(13).*c))...
        +fO.*exp(1i*(plotdata.Q.*zLa2CuO4(14).*c));
     end;
end;
%% CaCuO2 (CIF CaCuO2 tetragonal a=b=3.87328 c=3.20546):
if strcmp(plotdata.Material(k).Type,'CaCuO2'),
     if strcmp(plotdata.orientation,'(111)'),
         msgbox('This program doesn''t allow you to calculate CaCuO2-(111) oriented yet. If interested, please contact Celine.Lichtensteiger@unige.ch.');
         FLayer=StructureFactor(zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),zeros(size(plotdata.QQ)),plotdata.Q,zperovskite,c);
     else
    zCaCuO2(1)=0.5;                     %Ca
    zCaCuO2(2)=0;                       %Cu
    zCaCuO2(3)=0;                       %O
    zCaCuO2(4)=0;                       %O
    fCa=AtomicScatteringFactor('Ca',plotdata.QQ); 
    fCu=AtomicScatteringFactor('Cu',plotdata.QQ);
    FLayer=fCa.*exp(1i*(plotdata.Q.*zCaCuO2(1).*c))...
        +fCu.*exp(1i*(plotdata.Q.*zCaCuO2(2).*c))...
        +fO.*exp(1i*(plotdata.Q.*zCaCuO2(3).*c))...
        +fO.*exp(1i*(plotdata.Q.*zCaCuO2(4).*c));
     end;
end;
%% DoublePerovskites
if strcmp(plotdata.Material(k).Type,'La2NiMnO6'),
    if strcmp(plotdata.orientation,'(111)'),
        fLa=AtomicScatteringFactor('La',plotdata.QQ);
        fNi=AtomicScatteringFactor('Ni',plotdata.QQ); 
        fMn=AtomicScatteringFactor('Mn',plotdata.QQ);
        fO=AtomicScatteringFactor('O',plotdata.QQ);
        FLayer=fLa.*exp(1i*(plotdata.Q.*zdbperovskite(1).*c))...
        +fNi.*exp(1i*(plotdata.Q.*zdbperovskite(2).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(3).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(4).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(5).*c))...
        +fLa.*exp(1i*(plotdata.Q.*zdbperovskite(6).*c))...
        +fMn.*exp(1i*(plotdata.Q.*zdbperovskite(7).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(8).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(9).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(10).*c));
       else 
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=0.5*AtomicScatteringFactor('Ni',plotdata.QQ)+0.5*AtomicScatteringFactor('Mn',plotdata.QQ);
    fO=AtomicScatteringFactor('O',plotdata.QQ);
    FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite,c);
      end;
end;
if strcmp(plotdata.Material(k).Type,'Nd2NiMnO6'),
    if strcmp(plotdata.orientation,'(111)'),
        fNd=AtomicScatteringFactor('Nd',plotdata.QQ);
        fNi=AtomicScatteringFactor('Ni',plotdata.QQ); 
        fMn=AtomicScatteringFactor('Mn',plotdata.QQ);
        fO=AtomicScatteringFactor('O',plotdata.QQ);
        FLayer=fNd.*exp(1i*(plotdata.Q.*zdbperovskite(1).*c))...
        +fNi.*exp(1i*(plotdata.Q.*zdbperovskite(2).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(3).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(4).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(5).*c))...
        +fNd.*exp(1i*(plotdata.Q.*zdbperovskite(6).*c))...
        +fMn.*exp(1i*(plotdata.Q.*zdbperovskite(7).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(8).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(9).*c))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite(10).*c));
       else 
    fA=AtomicScatteringFactor('Nd',plotdata.QQ); 
    fB=0.5*AtomicScatteringFactor('Ni',plotdata.QQ)+0.5*AtomicScatteringFactor('Mn',plotdata.QQ);
    fO=AtomicScatteringFactor('O',plotdata.QQ);
    FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite,c);
      end;
end;