% Created by Celine Lichtensteiger
% Calculates the form factors of the different materials depending on their
% crystaollographic structures
% "z" are the out-of-plane atomic positions in the unit cells, expressed 
% in relative unit
%*********************************
function[FLayer]=FLayer(k,d)
global plotdata;
z=[]; f=[];
fA=zeros(size(plotdata.QQ));
fB=zeros(size(plotdata.QQ));
fO=AtomicScatteringFactor('O',plotdata.QQ);
FLayer=fA.*exp(1i*(plotdata.Q.*0));

%% Perovskite
%Glazer-Acta Crystallogr B-1978.pdf

zperovskite001(1)=0;                                                   %A
zperovskite001(2)=0.5+plotdata.Material(k).Polarization*0.162/4.156;    %B
zperovskite001(3)=0+plotdata.Material(k).Polarization*0.473/4.156;      %O1
zperovskite001(4)=0.5+plotdata.Material(k).Polarization*0.486/4.156;    %O2
zperovskite001(5)=0.5+plotdata.Material(k).Polarization*0.486/4.156;    %O3

zperovskite111(1)=0;   %A
zperovskite111(2)=0.5; %B
zperovskite111(3)=0;	%O1
zperovskite111(4)=0;	%O2
zperovskite111(5)=0;	%O3

%% Double perovskite
zdbperovskite111(1)=0;   %A-1
zdbperovskite111(2)=0.25;%B-1
zdbperovskite111(3)=0;   %O1
zdbperovskite111(4)=0;	  %O2
zdbperovskite111(5)=0;	  %O3
zdbperovskite111(6)=0.5; %A-2
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

%% none
if strcmp(plotdata.Material(k).Type,'none'),
    fO=zeros(size(plotdata.QQ));
    FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
end;

%% AlO2
if strcmp(plotdata.Material(k).Type,'AlO2'), %%BO2-type monolayer
       %plotdata.Vunitcell=54.052; % for Al2O4 https://materialsproject.org/materials/mp-1182858/
       %plotdata.Material(k).Sunitcell=plotdata.Material(k).Vunitcell/d;
       plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
       fB=AtomicScatteringFactor('Al',plotdata.QQ);
       FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
       %FLayer=Flayer/plotdata.Material(k).Sunitcell; %renormalized by the surface size to take into account the relative atomic density
       %if strained, the Sunitcell should be the same as the one of the
       %substrate
end;

%% BaBiO3
if strcmp(plotdata.Material(k).Type,'BaBiO3'),
    plotdata.Material(k).V=4.43647^3; %unit cell volume https://materialsproject.org/materials/mp-1227955/
    fA=AtomicScatteringFactor('Ba',plotdata.QQ); 
    fB=AtomicScatteringFactor('Bi',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with BaBiO3 orientation');
    end;
end;

%% BaO
if strcmp(plotdata.Material(k).Type,'BaO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('Ba',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% BaSnO3
if strcmp(plotdata.Material(k).Type,'BaSnO3'),
    plotdata.Material(k).V=4.188634^3; %unit cell volume https://materialsproject.org/materials/mp-3163/
    fA=AtomicScatteringFactor('Ba',plotdata.QQ); 
    fB=AtomicScatteringFactor('Sn',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with BaSnO3 orientation');
    end;
end;

%% (Ba_x,Sr_{1-x})TiO3
if strcmp(plotdata.Material(k).Type,'(Ba_x,Sr_{1-x})TiO3'),
    %BaTiO3: cubic a=4.036 https://materialsproject.org/materials/mp-2998/
    %SrTiO3: cubic a=3.945 https://materialsproject.org/materials/mp-5229/    
    plotdata.Material(k).V=plotdata.Material(k).xBST*4.036^3+(1-plotdata.Material(k).xBST)*3.945^3;%volume of the unit cell
    fA=plotdata.Material(k).xBST*AtomicScatteringFactor('Ba',plotdata.QQ)+(1-plotdata.Material(k).xBST)*AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
    %Sunitcell=Vunitcell/plotdata.Material(k).d; %in-plane area of the unit cell
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
        %FLayer=Flayer/Sunitcell; %renormalized by the surface size to take into account the relative atomic density
        %if strained, the Sunitcell should be the same as the one of the
        %substrate
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
        %FLayer=Flayer/Sunitcell; %renormalized by the surface size to take into account the relative atomic density
        %if strained, the Sunitcell should be the same as the one of the
        %substrate
    else
        msgbox('Error with (Ba_x,Sr_{1-x})TiO3 orientation');
    end;
end;

%% BaTiO3
if strcmp(plotdata.Material(k).Type,'BaTiO3'),
    plotdata.Material(k).V=4.036^3; %unit cell volume https://materialsproject.org/materials/mp-2998/
    fA=AtomicScatteringFactor('Ba',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with BaTiO3 orientation');
    end;
end;

%% BiFeO3
if strcmp(plotdata.Material(k).Type,'BiFeO3'),
    plotdata.Material(k).V=3.753^2*4.924; %unit cell volume https://materialsproject.org/materials/mp-1069079/
    fA=AtomicScatteringFactor('Bi',plotdata.QQ); 
    fB=AtomicScatteringFactor('Fe',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with BiFeO3 orientation');
    end;
end;

%% Bi2Te3
if strcmp(plotdata.Material(k).Type,'Bi2Te3'),
    % Bi2Te3: hexagonal with primitive unit cell Bi6Te9 a=b=4.44949, c=31.150873, alpha=beta=90?, gamma=120?
    a=4.44949; b=a; c=31.150873;
    plotdata.Material(k).V=sqrt(3)*a^2*c/2; %unit cell volume
    fBi=AtomicScatteringFactor('Bi',plotdata.QQ);
    fTe=AtomicScatteringFactor('Te',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'h(0001)'),
        % (0001): z along c-axis
        % inplane surface=sqrt(3)/2*a^2
        z(1)=0; f(1,:)=0.5*fTe; 
        z(2)=0.0659236741134029; f(2,:)=fBi;
        z(3)=0.121489339961676; f(3,:)=fTe;
        z(4)=0.211843982671047; f(4,:)=fTe;
        z(5)=0.26740968062115; f(5,:)=fBi;
        z(6)=0.333333354734553; f(6,:)=fTe;
        z(7)=0.399256996746127; f(7,:)=fBi;
        z(8)=0.4548226625944; f(8,:)=fTe;
        z(9)=0.5451773374056; f(9,:)=fTe;
        z(10)=0.600743003253873; f(10,:)=fBi;
        z(11)=0.666666709469105; f(11,:)=fTe;
        z(12)=0.732590287277021; f(12,:)=fBi;
        z(13)=0.788155953125294; f(13,:)=fTe;
        z(14)=0.878510627936495; f(14,:)=fTe;
        z(15)=0.934076293784768; f(15,:)=fBi;
        z(16)=1; f(16,:)=0.5*fTe;
        for j=1:16,
            FLayer=FLayer+f(j,:).*exp(1i*(plotdata.Q.*z(j).*d));
        end;
    elseif strcmp(plotdata.Material(k).Orientation,'h(11-20)'), 
        % (11-20): z along bisectrice a,b
        % inplane surface S=sqrt(3)*a*c
        % distance between planes = V(unit cell)/S = a/2
        % d=a/2;
        % TODO 
    dlg = dialog('Position',[300 300 250 300],'Name','Warning');
    txt = uicontrol('Parent',dlg,...
               'Style','text',...
               'Position',[20 80 210 200],...
               'String',['The ' plotdata.Material(k).Orientation ' orientation for ' plotdata.Material(k).Type ' is not available yet. Please contact Celine.Lichtensteiger@unige.ch if needed.']);
    btn = uicontrol('Parent',dlg,...
               'Position',[85 20 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');    
    plotdata.Material(k).Orientation='h(0001)';
    plotdata.Materuak(k).plotdata.Material(k).choosehOrientation.Value=1;
    elseif strcmp(plotdata.Material(k).Orientation,'h(10-10)'),
        % (10-10): z perpendicular to a
        % inplane surface S=ac
        % distance between planes = V(unit cell)/S = sqrt(3)a/2
        % d=a*sqrt(3)/2;
        % TODO;
    dlg = dialog('Position',[300 300 250 300],'Name','Warning');
    txt = uicontrol('Parent',dlg,...
               'Style','text',...
               'Position',[20 80 210 200],...
               'String',['The ' plotdata.Material(k).Orientation ' orientation for ' plotdata.Material(k).Type ' is not available yet. Please contact Celine.Lichtensteiger@unige.ch if needed.']);
    btn = uicontrol('Parent',dlg,...
               'Position',[85 20 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');    
    plotdata.Material(k).Orientation='h(0001)';
    plotdata.Materuak(k).plotdata.Material(k).choosehOrientation.Value=1;
    else
        msgbox('Error with Bi2Te3 orientation');
    end;
end;

%% CaCuO2 (CIF CaCuO2 tetragonal a=b=3.87328 c=3.20546):
if strcmp(plotdata.Material(k).Type,'CaCuO2'),
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
if strcmp(plotdata.Material(k).Type,'Ca2RuO4'),
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

%% (Ca_{3-x},Sr_x)Al2O6
if strcmp(plotdata.Material(k).Type,'(Ca_{3-x},Sr_x)Al2O6'),
    
    %% Using relative atomic positions from Sr3Al2O6 https://materialsproject.org/materials/mp-3393/ 
    %Sr is replaced by (Ca_{3-x},Sr_x)
    %Sr72Al48O144 in cubic unit cell a=b=c=15.999914A
    %In the (001) orientation, in-plane 16 unit cells of STO will fit in 1 unit cell of Sr72Al48O144. This is taken into account by divinding by 16 the StructureFactor.
    %Out-of-plane, I just consider half of the Sr72Al48O144 unit cell, i.e. with a bulk lattice parameter of 7.999957A (that the user can then modify). 
    a=15.999914;
    plotdata.Material(k).V=a^3/2; 
    zCSAO(1)=0;                      %2 (Ca_{3-x}Sr_x)/3
    zCSAO(2)=0.00283001521133176*2;  %2 O
    zCSAO(3)=0.0041260221773692*2;   %2 Al
    zCSAO(4)=0.0076369785487597*2;   %2 O
    zCSAO(5)=0.0159040229841235*2;   %2 Al
    zCSAO(6)=0.0168720281871515*2;   %2 O
    zCSAO(7)=0.0175839695138361*2;   %2 O
    zCSAO(8)=0.0199369821612791*2;   %2 O
    zCSAO(9)=0.0200889829782835*2;   %2 Al
    zCSAO(10)=0.0260960152660821*2;  %2 O
    zCSAO(11)=0.094133005964907*2;   %2 O
    zCSAO(12)=0.100253976365123*2;   %2 O
    zCSAO(13)=0.115953998252741*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(14)=0.119861019252978*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(15)=0.124367980977898*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(16)=0.124571982074404*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(17)=0.124742982993534*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(18)=0.125220985562797*2;   %2 O
    zCSAO(19)=0.129366007842292*2;   %2 O
    zCSAO(20)=0.131721020500485*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(21)=0.134670036351445*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(22)=0.153380011917564*2;   %2 O
    zCSAO(23)=0.153955015008206*2;   %2 O
    zCSAO(24)=0.220413997225235*2;   %2 O
    zCSAO(25)=0.222892010544557*2;   %2 O
    zCSAO(26)=0.230195987303432*2;   %2 O
    zCSAO(27)=0.233955007508165*2;   %2 Al
    zCSAO(28)=0.234410009953803*2;   %2 O
    zCSAO(29)=0.236916023423626*2;   %2 O
    zCSAO(30)=0.238600032475175*2;   %2 Al
    zCSAO(31)=0.244254000365252*2;   %2 O
    zCSAO(32)=0.247728019038102*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(33)=0.248376022521121*2;   %2 Al
    zCSAO(34)=0.251623977478879*2;   %2 Al
    zCSAO(35)=0.252272043462234*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(36)=0.255745999634748*2;   %2 O
    zCSAO(37)=0.261400030025161*2;   %2 Al
    zCSAO(38)=0.263083976576374*2;   %2 O
    zCSAO(39)=0.265589990046196*2;   %2 O
    zCSAO(40)=0.266044992491835*2;   %2 Al
    zCSAO(41)=0.269804012696568*2;   %2 O
    zCSAO(42)=0.277107989455443*2;   %2 O
    zCSAO(43)=0.279586002774765*2;   %2 O
    zCSAO(44)=0.346044984991794*2;   %2 O
    zCSAO(45)=0.346619988082436*2;   %2 O
    zCSAO(46)=0.365330026148891*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(47)=0.368279041999851*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(48)=0.370633992157708*2;   %2 O
    zCSAO(49)=0.374779014437203*2;   %2 O
    zCSAO(50)=0.37525695450613*2;    %2 (Ca_{3-x}Sr_x)/3
    zCSAO(51)=0.37542795542526*2;    %2 (Ca_{3-x}Sr_x)/3
    zCSAO(52)=0.375631956521766*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(53)=0.380138980747022*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(54)=0.384046001747259*2;   %2 (Ca_{3-x}Sr_x)/3
    zCSAO(55)=0.399744961129166*2;   %2 O
    zCSAO(56)=0.405866994035093*2;   %2 O
    zCSAO(57)=0.473904047234254*2;   %2 O
    zCSAO(58)=0.479911017021717*2;   %2 Al
    zCSAO(59)=0.480063017838721*2;   %2 O
    zCSAO(60)=0.482416030486164*2;   %2 O
    zCSAO(61)=0.483128034313184*2;   %2 O
    zCSAO(62)=0.484095977015876*2;   %2 Al
    zCSAO(63)=0.49236302145124*2;    %2 O
    zCSAO(64)=0.495873977822631*2;   %2 Al
    zCSAO(65)=0.497170047289004*2;   %2 O
    zCSAO(66)=0.5*2;                 %2 (Ca_{3-x}Sr_x)/3
    
    fCaSr=plotdata.Material(k).xCSA/3*AtomicScatteringFactor('Sr',plotdata.QQ)+(3-plotdata.Material(k).xCSA)/3*AtomicScatteringFactor('Ca',plotdata.QQ);
    fAl=AtomicScatteringFactor('Al',plotdata.QQ); 
    fO=AtomicScatteringFactor('O',plotdata.QQ);
    
    FLayer=+2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(1).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(2).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(3).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(4).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(5).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(6).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(7).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(8).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(9).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(10).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(11).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(12).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(13).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(14).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(15).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(16).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(17).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(18).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(19).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(20).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(21).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(22).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(23).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(24).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(25).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(26).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(27).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(28).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(29).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(30).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(31).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(32).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(33).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(34).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(35).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(36).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(37).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(38).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(39).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(40).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(41).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(42).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(43).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(44).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(45).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(46).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(47).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(48).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(49).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(50).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(51).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(52).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(53).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(54).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(55).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(56).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(57).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(58).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(59).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(60).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(61).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(62).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(63).*d))...
        +2*fAl.*exp(1i*(plotdata.Q.*zCSAO(64).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zCSAO(65).*d))...
        +2*fCaSr.*exp(1i*(plotdata.Q.*zCSAO(66).*d));
end;

%% CaTiO3
if strcmp(plotdata.Material(k).Type,'CaTiO3'),
    plotdata.Material(k).V=3.889^3; %unit cell volume https://materialsproject.org/materials/mp-5827/
    fA=AtomicScatteringFactor('Ca',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with CaTiO3 orientation');
    end;
end;

%% CaVO3
if strcmp(plotdata.Material(k).Type,'CaVO3'),
    plotdata.Material(k).V=3.830^3; %unit cell volume https://materialsproject.org/materials/mp-1016853/
    fA=AtomicScatteringFactor('Ca',plotdata.QQ); 
    fB=AtomicScatteringFactor('V',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with CaVO3 orientation');
    end;
end;

%% CoO (CIF CoO cubic a=b=d=4.263 all angles=90?) - idem MgO
if strcmp(plotdata.Material(k).Type,'CoO'),
    plotdata.Material(k).V=4.263^3; 
    fCo=AtomicScatteringFactor('Co',plotdata.QQ); 
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        zCoOc001(1)=0;                      % Co-O
        zCoOc001(2)=0.5;                    % 2Co-2O
        zCoOc001(3)=1;                      % CoO
        FLayer=(fCo+fO).*exp(1i*(plotdata.Q.*zCoOc001(1).*d))...
        +2*(fCo+fO).*exp(1i*(plotdata.Q.*zCoOc001(2).*d))...
        +(fCo+fO).*exp(1i*(plotdata.Q.*zCoOc001(3).*d));  %Structure factor F 
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
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

%% mu-Fe2O3
 if strcmp(plotdata.Material(k).Type,'mu-Fe2O3'),
        % 8xFe + 12xO
        % growth on (010) planes: a,c in-plane, b out-of-plane
        a=12.42035; b=3.07304; c=5.87522; alpha=90; beta=103.94; gamma=90;
        % d=b;                          %out-of-plane lattice parameter
        % S=a*c*sin(beta*pi/180);               %in-plane surface
        plotdata.Material(k).V=a*b*c*sin(beta*pi/180);            %volume of the unit cell
        fFe=AtomicScatteringFactor('Fe',plotdata.QQ);
        zmuFe2O3(1)=0;       % 3 O + 2 Fe
        zmuFe2O3(2)=0.5;     % 6 O + 4 Fe
        zmuFe2O3(3)=1;       % 3 O + 2 Fe
        FLayer=(3*fO+2*fFe).*exp(1i*(plotdata.Q.*zmuFe2O3(1).*d))...
                +(6*fO+4*fFe).*exp(1i*(plotdata.Q.*zmuFe2O3(2).*d))...
                +(3*fO+2*fFe).*exp(1i*(plotdata.Q.*zmuFe2O3(3).*d));
 end;

%% beta-Ga2O3
 if strcmp(plotdata.Material(k).Type,'beta-Ga2O3'),
        % 8xGa + 11xO +8xC
        % growth on (010) planes: a,c in-plane, b out-of-plane
        a=12.214; b=3.0371; c=5.7981; alpha=90; beta=103.83; gamma=90;
        % d=b;                          %out-of-plane lattice parameter
        % S=a*c*sin(beta*pi/180);               %in-plane surface
        plotdata.Material(k).V=a*b*c*sin(beta*pi/180);            %volume of the unit cell
        fGa=AtomicScatteringFactor('Ga',plotdata.QQ);
        fC=AtomicScatteringFactor('C',plotdata.QQ);
        zbetaGa2O3(1)=0;       % 3 O + 2 Ga
        zbetaGa2O3(2)=0.141618;% 2 C
        zbetaGa2O3(3)=0.358382;% 2 C
        zbetaGa2O3(4)=0.5;     % 5 O + 4 Ga
        zbetaGa2O3(5)=0.641618;% 2 C
        zbetaGa2O3(6)=0.858382;% 2 C
        zbetaGa2O3(7)=1;       % 3 O + 2 Ga
        FLayer=(3*fO+2*fGa).*exp(1i*(plotdata.Q.*zbetaGa2O3(1).*d))...
                +2*fC.*exp(1i*(plotdata.Q.*zbetaGa2O3(2).*d))...
                +2*fC.*exp(1i*(plotdata.Q.*zbetaGa2O3(3).*d))...
                +(5*fO+4*fGa).*exp(1i*(plotdata.Q.*zbetaGa2O3(4).*d))...
                +2*fC.*exp(1i*(plotdata.Q.*zbetaGa2O3(5).*d))...
                +2*fC.*exp(1i*(plotdata.Q.*zbetaGa2O3(6).*d))...
                +(3*fO+2*fGa).*exp(1i*(plotdata.Q.*zbetaGa2O3(7).*d));
 end;
        
 
%% GeTe
if strcmp(plotdata.Material(k).Type,'GeTe'),
    % GeTe: hexagonal with primitive unit cell Ge3Te3 a=b= 4.23067, c= 10.88958, alpha=beta=90, gamma=120
    a=4.23067; b=a; c=10.88958;
    plotdata.Material(k).V=sqrt(3)*a^2*c/2; %unit cell volume
    fGe=AtomicScatteringFactor('Ge',plotdata.QQ);
    fTe=AtomicScatteringFactor('Te',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'h(0001)'),
        % (0001): z along c-axis
        % inplane surface=sqrt(3)/2*a^2
        z(1)=0.00494004360131428; f(1,:)=fGe; 
        z(2)=0.14072700691854; f(2,:)=fTe;
        z(3)=0.338273009611023; f(3,:)=fGe;
        z(4)=0.474059972928249; f(4,:)=fTe;
        z(5)=0.671605975620731; f(5,:)=fGe;
        z(6)=0.807394040908832; f(6,:)=fTe;
        for j=1:6,
            FLayer=FLayer+f(j,:).*exp(1i*(plotdata.Q.*z(j).*d));
        end;
    elseif or(strcmp(plotdata.Material(k).Orientation,'h(11-20)'),strcmp(plotdata.Material(k).Orientation,'h(10-10)')), 
        % (11-20): z along bisectrice a,b
        % inplane surface S=sqrt(3)*a*c
        % distance between planes = V(unit cell)/S = a/2
        % d=a/2;
        % TODO 
        % (10-10): z perpendicular to a
        % inplane surface S=ac
        % distance between planes = V(unit cell)/S = sqrt(3)a/2
        % d=a*sqrt(3)/2;
        % TODO;
    dlg = dialog('Position',[300 300 250 300],'Name','Warning');
    txt = uicontrol('Parent',dlg,...
               'Style','text',...
               'Position',[20 80 210 200],...
               'String',['The ' plotdata.Material(k).Orientation ' orientation for ' plotdata.Material(k).Type ' is not available yet. Please contact Celine.Lichtensteiger@unige.ch if needed.']);
    btn = uicontrol('Parent',dlg,...
               'Position',[85 20 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');    
    plotdata.Material(k).Orientation='h(0001)';
    plotdata.Materuak(k).plotdata.Material(k).choosehOrientation.Value=1;
    else
        msgbox('Error with GeTe orientation');
    end;
end;

%% (Hf_{1-x},Zr_x)O2
if strcmp(plotdata.Material(k).Type,'(Hf_{1-x},Zr_x)O2'),
    %HfO2: pseudocubic a=5.1 fluorite
    %ZrO2: pseudocubic a=5.3 fluorite
    fA=(1-plotdata.Material(k).xHZO)*AtomicScatteringFactor('Hf',plotdata.QQ)+plotdata.Material(k).xHZO*AtomicScatteringFactor('Zr',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        plotdata.Material(k).V=plotdata.Material(k).xHZO*5.3^3+(1-plotdata.Material(k).xHZO)*5.1^3;%volume of the unit cell
        FLayer=StructureFactorFluorite001(fA,fO,plotdata.Q,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        plotdata.Material(k).V=3*(plotdata.Material(k).xHZO*5.3^3+(1-plotdata.Material(k).xHZO)*5.1^3);%volume of the unit cell
        FLayer=StructureFactorFluorite111(fA,fO,plotdata.Q,d);
    else
        msgbox('Error with (Hf_{1-x},Zr_x)O2 orientation');
    end;
end;

%% In2O3 %cubic bixbyite
if strcmp(plotdata.Material(k).Type,'In2O3'),
    plotdata.Material(k).V=10.29956^3; %mp-22598
    fIn=AtomicScatteringFactor('In',plotdata.QQ); 
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        z(1)=0; f(1,:)=2*fIn;
        z(2)=0.033444; f(2,:)=2*fIn;
        z(3)=0.110016; f(3,:)=4*fO;
        z(4)=0.117997; f(4,:)=4*fO;
        z(5)=0.154158; f(5,:)=4*fO;
        z(6)=0.250000; f(6,:)=8*fIn;
        z(7)=0.345842; f(7,:)=4*fO;
        z(8)=0.382003; f(8,:)=4*fO;
        z(9)= 0.389984; f(9,:)=4*fO;
        z(10)= 0.466556; f(10,:)=2*fIn;
        z(11)= 0.500000; f(11,:)=2*fIn;
        %Periodicity
        %z(1:11)=z(1:11)+0*0.5; %16 In + 24 O
        z(12:22)=z(1:11)+1*0.5; f(12:22,:)=f(1:11,:); %16 In + 24 O
        for j=1:22,
            FLayer=FLayer+f(j,:).*exp(1i*(plotdata.Q.*z(j).*d));
        end;
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        z(1)=0.0123092*6; f(1,:)=3*fO; 
        z(2)=0.0196382*6; f(2,:)=3*fO;
        z(3)=0.0243628*6; f(3,:)=3*fO;
        z(4)=0.0270232*6; f(4,:)=3*fO;
        z(5)=0.0360927*6; f(5,:)=6*fIn;
        z(6)=0.0416667*6; f(6,:)=4*fIn;
        z(7)=0.0472407*6; f(7,:)=6*fIn;
        z(8)=0.0563102*6; f(8,:)=3*fO;
        z(9)= 0.0589705*6; f(9,:)=3*fO;
        z(10)= 0.0636952*6; f(10,:)=3*fO;
        z(11)= 0.0710242*6; f(11,:)=3*fO;
        z(12)= 0.0956425*6; f(12,:)=3*fO;
        z(13)=0.1029715*6; f(13,:)=3*fO;
        z(14)=0.1076962*6; f(14,:)=3*fO;
        z(15)=0.1103565*6; f(15,:)=3*fO;
        z(16)=0.1194260*6; f(16,:)=6*fIn;
        z(17)=0.1250000*6; f(17,:)=4*fIn;
        z(18)=0.1305740*6; f(18,:)=6*fIn;
        z(19)=0.1396435*6; f(19,:)=3*fO;
        z(20)= 0.1423038*6; f(20,:)=3*fO;
        z(21)= 0.1470285*6; f(21,:)=3*fO;
        z(22)= 0.1543575*6; f(22,:)=3*fO;
       for j=1:22,
            FLayer=FLayer+f(j,:).*exp(1i*(plotdata.Q.*z(j).*d));
        end;
    else
        msgbox('Error with In2O3 orientation');
    end;
end;

%% LaAlO3
if strcmp(plotdata.Material(k).Type,'LaAlO3'),
    plotdata.Material(k).V=3.787^3; %Wikipedia... (room temperature pseudocubic)
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Al',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaAlO3 orientation');
    end;
end;

%% LaCoO3
if strcmp(plotdata.Material(k).Type,'LaCoO3'),
    plotdata.Material(k).V=3.816^3; %https://materialsproject.org/materials/mp-573180/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Co',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaCoO3 orientation');
    end;
end;

%% (La_x,Sr_{1-x})CoO3
if strcmp(plotdata.Material(k).Type,'(La_x,Sr_{1-x})CoO3'),
    %LaCoO3: %https://materialsproject.org/materials/mp-573180/
    %plotdata.Material(k).V=3.816^3; 
    %SrCoO3: https://materialsproject.org/materials/mp-505766/ alpha=beta=gamma=90
    %a=3.860; b=3.853; c=b; plotdata.Material(k).V=a*b*c;
    plotdata.Material(k).V=plotdata.Material(k).xLSCO*(3.816^3)+(1-plotdata.Material(k).xLSCO)*3.860*3.853^2;%volume of the unit cell   
    fA=plotdata.Material(k).xLSCO*AtomicScatteringFactor('La',plotdata.QQ)+(1-plotdata.Material(k).xLSCO)*AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Co',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with (La_x,Sr_{1-x})CoO3 orientation');
    end;
end;

%% LaCrO3
if strcmp(plotdata.Material(k).Type,'LaCrO3'),
    %pseudocubic a=3.885 (orthorhombic c-axis/2) a=5.513 A, b=5.476 A, c=7.759 A
    plotdata.Material(k).V=3.885^3;%volume of the unit cell    
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Cr',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaCrO3 orientation');
    end;
end;

%% La2CuO4 (CIF orthorhombic Cmca a=5.35 b=13.148 c=5.398): 
if strcmp(plotdata.Material(k).Type,'La2CuO4'),
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
if strcmp(plotdata.Material(k).Type,'LaFeO3'),
    plotdata.Material(k).V=3.959^3; %https://materialsproject.org/materials/mp-552676/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Fe',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaFeO3 orientation');
    end;
end;

%% LaMnO3
if strcmp(plotdata.Material(k).Type,'LaMnO3'),
    plotdata.Material(k).V=3.945^3; %https://materialsproject.org/materials/mp-19025/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Mn',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaMnO3 orientation');
    end;
end;

%% LaNiO3
if strcmp(plotdata.Material(k).Type,'LaNiO3'),
    plotdata.Material(k).V=3.857^3; %https://materialsproject.org/materials/mp-1075921/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaNiO3 orientation');
    end;
end;

%% La_xNi_yO3
if strcmp(plotdata.Material(k).Type,'La_xNi_yO3'),
    plotdata.Material(k).V=3.857^3; %https://materialsproject.org/materials/mp-1075921/
    fA=plotdata.Material(k).xLNO*AtomicScatteringFactor('La',plotdata.QQ); 
    fB=plotdata.Material(k).yLNO*AtomicScatteringFactor('Ni',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with La_xNi_yO3 orientation');
    end;
end;

%% La2NiMnO6 (double perovskites)
if strcmp(plotdata.Material(k).Type,'La2NiMnO6'),
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        plotdata.Material(k).V=3.871^3; %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5802844/
        fA=AtomicScatteringFactor('La',plotdata.QQ); 
        fB=0.5*AtomicScatteringFactor('Ni',plotdata.QQ)+0.5*AtomicScatteringFactor('Mn',plotdata.QQ);
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
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

%% La2NiO4:
% La8 Ni4 O16 (28 atoms)
% a = b = 5.60789 Å
% c = 12.42323 Å
% alpha = beta = gamma = 90.0000°
% Unit-cell volume = abc = 390.691341 Å^3
if strcmp(plotdata.Material(k).Type,'La2NiO4'),
    plotdata.Material(k).V=390.691341;
    zLa2NiO4(1)=0.07067; %2O
    zLa2NiO4(2)=0.11298; %2La
    zLa2NiO4(3)=0.22520; %O
    zLa2NiO4(4)=0.25000; %2O+2Ni
    zLa2NiO4(5)=0.27480; %O
    zLa2NiO4(6)=0.38702; %2La
    zLa2NiO4(7)=0.42933; %2O
    zLa2NiO4(8)=0.57067; %2O
    zLa2NiO4(9)=0.61298; %2La
    zLa2NiO4(10)=0.72520; %O
    zLa2NiO4(11)=0.75000;%2O+2Ni
    zLa2NiO4(12)=0.77480;%O
    zLa2NiO4(13)=0.88702;%2La
    zLa2NiO4(14)=0.92933;%2O

    fLa=AtomicScatteringFactor('La',plotdata.QQ); 
    fNi=AtomicScatteringFactor('Ni',plotdata.QQ);
    FLayer=2*fO.*exp(1i*(plotdata.Q.*zLa2NiO4(1).*d))...
        +2*fLa.*exp(1i*(plotdata.Q.*zLa2NiO4(2).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2NiO4(3).*d))...
        +(2*fO+2*fNi).*exp(1i*(plotdata.Q.*zLa2NiO4(4).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2NiO4(5).*d))...
        +2*fLa.*exp(1i*(plotdata.Q.*zLa2NiO4(6).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zLa2NiO4(7).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zLa2NiO4(8).*d))...
        +2*fLa.*exp(1i*(plotdata.Q.*zLa2NiO4(9).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2NiO4(10).*d))...
        +(2*fO+2*fNi).*exp(1i*(plotdata.Q.*zLa2NiO4(11).*d))...
        +fO.*exp(1i*(plotdata.Q.*zLa2NiO4(12).*d))...
        +2*fLa.*exp(1i*(plotdata.Q.*zLa2NiO4(13).*d))...
        +2*fO.*exp(1i*(plotdata.Q.*zLa2NiO4(14).*d));
end;

%% LaO
if strcmp(plotdata.Material(k).Type,'LaO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('La',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% LaTiO3
if strcmp(plotdata.Material(k).Type,'LaTiO3'),
    plotdata.Material(k).V=3.959^3; %https://materialsproject.org/materials/mp-8020/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaTiO3 orientation');
    end;
end;

%% LaVO3
if strcmp(plotdata.Material(k).Type,'LaVO3'),
    plotdata.Material(k).V=3.951^3; %https://materialsproject.org/materials/mp-19053/
    fA=AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('V',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LaVO3 orientation');
    end;
end;

%% LSMO
if strcmp(plotdata.Material(k).Type,'LSMO'),
    plotdata.Material(k).V=3.875^3;%http://ematweb.cmi.ua.ac.be/emat/pdf/1214.pdf
    fA=0.67*AtomicScatteringFactor('La',plotdata.QQ)+0.33*AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Mn',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with LSMO orientation');
    end;
end;

%% Mg3N2
% Ia-3 a=b=c=9.9528 alpha=bet=gamma=90 80 atoms
if strcmp(plotdata.Material(k).Type,'Mg3N2'),
    fMg=AtomicScatteringFactor('Mg',plotdata.QQ);
    fN=AtomicScatteringFactor('N',plotdata.QQ);
    a=9.9528; b=a; c=a;
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
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
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
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
if strcmp(plotdata.Material(k).Type,'MgO'),
    % cubic a=4.212 A
    % (Crystec Datasheets http://www.crystec.de/daten/mgo.pdf)
    plotdata.Material(k).V=4.212^3; 
    fMg=AtomicScatteringFactor('Mg',plotdata.QQ); 
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        zMgOc001(1)=0;                      % Mg-O
        zMgOc001(2)=0.5;                    % 2Mg-2O
        zMgOc001(3)=1;                      % MgO
        FLayer=(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(1).*d))...
        +2*(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(2).*d))...
        +(fMg+fO).*exp(1i*(plotdata.Q.*zMgOc001(3).*d));  %Structure factor F 
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
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
if strcmp(plotdata.Material(k).Type,'MnO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('Mn',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% MnO2
if strcmp(plotdata.Material(k).Type,'MnO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('Mn',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% MnTiO3
if strcmp(plotdata.Material(k).Type,'MnTiO3'),
    plotdata.Material(k).V=3.832^3; %https://materialsproject.org/materials/mp-19082/ V=112.583 for double unit cell
    fA=AtomicScatteringFactor('Mn',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with MnTiO3 orientation');
    end;
end;

%% MoS2 (P 63/m m c)
if strcmp(plotdata.Material(k).Type,'MoS2(P.63/m.m.c)'), 
    %hexagonal with primitive unit cell Mo2S4 a=b=3.16040, c=12.295, alpha=beta=90, gamma=120
    a=3.16040; b=a; c=12.295;
    plotdata.Material(k).V=sqrt(3)*a^2*c/2; %unit cell volume
    if strcmp(plotdata.Material(k).Orientation,'h(0001)'),
        % (0001): z along c-axis
        % inplane surface=sqrt(3)/2*a^2
        % V(unit cell)=sqrt(3)*a^2*c/2
        zMoS2h0001(1)=0.129;    %S
        zMoS2h0001(2)=0.250;    %Mo
        zMoS2h0001(3)=0.371;    %S
        zMoS2h0001(4)=0.629;    %S
        zMoS2h0001(5)=0.750;    %Mo
        zMoS2h0001(6)=0.871;    %S
        fMo=AtomicScatteringFactor('Mo',plotdata.QQ);
        fS=AtomicScatteringFactor('S',plotdata.QQ);
        FLayer=fS.*exp(1i*(plotdata.Q.*zMoS2h0001(1).*d))...
            +fMo.*exp(1i*(plotdata.Q.*zMoS2h0001(2).*d))...
            +fS.*exp(1i*(plotdata.Q.*zMoS2h0001(3).*d))...
            +fS.*exp(1i*(plotdata.Q.*zMoS2h0001(4).*d))...
            +fMo.*exp(1i*(plotdata.Q.*zMoS2h0001(5).*d))...
            +fS.*exp(1i*(plotdata.Q.*zMoS2h0001(6).*d));
       elseif or(strcmp(plotdata.Material(k).Orientation,'h(11-20)'),strcmp(plotdata.Material(k).Orientation,'h(10-10)')), 
        % (11-20): z along bisectrice a,b
        % inplane surface S=sqrt(3)*a*c
        % distance between planes = V(unit cell)/S = a/2
        % d=a/2;
        % TODO 
        % (10-10): z perpendicular to a
        % inplane surface S=ac
        % distance between planes = V(unit cell)/S = sqrt(3)a/2
        % d=a*sqrt(3)/2;
        % TODO;
    dlg = dialog('Position',[300 300 250 300],'Name','Warning');
    txt = uicontrol('Parent',dlg,...
               'Style','text',...
               'Position',[20 80 210 200],...
               'String',['The ' plotdata.Material(k).Orientation ' orientation for ' plotdata.Material(k).Type ' is not available yet. Please contact Celine.Lichtensteiger@unige.ch if needed.']);
    btn = uicontrol('Parent',dlg,...
               'Position',[85 20 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');    
    plotdata.Material(k).Orientation='h(0001)';
    plotdata.Materuak(k).plotdata.Material(k).choosehOrientation.Value=1;
    else
        msgbox('Error with MoS2 orientation');
    end;
end;

%% MoS2 (P -3 m 1)
if strcmp(plotdata.Material(k).Type,'MoS2(P.-3.m.1)'), 
    %hexagonal with primitive unit cell MoS2 a=b=3.19000, c=5.94500, alpha=beta=90, gamma=120
    a=3.19000; b=a; c=5.94500;
    plotdata.Material(k).V=sqrt(3)*a^2*c/2; %unit cell volume
    if strcmp(plotdata.Material(k).Orientation,'h(0001)'),
        % (0001): z along c-axis
        % inplane surface=sqrt(3)/2*a^2
        % V(unit cell)=sqrt(3)*a^2*c/2
        fMo=AtomicScatteringFactor('Mo',plotdata.QQ);
        fS=AtomicScatteringFactor('S',plotdata.QQ);
        z(1)=0; f(1,:)=fMo; 
        z(2)=0.25; f(2,:)=fS;
        z(3)=0.75; f(3,:)=fS;
        z(4)=1; f(4,:)=fMo;
        for j=1:4,
            FLayer=FLayer+f(j,:).*exp(1i*(plotdata.Q.*z(j).*d));
        end;
   elseif sum(strcmp(plotdata.Material(k).Orientation,{'h(11-20)','h(10-10)'})), 
        % (11-20): z along bisectrice a,b
        % inplane surface S=sqrt(3)*a*c
        % distance between planes = V(unit cell)/S = a/2
        % d=a/2;
        % TODO 
        % (10-10): z perpendicular to a
        % inplane surface S=ac
        % distance between planes = V(unit cell)/S = sqrt(3)a/2
        % d=a*sqrt(3)/2;
        % TODO;
    dlg = dialog('Position',[300 300 250 300],'Name','Warning');
    txt = uicontrol('Parent',dlg,...
               'Style','text',...
               'Position',[20 80 210 200],...
               'String',['The ' plotdata.Material(k).Orientation ' orientation for ' plotdata.Material(k).Type ' is not available yet. Please contact Celine.Lichtensteiger@unige.ch if needed.']);
    btn = uicontrol('Parent',dlg,...
               'Position',[85 20 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');    
    plotdata.Material(k).Orientation='h(0001)';
    plotdata.Materuak(k).plotdata.Material(k).choosehOrientation.Value=1;
    else
        msgbox('Error with MoS2 orientation');
    end;
end;

%% (Nd_{1-x},Ca_x)MnO3
if strcmp(plotdata.Material(k).Type,'(Nd_{1-x},Ca_x)MnO3'),
    %NdMnO3: pseudocubic a=3.834 (orthorhombic c-axis/2) https://materialsproject.org/materials/mp-20852/
    %CaMnO3: cubic a=3.802; https://materialsproject.org/materials/mp-1017467/
    plotdata.Material(k).V=plotdata.Material(k).xNCM*3.802^3+(1-plotdata.Material(k).xNCM)*3.834^3;%volume of the unit cell    
    fA=(1-plotdata.Material(k).xNCM)*AtomicScatteringFactor('Nd',plotdata.QQ)+plotdata.Material(k).xNCM*AtomicScatteringFactor('Ca',plotdata.QQ); 
    fB=AtomicScatteringFactor('Mn',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with (Nd_{1-x},Ca_x)MnO3 orientation');
    end;
end;

%% NdCaMn2O6 (double perovskites)
if strcmp(plotdata.Material(k).Type,'NdCaMn2O6'),
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        plotdata.Material(k).V=231.502/4; %https://materialsproject.org/materials/mp-1227144/
        fA=0.5*AtomicScatteringFactor('Nd',plotdata.QQ)+0.5*AtomicScatteringFactor('Ca',plotdata.QQ);
        fB=AtomicScatteringFactor('Mn',plotdata.QQ); 
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        plotdata.Material(k).V=231.502/2; %https://materialsproject.org/materials/mp-1227144/ 
        fNd=AtomicScatteringFactor('Nd',plotdata.QQ);
        fCa=AtomicScatteringFactor('Ca',plotdata.QQ); 
        fMn=AtomicScatteringFactor('Mn',plotdata.QQ);
        FLayer=fNd.*exp(1i*(plotdata.Q.*zdbperovskite111(1).*d))...
        +fMn.*exp(1i*(plotdata.Q.*zdbperovskite111(2).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(3).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(4).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(5).*d))...
        +fCa.*exp(1i*(plotdata.Q.*zdbperovskite111(6).*d))...
        +fMn.*exp(1i*(plotdata.Q.*zdbperovskite111(7).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(8).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(9).*d))...
        +fO.*exp(1i*(plotdata.Q.*zdbperovskite111(10).*d));
    else
        msgbox('Error with NdCaMn2O6 orientation');
    end;
end;

%% (Nd_x,La_{1-x})NiO3
if strcmp(plotdata.Material(k).Type,'(Nd_x,La_{1-x})NiO3'),
    %NdNiO3: pseudocubic a=3.861 cubicroot(230.143/4) https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
    %LaNiO3: cubic a=3.857; %https://materialsproject.org/materials/mp-1075921/
    plotdata.Material(k).V=plotdata.Material(k).xNLN*3.861^3+(1-plotdata.Material(k).xNLN)*3.857^3;%volume of the unit cell    
    fA=plotdata.Material(k).xNLN*AtomicScatteringFactor('Nd',plotdata.QQ)+(1-plotdata.Material(k).xNLN)*AtomicScatteringFactor('La',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with (Nd_x,La_{1-x})NiO3 orientation');
    end;
end;

%% NdNiO2 (cf CaCuO2) 
if strcmp(plotdata.Material(k).Type,'NdNiO2'), 
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
if strcmp(plotdata.Material(k).Type,'NdNiO3'),
    plotdata.Material(k).V=3.861^3; %https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
    fA=AtomicScatteringFactor('Nd',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with NdNiO3 orientation');
    end;
end;

%% Nd_xNi_yO3
if strcmp(plotdata.Material(k).Type,'Nd_xNi_yO3'),
    plotdata.Material(k).V=3.861^3; %https://materialsproject.org/materials/mp-22106/ V=230.143A 4 unit cells
    fA=plotdata.Material(k).xNNO*AtomicScatteringFactor('Nd',plotdata.QQ); 
    fB=plotdata.Material(k).yNNO*AtomicScatteringFactor('Ni',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with Nd_xNi_yO3 orientation');
    end;
end;

%% Nd2NiMnO6 (double perovskite)
if strcmp(plotdata.Material(k).Type,'Nd2NiMnO6'),%double perovskite
    %https://aip.scitation.org/doi/full/10.1063/1.4906989?ver=pdfcov
    %a=5.4162, b=5.4863, c=7.6718 V=a*b*c 4 unit cells
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        plotdata.Material(k).V=3.851^3;
        fA=AtomicScatteringFactor('Nd',plotdata.QQ); 
        fB=0.5*AtomicScatteringFactor('Ni',plotdata.QQ)+0.5*AtomicScatteringFactor('Mn',plotdata.QQ);
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
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
if strcmp(plotdata.Material(k).Type,'NdO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('Nd',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% (Nd_{1-x},Sr_x)MnO3
if strcmp(plotdata.Material(k).Type,'(Nd_{1-x},Sr_x)MnO3'),
    %NdMnO3: pseudocubic a=3.834 (orthorhombic c-axis/2) https://materialsproject.org/materials/mp-20852/
    %SrMnO3: pseudocubic a=3.869 (orthorhombic (V/4)^1/3) https://materialsproject.org/materials/mp-1017466/
    plotdata.Material(k).V=plotdata.Material(k).xNSM*3.869^3+(1-plotdata.Material(k).xNSM)*3.834^3;%volume of the unit cell    
    fA=(1-plotdata.Material(k).xNSM)*AtomicScatteringFactor('Nd',plotdata.QQ)+plotdata.Material(k).xNSM*AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Mn',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with (Nd_{1-x},Sr_x)MnO3 orientation');
    end;
end;

%% NiO2
if strcmp(plotdata.Material(k).Type,'NiO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% PbO
if strcmp(plotdata.Material(k).Type,'PbO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('Pb',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% PbNiO3
if strcmp(plotdata.Material(k).Type,'PbNiO3'),
    plotdata.Material(k).V=55.330; %https://materialsproject.org/materials/mp-974108/ V=55.330
    fA=AtomicScatteringFactor('Pb',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with PbNiO3 orientation');
    end;
end;

%% (Pb_x,Sr_{1-x})TiO3
if strcmp(plotdata.Material(k).Type,'(Pb_x,Sr_{1-x})TiO3'),
    %PbTiO3: known bulk values a=b=3.904, c=4.152
    %SrTiO3: known bulk values a=b=c=3.905
    plotdata.Material(k).V=plotdata.Material(k).xPST*(3.904^2*4.152)+(1-plotdata.Material(k).xPST)*3.905^3;%volume of the unit cell   
    fA=plotdata.Material(k).xPST*AtomicScatteringFactor('Pb',plotdata.QQ)+(1-plotdata.Material(k).xPST)*AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with (Pb_x,Sr_{1-x})TiO3 orientation');
    end;
end;

%% PbTiO3
if strcmp(plotdata.Material(k).Type,'PbTiO3'),
    plotdata.Material(k).V=3.904^2*4.152;
    fA=AtomicScatteringFactor('Pb',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with PbTiO3 orientation');
    end;
end;

%% Pb(Zr_x,Ti_{1-x})O3
if strcmp(plotdata.Material(k).Type,'Pb(Zr_x,Ti_{1-x})O3'),
    %PTO: 3.904^2*4.152 PZO:4.209^3 https://materialsproject.org/materials/mp-1068577/
    plotdata.Material(k).V=plotdata.Material(k).xPZT*(4.209^3)+(1-plotdata.Material(k).xPZT)*(3.904^2*4.152); %volume of the unit cell   
    fA=AtomicScatteringFactor('Pb',plotdata.QQ); 
    fB=plotdata.Material(k).xPZT*AtomicScatteringFactor('Zr',plotdata.QQ)+(1-plotdata.Material(k).xPZT)*AtomicScatteringFactor('Ti',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with Pb(Zr_x,Ti_{1-x})O3 orientation');
    end;
end;

%% PrBa2Cu3O7
if strcmp(plotdata.Material(k).Type,'PrBa2Cu3O7'),
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
if strcmp(plotdata.Material(k).Type,'PrNiO2'), 
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
if strcmp(plotdata.Material(k).Type,'PrNiO3'),
    plotdata.Material(k).V=3.872^3; %volume of the unit cell - https://materialsproject.org/materials/mp-22280/ V=232.262 4 unit cells
    fA=AtomicScatteringFactor('Pr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with PrNiO3 orientation');
    end;
end;

%% PrVO3
if strcmp(plotdata.Material(k).Type,'PrVO3'),
    plotdata.Material(k).V=3.936^3; %volume of the unit cell - https://materialsproject.org/materials/mp-1069346/
    fA=AtomicScatteringFactor('Pr',plotdata.QQ); 
    fB=AtomicScatteringFactor('V',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with PrVO3 orientation');
    end;
end;

%% RuO2
if strcmp(plotdata.Material(k).Type,'RuO2'), %%BO2-type monolayer
       plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
       fB=AtomicScatteringFactor('Ru',plotdata.QQ);
       FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% SmNiO3
if strcmp(plotdata.Material(k).Type,'SmNiO3'),
    plotdata.Material(k).V=3.794^3; %https://materialsproject.org/materials/mp-1099668/
    fA=AtomicScatteringFactor('Sm',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SmNiO3 orientation');
    end;
end;

%% Sr3Al2O6 https://materialsproject.org/materials/mp-3393/  
%Sr72Al48O144 in cubic unit cell a=b=c=15.999914A
if strcmp(plotdata.Material(k).Type,'Sr3Al2O6'),
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
if strcmp(plotdata.Material(k).Type,'SrCoO2.5'),
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
if strcmp(plotdata.Material(k).Type,'SrCoO3'),
    %https://materialsproject.org/materials/mp-505766/ alpha=beta=gamma=90
    a=3.860; b=3.853; c=b; 
    plotdata.Material(k).V=a*b*c;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Co',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrCoO3 orientation');
    end;
end;

%% SrCrO3
if strcmp(plotdata.Material(k).Type,'SrCrO3'),
    plotdata.Material(k).V=3.8185^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Cr',plotdata.QQ);    
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrCrO3 orientation');
    end;
end;

%% SrCuO2 (cf CaCuO2) 
if strcmp(plotdata.Material(k).Type,'SrCuO2'), 
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
if strcmp(plotdata.Material(k).Type,'SrIrO3'),
    %https://materialsproject.org/materials/mp-1016848/
    %cubic
    a=3.998;
    plotdata.Material(k).V=a^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ir',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrIrO3 orientation');
    end;
end;

%% SrMoO3
if strcmp(plotdata.Material(k).Type,'SrMoO3'),
    %https://materialsproject.org/materials/mp-18747/
    %cubic
    a=4.082;
    plotdata.Material(k).V=a^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Mo',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrMoO3 orientation');
    end;
end;

%% SrO
if strcmp(plotdata.Material(k).Type,'SrO'),
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fA=AtomicScatteringFactor('Sr',plotdata.QQ);
    FLayer=fA.*exp(1i*(plotdata.Q.*zAO(1).*d))+fO.*exp(1i*(plotdata.Q.*zAO(2).*d));
end;

%% SrO2
if strcmp(plotdata.Material(k).Type,'SrO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('Sr',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% SrRuO3
if strcmp(plotdata.Material(k).Type,'SrRuO3'),
    %https://materialsproject.org/materials/mp-4346/
    %cubic
    a=3.985;
    plotdata.Material(k).V=a^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ru',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrRuO3 orientation');
    end;
end;

%% SrTiO3
if strcmp(plotdata.Material(k).Type,'SrTiO3'),
    %cubic
    a=3.905;
    plotdata.Material(k).V=a^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrTiO3 orientation');
    end;
end;

%% SrVO3
if strcmp(plotdata.Material(k).Type,'SrVO3'),
    %https://materialsproject.org/materials/mp-18717/
    %cubic
    a=3.901;
    plotdata.Material(k).V=a^3;
    fA=AtomicScatteringFactor('Sr',plotdata.QQ); 
    fB=AtomicScatteringFactor('V',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with SrVO3 orientation');
    end;
end;

%% TiO2
if strcmp(plotdata.Material(k).Type,'TiO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('Ti',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% Tm3Fe5O12
if strcmp(plotdata.Material(k).Type,'Tm3Fe5O12'),
    % Thulium Iron garnet (TmIG)
    a=12.325; b=a; c=a; %Landolt-Bornstein
    plotdata.Material(k).V=a*b*c/4; %volume of the "new" unit cell (see below)
    fTm=AtomicScatteringFactor('Tm',plotdata.QQ); 
    fFe=AtomicScatteringFactor('Fe',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
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
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
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
if strcmp(plotdata.Material(k).Type,'VO2'), %%BO2-type monolayer
    plotdata.Material(k).V=4*4*2; %unit cell volume (ARBITRARY AS THIS IS "ARTIFICIAL")
    fB=AtomicScatteringFactor('V',plotdata.QQ);
    FLayer=fB.*exp(1i*(plotdata.Q.*zBO2(1).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(2).*d))+fO.*exp(1i*(plotdata.Q.*zBO2(3).*d));        
end;

%% WS2
if strcmp(plotdata.Material(k).Type,'WS2'), 
    %hexagonal with primitive unit cell W2S4 a=b=3.15320, c=12.32300 , alpha=beta=90, gamma=120
    a=3.15320; b=a; c=12.32300;
    plotdata.Material(k).V=sqrt(3)*a^2*c/2; %unit cell volume
    if strcmp(plotdata.Material(k).Orientation,'h(0001)'),
        % (0001): z along c-axis
        % inplane surface=sqrt(3)/2*a^2
        % V(unit cell)=sqrt(3)*a^2*c/2
        zWS2h0001(1)=0.122; %S
        zWS2h0001(2)=0.250; %W
        zWS2h0001(3)=0.377; %S
        zWS2h0001(4)=0.623; %S
        zWS2h0001(5)=0.750; %W
        zWS2h0001(6)=0.878; %S
        fW=AtomicScatteringFactor('W',plotdata.QQ);
        fS=AtomicScatteringFactor('S',plotdata.QQ);
        FLayer=fS.*exp(1i*(plotdata.Q.*zWS2h0001(1).*d))...
            +fW.*exp(1i*(plotdata.Q.*zWS2h0001(2).*d))...
            +fS.*exp(1i*(plotdata.Q.*zWS2h0001(3).*d))...
            +fS.*exp(1i*(plotdata.Q.*zWS2h0001(4).*d))...
            +fW.*exp(1i*(plotdata.Q.*zWS2h0001(5).*d))...
            +fS.*exp(1i*(plotdata.Q.*zWS2h0001(6).*d));
    else
        msgbox('Error with WS2 orientation');
    end;
end;

%% YBCO:
if strcmp(plotdata.Material(k).Type,'YBa2Cu3O7'),
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
if strcmp(plotdata.Material(k).Type,'YBiO3'),
    plotdata.Material(k).V=4.408542^3; %unit cell volume https://materialsproject.org/materials/mp-13598/
    fA=AtomicScatteringFactor('Y',plotdata.QQ); 
    fB=AtomicScatteringFactor('Bi',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with YBiO3 orientation');
    end;
end;

%% Y3Fe5O12
if strcmp(plotdata.Material(k).Type,'Y3Fe5O12'),
    % Ytrium Iron garnet (YIG)
    a=12.376; b=a; c=a; %Landolt-Bornstein
    plotdata.Material(k).V=a*b*c/4; %volume of the "new" unit cell (see below)
    fY=AtomicScatteringFactor('Y',plotdata.QQ); 
    fFe=AtomicScatteringFactor('Fe',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
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
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
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
if strcmp(plotdata.Material(k).Type,'(Y_xTm_{3-x})Fe5O12'),
    %YFe5O12: pseudo-cubic a=12.376 Landolt/Bornstein
    %TmFe5O12: pseudo- cubic a=12.325 Landolt/Bornstein
    aTmIG=12.325; aYIG=12.376; %Landolt-Bornstein
    VTmIG=aTmIG^3/4; VYIG=aYIG^3/4;
    plotdata.Material(k).V=(plotdata.Material(k).xYTmIG*VYIG+(3-plotdata.Material(k).xYTmIG)*VTmIG)/3; %volume of the "new" unit cell
    fYTm=(1/3)*(plotdata.Material(k).xYTmIG*AtomicScatteringFactor('Y',plotdata.QQ)+(3-plotdata.Material(k).xYTmIG)*AtomicScatteringFactor('Tm',plotdata.QQ)); 
    fFe=AtomicScatteringFactor('Fe',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
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
        zYTmIG001=(plotdata.Material(k).xYTmIG*zYIG001+(3-plotdata.Material(k).xYTmIG)*zTmIG001)/3;
        FLayer=(4*fFe+2*fYTm).*exp(1i*(plotdata.Q.*zYTmIG001(1).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(2).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(3).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(4).*d))...
        +(2*fFe+2*fYTm).*exp(1i*(plotdata.Q.*zYTmIG001(5).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(6).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(7).*d))...
        +4*fO.*exp(1i*(plotdata.Q.*zYTmIG001(8).*d))...
        +(4*fFe+2*fYTm).*exp(1i*(plotdata.Q.*zYTmIG001(9).*d));
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
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
        zYTmIG111=(plotdata.Material(k).xYTmIG*zYIG111+(3-plotdata.Material(k).xYTmIG)*zTmIG111)/3;
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
if strcmp(plotdata.Material(k).Type,'YNiO3'),
    %https://materialsproject.org/materials/mvc-15448/
    %a=3.753 b=c=3.759 alpha=89.967 beta=gamma=90 V=53.025
    plotdata.Material(k).V=53.025; %unit cell volume
    fA=AtomicScatteringFactor('Y',plotdata.QQ); 
    fB=AtomicScatteringFactor('Ni',plotdata.QQ);
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite001,d);
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        FLayer=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite111,d);
    else
        msgbox('Error with YNiO3 orientation');
    end;
end;

%% ZnMgO
if strcmp(plotdata.Material(k).Type,'(Zn_x,Mg_{1-x})O'),
    % same structure as ZnO but with Mg substituting Zn for Mg<30%
    % ZnO: hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
    % substituting with Mg increases out-of-plane lattice parameter
    a=3.2494; b=a; c=5.20380;
    plotdata.Material(k).V=sqrt(3)*a^2*c/2; %unit cell volume
    fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
    fMg=AtomicScatteringFactor('Mg',plotdata.QQ);
    x=plotdata.Material(k).xZMO;
    if strcmp(plotdata.Material(k).Orientation,'h(0001)'),
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
    elseif strcmp(plotdata.Material(k).Orientation,'h(11-20)'), 
        %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
        % (11-20): z along bisectrice a,b
        % inplane surface S=sqrt(3)*a*c
        % distance between planes = V(unit cell)/S = a/2
        zZnMgOh11bar20(1)=0;     %2 Zn + 2 O
        %a=3.2494;
        %d=a/2;
        FLayer=(2*(x*fZn+(1-x)*fMg)+2*fO).*exp(1i*(plotdata.Q.*zZnMgOh11bar20(1).*d));
    elseif strcmp(plotdata.Material(k).Orientation,'h(10-10)'),
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
% Ia-3 a=b=c=9.76910 alpha=bet=gamma=90 80 atomes
if strcmp(plotdata.Material(k).Type,'Zn3N2'),
    fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
    fN=AtomicScatteringFactor('N',plotdata.QQ);
    a=9.76910; b=a; c=a;
    if strcmp(plotdata.Material(k).Orientation,'pc(001)'),
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
    elseif strcmp(plotdata.Material(k).Orientation,'pc(111)'),
        %need to use hexagonal unit cell (6*V) -> result in 80*6=480 atoms divided in 12 sub-layers (12Zn/16N/12Zn)
        % -> new unit cell is 2*a*sqrt(3)/12=a/(2*sqrt(3))
        plotdata.Material(k).V=a^3/2;
        zZn3N2(1)=0.15359981984;    %3Zn
        zZn3N2(2)=0.24720025386;    %3Zn
        zZn3N2(3)=0.25639987307;    %3Zn
        zZn3N2(4)=0.34280005323;    %3Zn
        zZn3N2(5)=0.45679990992;    %6N
        zZn3N2(6)=0.50000000000;    %4N
        zZn3N2(7)=0.54320009008;    %6N
        zZn3N2(8)=0.65719994677;    %3Zn
        zZn3N2(9)=0.74359992220;    %3Zn
        zZn3N2(10)=0.75279974614;   %3Zn
        zZn3N2(11)=0.84640018016;   %3Zn
        FLayer=3*fZn.*exp(1i*(plotdata.Q.*zZn3N2(1).*d))...
        +3*fZn.*exp(1i*(plotdata.Q.*zZn3N2(2).*d))...
        +3*fZn.*exp(1i*(plotdata.Q.*zZn3N2(3).*d))...
        +3*fZn.*exp(1i*(plotdata.Q.*zZn3N2(4).*d))...
        +6*fN.*exp(1i*(plotdata.Q.*zZn3N2(5).*d))...
        +4*fN.*exp(1i*(plotdata.Q.*zZn3N2(6).*d))...
        +6*fN.*exp(1i*(plotdata.Q.*zZn3N2(7).*d))...
        +3*fZn.*exp(1i*(plotdata.Q.*zZn3N2(8).*d))...
        +3*fZn.*exp(1i*(plotdata.Q.*zZn3N2(9).*d))...
        +3*fZn.*exp(1i*(plotdata.Q.*zZn3N2(10).*d))...
        +3*fZn.*exp(1i*(plotdata.Q.*zZn3N2(11).*d)); 
    else
        msgbox('Error with Zn3N2 orientation');
    end;
end;

%% ZnO
if strcmp(plotdata.Material(k).Type,'ZnO'), 
    %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
    a=3.2494; b=a; c=5.20380;
    plotdata.Material(k).V=sqrt(3)*a^2*c/2; %unit cell volume
    if strcmp(plotdata.Material(k).Orientation,'h(0001)'),
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
    elseif strcmp(plotdata.Material(k).Orientation,'h(11-20)'), 
        %hexagonal with primitive unit cell Zn2O2 a=b=3.2494, c=5.20380, alpha=beta=90?, gamma=120?
        % (11-20): z along bisectrice (a,b)
        % inplane surface S=sqrt(3)*a*c
        % distance between planes = V(unit cell)/S = a/2
        zZnOh11bar20(1)=0;     %2 Zn + 2 O
        %a=3.2494;
        %d=a/2;
        fZn=AtomicScatteringFactor('Zn',plotdata.QQ);
        FLayer=(2*fZn+2*fO).*exp(1i*(plotdata.Q.*zZnOh11bar20(1).*d));
    elseif strcmp(plotdata.Material(k).Orientation,'h(10-10)'),
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

FLayer=FLayer*d/plotdata.Material(k).V; %renormalized by the in-plane surface of the unit cell