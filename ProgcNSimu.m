%Created by Celine Lichtensteiger
%Calculates the XRD intensity and c-axis profile of the heterostructure
%*********************************
global fitdata plotdata;
plotdata.ThicknessBelow=plotdata.Substrate.thickness;
plotdata.d=[];
for k=1:6,
    plotdata.Material(k).expzm=zeros(size(plotdata.Q));
end;
expzm(1,0);
for nRepetition=1:fitdata.Repetition.N,
    for k=2:5,
        expzm(k,nRepetition);
    end;
end;
expzm(6,0);
%Scattering amplitude g:
g=plotdata.Substrate.g;
for k=1:6,
    g=g+plotdata.Material(k).expzm;
end,
%scattered intensity I:
I=g.*conj(g);
plotdata.fit.y=I/max(I);
PlotFitAndData;