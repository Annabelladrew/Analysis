%Created by Celine Lichtensteiger
%Calculates the XRD intensity and c-axis profile of the heterostructure
%*********************************
global plotdata;
plotdata.ThicknessBelow=plotdata.Substrate.N001*plotdata.Substrate.c;
plotdata.c=[];
for k=1:5,
    plotdata.Material(k).expzm=zeros(size(plotdata.Q));
end;
expzm(1);
for nRepetition=1:plotdata.Repetition.N,
    for k=2:4,
        expzm(k);
    end;
end;
expzm(5);
%Scattering amplitude g:
g=plotdata.Substrate.g;
for k=1:5,
    g=g+plotdata.Material(k).expzm;
end,
%scattered intensity I:
I=g.*conj(g);
plotdata.fit.y=I/max(I);
PlotFitAndData;