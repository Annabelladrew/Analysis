% Created by Celine Lichtensteiger
% Calculates the contribution from each layer material to:
% plotdata.Material(k).expzm
% plotdata.ThicknessBelow
% plotdata.d
%*********************************

function expzm(k,nRepetition)
global fitdata plotdata;

    if strcmp(fitdata.Material(k).dDistribution,'constant')
        plotdata.Material(k).Thickness=fitdata.Material(k).N*fitdata.Material(k).d;
        plotdata.Material(k).F=FLayer(k,fitdata.Material(k).d);
        for zm=1:fitdata.Material(k).N,
            plotdata.ThicknessBelow=plotdata.ThicknessBelow+fitdata.Material(k).d;
            plotdata.d=[plotdata.d fitdata.Material(k).d];
            ThicknessAbove=plotdata.TotalThickness-plotdata.ThicknessBelow;
            plotdata.Material(k).expzm=plotdata.Material(k).expzm+plotdata.Material(k).F.*exp(1i*plotdata.Q.*plotdata.ThicknessBelow).*exp(-(ThicknessAbove)./plotdata.mu);
        end;
        
    elseif strcmp(fitdata.Material(k).dDistribution,'exp')
        plotdata.Material(k).Thickness=0;
        fitdata.Material(k).Cexp=fitdata.Material(k).Meandexp-fitdata.Material(k).Aexp/fitdata.Material(k).N*(1-exp(fitdata.Material(k).N/fitdata.Material(k).Bexp))/(1-exp(1/fitdata.Material(k).Bexp));
        for zm=1:fitdata.Material(k).N,
            dzm=fitdata.Material(k).Aexp*exp((zm-1)/fitdata.Material(k).Bexp)+fitdata.Material(k).Cexp;
            plotdata.ThicknessBelow=plotdata.ThicknessBelow+dzm;
            ThicknessAbove=plotdata.TotalThickness-plotdata.ThicknessBelow;
            plotdata.Material(k).F=FLayer(k,dzm);
            plotdata.Material(k).expzm=plotdata.Material(k).expzm+plotdata.Material(k).F.*exp(1i*plotdata.Q.*plotdata.ThicknessBelow).*exp(-(ThicknessAbove)./plotdata.mu);
            plotdata.Material(k).Thickness=plotdata.Material(k).Thickness+dzm;
            plotdata.d=[plotdata.d dzm];
        end;
        
    elseif strcmp(fitdata.Material(k).dDistribution,'Jacobi')
        plotdata.Material(k).Thickness=0;
        m=fitdata.Material(k).m;
        if m==1, %Jacobi = step function (no intermixing)
            for zm=1:fitdata.Material(k).N,
                dzm=fitdata.Material(k).dCenter;
                plotdata.ThicknessBelow=plotdata.ThicknessBelow+dzm;
                ThicknessAbove=plotdata.TotalThickness-plotdata.ThicknessBelow;
                plotdata.Material(k).F=FLayer(k,dzm);
                plotdata.Material(k).expzm=plotdata.Material(k).expzm+plotdata.Material(k).F.*exp(1i*plotdata.Q.*plotdata.ThicknessBelow).*exp(-(ThicknessAbove)./plotdata.mu);
                plotdata.Material(k).Thickness=plotdata.Material(k).Thickness+dzm;
                plotdata.d=[plotdata.d dzm];
            end;
        else
            %u=[0:1:fitdata.Material(k).N];
            u=[0:1:fitdata.Material(k).N*2];
            %[SN,CN,DN]=ellipj(u*ellipke(m)/(fitdata.Material(k).N-1)*2,m);
            [SN,CN,DN]=ellipj(u*ellipke(m)/fitdata.Material(k).N,m);
            disp(fitdata.Material(k).Type);
            for zm=1:fitdata.Material(k).N,
                SNzm=SN(2*zm);
                dzm=SNzm*(fitdata.Material(k).dCenter-fitdata.Material(k).dBorder)+fitdata.Material(k).dBorder;
                plotdata.ThicknessBelow=plotdata.ThicknessBelow+dzm;
                ThicknessAbove=plotdata.TotalThickness-plotdata.ThicknessBelow;
                
                
                
                
                
                
                %Add atomic intermixing: proportion given by Jacobi function
                %atoms depend on previous or next layer type
                if zm<(fitdata.Material(k).N+1)/2,
                    %zm=zm,
                    %disp('case mix with previous');
                    if k==1,
                        plotdata.Material(k).F=intermixing(SNzm,dzm,1,0);
                    
                    elseif k==2, 
                        if nRepetition==1,
                            if strcmp(fitdata.Material(1).Type,'none'),
                                plotdata.Material(k).F=intermixing(SNzm,dzm,2,0);
                            else plotdata.Material(k).F=intermixing(SNzm,dzm,2,1);
                            end;
                        else
                            if strcmp(fitdata.Material(5).Type,'none'),
                                if strcmp(fitdata.Material(4).Type,'none'),
                                    if strcmp(fitdata.Material(3).Type,'none'),
                                        plotdata.Material(k).F=FLayer(k,dzm);
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,2,3);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,2,4);
                                end;
                            else plotdata.Material(k).F=intermixing(SNzm,dzm,2,5);
                            end;
                        end;
                        
                    elseif k==3, 
                        if nRepetition==1,
                            if strcmp(fitdata.Material(2).Type,'none'),
                                if strcmp(fitdata.Material(1).Type,'none'),
                                    plotdata.Material(k).F=intermixing(SNzm,dzm,3,0);
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,3,1);
                                end;
                            else plotdata.Material(k).F=intermixing(SNzm,dzm,3,2);
                            end;
                        else
                            if strcmp(fitdata.Material(2).Type,'none'),
                                if strcmp(fitdata.Material(5).Type,'none'),
                                    if strcmp(fitdata.Material(4).Type,'none'),
                                        plotdata.Material(k).F=FLayer(k,dzm);
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,3,4);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,3,5);
                                end;
                            else plotdata.Material(k).F=intermixing(SNzm,dzm,3,2);
                            end
                        end;
                        
                        
                    elseif k==4,
                        if nRepetition==1,
                            if strcmp(fitdata.Material(3).Type,'none'),
                                if strcmp(fitdata.Material(2).Type,'none'),
                                    if strcmp(fitdata.Material(1).Type,'none'),
                                        plotdata.Material(k).F=intermixing(SNzm,dzm,4,0);
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,4,1);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,4,2);
                                end;
                            else plotdata.Material(k).F=intermixing(SNzm,dzm,4,3);
                            end;
                        else
                            if strcmp(fitdata.Material(3).Type,'none'),
                                if strcmp(fitdata.Material(2).Type,'none'),
                                    if strcmp(fitdata.Material(5).Type,'none');
                                        plotdata.Material(k).F=FLayer(k,dzm);
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,4,5);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,4,2);
                                end;
                            else plotdata.Material(k).F=intermixing(SNzm,dzm,4,5);
                            end;
                        end;
                        
                    elseif k==5,
                        if nRepetition==1,
                            if strcmp(fitdata.Material(4).Type,'none'),
                                if strcmp(fitdata.Material(3).Type,'none'),
                                    if strcmp(fitdata.Material(2).Type,'none'),
                                        if strcmp(fitdata.Material(1).Type,'none'),
                                            plotdata.Material(k).F=intermixing(SNzm,dzm,5,0);
                                        else plotdata.Material(k).F=intermixing(SNzm,dzm,5,1);
                                        end;
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,5,2);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,5,3);
                                end;
                            else plotdata.Material(k).F=intermixing(SNzm,dzm,5,4);
                            end;
                        else
                            if strcmp(fitdata.Material(4).Type,'none'),
                                if strcmp(fitdata.Material(3).Type,'none'),
                                    if strcmp(fitdata.Material(2).Type,'none'),
                                        plotdata.Material(k).F=FLayer(k,dzm);
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,5,2);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,5,3);
                                end;
                            else plotdata.Material(k).F=intermixing(SNzm,dzm,5,4);
                            end;
                        end;
                        
                    elseif k==6,
                        if strcmp(fitdata.Material(5).Type,'none'),
                            if strcmp(fitdata.Material(4).Type,'none'),
                                if strcmp(fitdata.Material(3).Type,'none'),
                                    if strcmp(fitdata.Material(2).Type,'none'),
                                        if strcmp(fitdata.Material(1).Type,'none'),
                                            plotdata.Material(k).F=intermixing(SNzm,dzm,6,0);
                                        else plotdata.Material(k).F=intermixing(SNzm,dzm,6,1);
                                        end;
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,6,2);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,6,3);
                                end;
                            else plotdata.Material(k).F=intermixing(SNzm,dzm,6,4);
                            end;
                        else plotdata.Material(k).F=intermixing(SNzm,dzm,6,5);
                        end;
                    end;
                    
                    elseif zm==(fitdata.Material(k).N+1)/2,
                        %zm=zm,
                        %disp('case no mix (center)');
                        plotdata.Material(k).F=FLayer(k,dzm);

                    elseif zm>(fitdata.Material(k).N+1)/2,
                        %zm=zm,
                        %disp('case mix with next');
                        if k==1,
                            if strcmp(fitdata.Material(2).Type,'none'),
                                if strcmp(fitdata.Material(3).Type,'none'),
                                    if strcmp(fitdata.Material(4).Type,'none'),
                                        if strcmp(fitdata.Material(5).Type,'none'),
                                            if strcmp(fitdata.Material(6).Type,'none'),
                                                plotdata.Material(k).F=FLayer(k,dzm);
                                            else plotdata.Material(k).F=intermixing(SNzm,dzm,1,6);
                                            end;
                                        else plotdata.Material(k).F=intermixing(SNzm,dzm,1,5);
                                        end;
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,1,4);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,1,3);
                                end;
                            else plotdata.Material(k).F=intermixing(SNzm,dzm,1,2);
                            end;
                        elseif k==2,
                            if nRepetition==fitdata.Repetition.N,
                                if strcmp(fitdata.Material(3).Type,'none'),
                                    if strcmp(fitdata.Material(4).Type,'none'),
                                        if strcmp(fitdata.Material(5).Type,'none'),
                                            if strcmp(fitdata.Material(6).Type,'none'),
                                                plotdata.Material(k).F=FLayer(k,dzm);
                                            else plotdata.Material(k).F=intermixing(SNzm,dzm,2,6);
                                            end;
                                        else plotdata.Material(k).F=intermixing(SNzm,dzm,2,5);
                                        end;
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,2,4);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,2,3);
                                end;
                            else
                                if strcmp(fitdata.Material(3).Type,'none'),
                                    if strcmp(fitdata.Material(4).Type,'none'),
                                        if strcmp(fitdata.Material(5).Type,'none'),
                                            plotdata.Material(k).F=FLayer(k,dzm);
                                        else plotdata.Material(k).F=intermixing(SNzm,dzm,2,5);
                                        end;
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,2,4);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,2,3);
                                end;
                            end;
                        elseif k==3,
                            if nRepetition==fitdata.Repetition.N,
                                if strcmp(fitdata.Material(4).Type,'none'),
                                    if strcmp(fitdata.Material(5).Type,'none'),
                                        if strcmp(fitdata.Material(6).Type,'none'),
                                            plotdata.Material(k).F=FLayer(k,dzm);
                                        else plotdata.Material(k).F=intermixing(SNzm,dzm,3,6);
                                        end;
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,3,5);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,3,4);
                                end;
                            else
                                if strcmp(fitdata.Material(4).Type,'none'),
                                    if strcmp(fitdata.Material(5).Type,'none'),
                                        if strcmp(fitdata.Material(2).Type,'none'),
                                            plotdata.Material(k).F=FLayer(k,dzm);
                                            else plotdata.Material(k).F=intermixing(SNzm,dzm,3,2);
                                        end;
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,3,5);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,3,4);
                                end;
                            end;
                        elseif k==4,
                            if nRepetition==fitdata.Repetition.N,
                                if strcmp(fitdata.Material(5).Type,'none'),
                                    if strcmp(fitdata.Material(6).Type,'none'),
                                        plotdata.Material(k).F=FLayer(k,dzm);
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,4,6);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,4,5);
                                end;
                            else
                                if strcmp(fitdata.Material(5).Type,'none'),
                                    if strcmp(fitdata.Material(2).Type,'none'),
                                        if strcmp(fitdata.Material(3).Type,'none'),
                                            plotdata.Material(k).F=FLayer(k,dzm);
                                        else plotdata.Material(k).F=intermixing(SNzm,dzm,4,3);
                                        end;
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,4,2);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,4,5);
                                end;
                            end;
                            
                         elseif k==5,
                            if nRepetition==fitdata.Repetition.N,
                                if strcmp(fitdata.Material(6).Type,'none'),
                                    plotdata.Material(k).F=FLayer(k,dzm);
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,5,6);
                                end;
                            else
                                if strcmp(fitdata.Material(2).Type,'none'),
                                    if strcmp(fitdata.Material(3).Type,'none'),
                                        if strcmp(fitdata.Material(4).Type,'none'),
                                            plotdata.Material(k).F=FLayer(k,dzm);
                                        else plotdata.Material(k).F=intermixing(SNzm,dzm,5,4);
                                        end;
                                    else plotdata.Material(k).F=intermixing(SNzm,dzm,5,3);
                                    end;
                                else plotdata.Material(k).F=intermixing(SNzm,dzm,5,2);
                                end;
                            end;
                        elseif k==6,
                            plotdata.Material(k).F=FLayer(k,dzm);
                        end;
                    end;
                    
                    plotdata.Material(k).expzm=plotdata.Material(k).expzm+plotdata.Material(k).F.*exp(1i*plotdata.Q.*plotdata.ThicknessBelow).*exp(-(ThicknessAbove)./plotdata.mu);
                    plotdata.Material(k).Thickness=plotdata.Material(k).Thickness+dzm;
                    plotdata.d=[plotdata.d dzm];
                    %disp(dzm);
            end;
        end;
    elseif strcmp(fitdata.Material(k).dDistribution,'upload')
        plotdata.Material(k).Thickness=0;
        for zm=1:fitdata.Material(k).N,
            dzm=plotdata.Material(k).dupload(zm);
            plotdata.ThicknessBelow=plotdata.ThicknessBelow+dzm;
            ThicknessAbove=plotdata.TotalThickness-plotdata.ThicknessBelow;
            plotdata.Material(k).F=FLayer(k,dzm);
            plotdata.Material(k).expzm=plotdata.Material(k).expzm+plotdata.Material(k).F.*exp(1i*plotdata.Q.*plotdata.ThicknessBelow).*exp(-(ThicknessAbove)./plotdata.mu);
            plotdata.Material(k).Thickness=plotdata.Material(k).Thickness+dzm;
            plotdata.d=[plotdata.d dzm];
        end;
    else
        warning('Error with toggle function for dDistribution in expzm')
    end;

    function y=intermixing(SNzm,d,ki,kj)
        global fitdata plotdata;
        %SNzm=SNzm,
        %d=d,
        %ki=ki,
        %kj=kj,
    if ki==0,%substrate
        if strcmp(fitdata.Substrate.Type,'SrTiO3'), %intermixing with STO substrate
            fA=AtomicScatteringFactor('Sr',plotdata.QQ); fB=AtomicScatteringFactor('Ti',plotdata.QQ); fO=AtomicScatteringFactor('O',plotdata.QQ);
            zperovskite(1)=0; zperovskite(2)=0.5; zperovskite(3)=0; zperovskite(4)=0.5; zperovskite(5)=0.5; FLayerSTO=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite,d);
            y=(0.5+SNzm/2)*FLayerSTO+(0.5-SNzm/2)*FLayer(kj,d);
            %disp(['Intermixing: ',0.5+SNzm/2,'* STO substrate +',(0.5-SNzm/2),'*',plotdata.Material(kj).Type]);
        else y=FLayer(kj,c); %no intermixing with substrate if different than STO
        end
    elseif kj==0,%substrate
        if strcmp(fitdata.Substrate.Type,'SrTiO3'), %intermixing with STO substrate
            fA=AtomicScatteringFactor('Sr',plotdata.QQ); fB=AtomicScatteringFactor('Ti',plotdata.QQ); fO=AtomicScatteringFactor('O',plotdata.QQ);
            zperovskite(1)=0; zperovskite(2)=0.5; zperovskite(3)=0; zperovskite(4)=0.5; zperovskite(5)=0.5; FLayerSTO=StructureFactor(fA,fB,fO,plotdata.Q,zperovskite,d);
            y=(0.5+SNzm/2)*FLayer(ki,c)+(0.5-SNzm/2)*FLayerSTO;
            %disp(['Intermixing: ',0.5+SNzm/2,' * ',fitdata.Material(ki).Type,' + ',(0.5-SNzm/2),'* STO substrate']);
        else y=FLayer(ki,d);
        end
    else
        y=(0.5+SNzm/2)*FLayer(ki,d)+(0.5-SNzm/2)*FLayer(kj,d);
        disp(['Intermixing: ',0.5+SNzm/2,' * ',fitdata.Material(ki).Type,' + ',(0.5-SNzm/2),' * ',fitdata.Material(kj).Type]);
    end;