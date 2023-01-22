% Created by Celine Lichtensteiger
% Calculates the structure factor for materials with a fluorite structure
% AB_2 oriented 001
%% Fluorite 4xAB2 => 12 atoms => A4B8
%*********************************
function[F]=StructureFactorFluorite001(fA,fB,Q,d)

z(1)=0;                     % 2 A
z(2)=0.25;                  % 4 B
z(3)=0.5;                   % 2 A 
z(4)=0.75;                  % 4 B 
                
F=2*fA.*exp(1i*(Q.*z(1).*d))...
    +4*fB.*exp(1i*(Q.*z(2).*d))...
    +2*fA.*exp(1i*(Q.*z(3).*d))...
    +4*fB.*exp(1i*(Q.*z(4).*d));