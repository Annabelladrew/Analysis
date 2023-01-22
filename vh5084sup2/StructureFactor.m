% Created by Celine Lichtensteiger
% Calculates the structure factor for materials with a perovskite structure
%*********************************
function[F]=StructureFactor(fA,fB,fO,Q,z,c)
F=fA.*exp(1i*(Q.*z(1).*c))+fB.*exp(1i*(Q.*z(2).*c))+fO.*exp(1i*(Q.*z(3).*c))+fO.*exp(1i*(Q.*z(4).*c))+fO.*exp(1i*(Q.*z(5).*c));