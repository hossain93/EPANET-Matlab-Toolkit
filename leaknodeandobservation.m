%clc;
%clear;
d=epanet('ACOANT.inp');
ID=d.getNodeJunctionNameID;
[rID,cID]=size(ID);

Demand=75 ;
d.setNodeBaseDemands(21,Demand);
d.solveCompleteHydraulics
Pe=d.getNodePressure;
[rP,cP]=size(Pe);
 Pe(:,[cID+1,cP])=[];% hazf fesharhae ezafi manand tank va reservoir
  h=d.getNodeBaseDemands;
    pressure=cell(cID,2);
for i=1:cID
   pressure{i,2}=Pe(i);
   pressure{i,1}=ID{1,i};
end