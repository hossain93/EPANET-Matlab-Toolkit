%clc;
%clear;
d=epanet('ACOANT.inp');
ID=d.getNodeJunctionNameID;
[rID,cID]=size(ID);
 for i=1:cID
    DEMAND_haye_entekhabi_baraye_NN_NODE(1,i)=DEMAND_haye_entekhabi_baraye_NN_NODE(i,1);
end

for i=1:cID
Demand= DEMAND_haye_entekhabi_baraye_NN_NODE(1,i);
d.setNodeBaseDemands(i,Demand);
end
d.solveCompleteHydraulics
P=d.getNodePressure;
[rP,cP]=size(P);
 P(:,[cID+1,cP])=[];% hazf fesharhae ezafi manand tank va reservoir
  h=d.getNodeBaseDemands;