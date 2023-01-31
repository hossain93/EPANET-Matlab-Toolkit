function [ pressure ] = get_pressure( DEMAND_haye_entekhabi_baraye_NN_NODE)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
d=epanet('ACOANT.inp');
ID=d.getNodeJunctionNameID;
[~,cID]=size(ID);% in alamat ~ mige oon motagheer ra hesab nakon
for i=1:cID
    % demandha bayad dar eyk satr bashad
  DEMAND_haye_entekhabi_baraye_NN_NODE(1,i)=DEMAND_haye_entekhabi_baraye_NN_NODE(i,1);
end
for i=1:cID
Demand= DEMAND_haye_entekhabi_baraye_NN_NODE(1,i);
d.setNodeBaseDemands(i,Demand);
end
d.solveCompleteHydraulics
p=d.getNodePressure;
[~,cP]=size(p);
p(:,[cID+1,cP])=[];% hazf fesharhae ezafi manand tank va reservoir va ...
pressure=p;
end
