


d = epanet('ACOANT.inp');
d.solveCompleteHydraulics
ID=d.getNodeJunctionNameID;
[rID,cID]=size(ID);
P=d.getNodePressure;
[rP,cP]=size(P);
 P(:,[cID+1,cP])=[];% hazf fesharhae ezafi manand tank va reservoir
 %ID=strjoin(ID);% tabdil cell be double
 %P=mat2cell(P);
  pressure=cell(cID,2);
for i=1:cID
   pressure{i,2}=P(i);
   pressure{i,1}=ID{1,i};
end