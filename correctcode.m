%leakage in node 21 or J-19.........leak=25 l/s

start_toolkit()
clc;
clear all;
format short g
time1=clock; % clock = [year month day hour minute seconds]
hr1= time1(1,4);
min1=time1(1,5);
sec11=time1(1,6);
sec1=60-sec11;
%in bakhsh faghat zaman ra midahad ke saat chand ast




%#######################################################
minv=1;
maxv=1.5;   
%baze afzaesh debiha
increment=.1;
AN=(maxv-minv)/increment+1; % AN bayad zoj bashad-- in bakhsh tedad bazeha afzaesh debi ra midahad
NN=30;%tedad gerehha ya tedad debiha
constNN=NN;
constNN2=NN;
iNN=1;
t=1;
im=zeros(1,NN);
tm=zeros(1,NN);
in=zeros(AN,NN);
b=0;
cNN=0;
%percentffco=zeros(AN-1,NN);


%######################################################
% masaref bardashti az kountourha %
%$$####################################################
DEMAND_BASE=zeros(NN,1);
first_DEMAND=50; % first_DEMAND= avalin meghdare DEMAND
DEMAND_BASE(1,1)=first_DEMAND; % avalin meghdare DEMAND
step_DEMAND=.0; % step_c= mizan afzayeshe c
for j=2:NN
DEMAND_BASE(j,1)=DEMAND_BASE(j-1,1)+ step_DEMAND ;
end
ccc_demand_haye_paye=DEMAND_BASE;




%################# baze taghirat demand ################
Bazeye_zaribe_demand_ha=zeros(AN,1);
Bazeye_zaribe_demand_ha(1,1)=minv;%zarib avali eak ast
for i=2:AN
Bazeye_zaribe_demand_ha(i,1)=Bazeye_zaribe_demand_ha(i-1,1)+ increment ;
end
Bazeye_zaribe_demand_ha=Bazeye_zaribe_demand_ha;%zaribi ke gharare dar demandha zarb bshe




%############### kol entekhabhae mojod barae demand ################
matrice_Demand=zeros(AN,NN);
for j=1:NN
for i=1:AN
if DEMAND_BASE(j,1)~=0%mige agar demand sefr nabod bro zaribsho pida kon va zarb in debi kon
matrice_Demand(i,j)=DEMAND_BASE(j,1)*Bazeye_zaribe_demand_ha(i,1);
else matrice_Demand(i,j)=...% va dar ghir in sorat hamon demand ba zaribesh menhae kamtrin zarib kon
DEMAND_BASE(j,1)+Bazeye_zaribe_demand_ha(i,1)-minv;
end
end
end
matrice_Demand=matrice_Demand;
matrice_Demand2=matrice_Demand;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% moshahedat %%%%%%%%%%%%%%%%%%%%%%%%
% dar inja fesharhae moshahedati grefte mishe ke mitavanad be sorat autumat
% az web GIS grefte beshe
Ho=zeros(NN,1);
HHo=zeros(NN,1);
Ho(24,1)=54.579;% or J-20
Ho(17,1)=40.15;% or J-25
Ho(30,1)=52.27;% or J-26

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#############################

while t<=2
 for i=1:constNN2
disp('A')
 DEMAND_haye_entekhabi_baraye_NN_NODE(:,1)=matrice_Demand2(1,:);
disp('B')
    if im(1,s)>0
DEMAND_haye_entekhabi_baraye_NN_NODE(s,1)=matrice_Demand2(1,s);     
        continue
    end
 %###########################################
 
 %hazfe nodhaye ba ekhtelaf zeyad
 %matrice_Demand=matrice_Demand2;
%c=1;
 %for j=1:constNN
% if im(1,j)~=0
 %DEMAND_haye_entekhabi_baraye_NN_NODE(c,:)=[];
  %matrice_Demand(:,j)=[];
 % c=c-1;
 %constNN=constNN-1;
 %end 
 %c=c+1;
 %end
 %constNN=constNN2;
 %disp('C')
 
 %###########################################
 
 %jabejaee nodhaye hazf shode va ekhtesas debi nodhaye hazf shode
%c=NN;
 %for p=1:constNN
 %if im(1,p)>0 
 %for tr=p:c
 %DEMAND_haye_entekhabi_baraye_NN_NODE(tr+1,1)=DEMAND_haye_entekhabi_baraye_NN_NODE(tr,1);       
 
 %end
 %c=c+1;
 % DEMAND_haye_entekhabi_baraye_NN_NODE(p,1)=matrice_Demand2(1,p);
 %end
 %end
 disp('D')
 %###########################################
 % ek halghe dige lazem darim ke ta enteha bere
 % ekhtesase debi badi ba dar nazar gereftan nodhaye hazf shode
 c=iNN;
 if im(1,c)>0 
  c=c+1;
 DEMAND_haye_entekhabi_baraye_NN_NODE(c,1)=matrice_Demand2(t+1,c);
 else
 DEMAND_haye_entekhabi_baraye_NN_NODE(c,1)=matrice_Demand2(t+1,c);    
 end
disp('F')

  %###########################################
 
HSS=get_pressure( DEMAND_haye_entekhabi_baraye_NN_NODE);
for j=1:constNN
if HSS(1,j)<0
HSS(1,j)=0;
end
end
for j=1:constNN
HS(j,1)=HSS(1,j);
kole_fesharhaye_tolid_shode(j,t)=HS(j,1);
end
for j=1:constNN
if Ho(j,1)==0
HHo(j,1)=HS(j,1);
end
HHo2=HHo+Ho;
ff(j,1)=(HS(j,1)-HHo2(j,1))^2;
end
ff_moraba_e_ekhtelaf_Ho_va_HS=ff;

tNN=iNN;

if im(1,iNN)>0
    b=b+1;
tNN=tNN+b;
cNN=tNN;
f_jame_sotun_haye_ff(t,tNN+1)=sum(ff_moraba_e_ekhtelaf_Ho_va_HS);
else
 f_jame_sotun_haye_ff(t,cNN)=sum(ff_moraba_e_ekhtelaf_Ho_va_HS);   
end
disp('G')

%###########################################

% pida kardan minimm ff
MINIMUM_f_dar_LOOP=f_jame_sotun_haye_ff(1,1);
for j=1:iNN
    for i=1:t
if MINIMUM_f_dar_LOOP>f_jame_sotun_haye_ff(i,j)
   MINIMUM_f_dar_LOOP=f_jame_sotun_haye_ff(i,j); 
end
 if im(1,j)>0
     continue
 end
    
    
    end
end

disp('H')
  %###########################################
  %namayesh gerafiki
%fitness(1,t)=MINIMUM_f_dar_LOOP;
disp(['tekrare= ',num2str(t),'   and   ','payan','----->','sooton= ',num2str(iNN)])
iNN=iNN+1;
disp(['shrooe','----->','sooton= ',num2str(iNN)])
 end 
 
 %###########################################
 disp('I')
 % dar inja stoon minimm moraba fesharha ra bedast my avarim
 % baraye in ke stoonhaye hazf shode tasir ndashte bashad anha ra hazf
 % mikonim
 jam_moraba_ekhtelafha=f_jame_sotun_haye_ff;
 c=1;
 for j=1:constNN
 if im(1,j)>0
 jam_moraba_ekhtelafha(:,c)=[];
 c=c-1;
 constNN=constNN-1;
 end 
 c=c+1;
 end
 c=1;
 disp('J')
 constNN=constNN2;
 %f=0;
ffmin=min(min(jam_moraba_ekhtelafha));
[~,n]=find(f_jame_sotun_haye_ff==ffmin);
disp('K')
in(t,n)=n;
%for j=1:constNN
  % if j>=n
    %  n=n+1+f;
    %  break
  % end
  % if im(1,j)>0 
  %   f=f+1;  
 %  end
%end
%disp('L')
%###########################################

for i=1:constNN % in halghe mohasebe ekhtelaf moraba ekhtelafha va darsad ekhtelaf anha ra anjam midahad
  ffsumco(t,i)=sum(f_jame_sotun_haye_ff(:,n)-f_jame_sotun_haye_ff(:,i));
  
 for j=1:t % baraye mohasebe darsad ekhtelaf
    percentffco(j,i)=(f_jame_sotun_haye_ff(j,n)-f_jame_sotun_haye_ff(j,i))./f_jame_sotun_haye_ff(j,n);
 end
  
end
disp('M')
%#############################

for j=1:t% in halghe baraye mosbat kardan ekhtelaf darsadha ast
    for i=1:constNN
      if percentffco(j,i)<0
       percentffco(j,i)=percentffco(j,i)*(-1);
      end
    end
end
disp('N')
%#############################

for i=1:constNN
percentffsumco(1,i)=sum( percentffco(:,i));
end
disp('O')
%#############################

percentffsumco2=percentffsumco;
c=1;
 for j=1:constNN
 if im(1,j)>0
 percentffsumco2(:,c)=[];
 c=c-1;
 constNN=constNN-1;
 end  
 c=c+1;
 end
 disp('P')
 constNN=constNN2;
percentff=max(max(percentffsumco2));
[~,m]=find(percentffsumco==percentff);
c=1;
%disp('Q')
%f=0;
%for j=1:constNN
    %if j>=m
    %    m=m+1+f; 
    %    break
   % end      
   %if tm(1,j)>0 
 %    f=f+1;
 %  end
%end
im(1,m)=m;
tm(1,:)=im(1,:);
NN=NN-1;
iNN=1;
disp('R')
 %###########################################
 
 
%disp(['tekrare= ',num2str(t)])% dar disp bayad hame reshte bashad va tabe num2str hamin kar ra mikonad
disp(['akhare while for tekrare= ',num2str(t)])
t=t+1;
disp(['shrooe tekrare= ',num2str(t)])

end




%%
 kamtrinffco=zeros(AN-1,NN);
 percentffsumco=zeros(AN-1,NN);
ffmin=min(min(f_jame_sotun_haye_ff));
[~,n]=find(f_jame_sotun_haye_ff==ffmin);

for i=1:NN % in halghe mohasebe ekhtelaf moraba ekhtelafha va darsad ekhtelaf anha ra anjam midahad
  ffsumco(1,i)=sum(f_jame_sotun_haye_ff(:,n)-f_jame_sotun_haye_ff(:,i));
   for j=1:AN-1% baraye mohasebe darsad ekhtelaf
    percentffsumco(j,i)=(f_jame_sotun_haye_ff(j,n)-f_jame_sotun_haye_ff(j,i))./f_jame_sotun_haye_ff(j,n);
   end
  
end
%#############################
for j=1:AN-1% in halghe baraye mosbat kardan ekhtelaf darsadha ast
    for i=1:NN
      if percentffsumco(j,i)<0
       percentffsumco(j,i)=percentffsumco(j,i)*(-1);
      end
    end
end
%#############################


kamtrinffco(:,n)=f_jame_sotun_haye_ff(:,n);
%%
for it=1:NN
  if ffsumco(1,it)>=ffsumco(1,n)
      
   kamtrinff(:,it)=f_jame_sotun_haye_ff(:,it);   
    
  end 
end

%%
g=1:1:NN;
legendCell = cellstr(num2str(g', 'J-%-d'));
fg=matrice_Demand;
fg(1,:)=[];
plot(fg,kamtrinff)
legend(legendCell)
%%
hold on
plot(fg(:,1),f_jame_sotun_haye_ff(:,21),'*')
 



