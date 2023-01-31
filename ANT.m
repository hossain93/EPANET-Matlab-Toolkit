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
%#######################################################
minv=1;
maxv=1.5;   
%baze afzaesh debiha
increment=.01;
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

%############################################
%******************** maghadir marboot be tekrarha *******************
TEK=2000; % tedad tekrare koli
TEK_MIN=25; % TEK_MIN<TEK


%%********************************** pheromone update*******************
update_foromon_loope_asli_fit1=2*AN; %TEKRAR_e_update_shodan;% mohemtarin %%%
zaribe_update=1;
TABKHIR=.97;%.997;
fitness_e_morede_nazar=.00003*ones(1,TEK);
shorue_yek_kardane_pheromon=3200;
shorue_yek_kardane_pheromon_final=44900;
tekrari_ke_moghayese_demande_jadid_shoroo_mishavad=TEK_MIN+1;
moghayese_demandha_az_in_shomare_demand_be_bad_shoroo_shavad=TEK_MIN;


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
%zaribi ke gharare dar demandha zarb bshe


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% moshahedat %%%%%%%%%%%%%%%%%%%%%%%%
% dar inja fesharhae moshahedati grefte mishe ke mitavanad be sorat autumat
% az web GIS grefte beshe
Ho=zeros(NN,1);
HHo=zeros(NN,1);
Ho(24,1)=54.579;% or J-20
Ho(17,1)=40.15;% or J-25
Ho(30,1)=52.27;% or J-26



%##################################################################
feshar_haye_gerehyi_shabih_sazi_shode=0;
BEHTARIN_DEMAND_HA=0;
BEST_TEKRAR=1;
disp('maghadir avalea')
t=1;
BASE=first_DEMAND-4*.05*first_DEMAND;% barae shroe kar debihae avalea dade bedon mabna khasi
%meghdar demande tamame gere ha dar avalin mohasebe
DEMAND(:,t)=BASE.*ones(NN,t);
for j=1:NN
DEMAND_haye_entekhabi_baraye_NN_NODE(j,1)=DEMAND(j,t);
end
DEMAND_haye_entekhabi_baraye_NN_NODE=DEMAND_haye_entekhabi_baraye_NN_NODE;
%HSS=get_pressure( DEMAND_haye_entekhabi_baraye_NN_NODE);% matris demandha be sorat ofoghi hast in ja va dar tabe get_pressur bayad amoodi bashad



%%


matrice_Demand2=matrice_Demand;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#############################

while t<=AN-1
 for i=1:constNN2
disp('A')
 DEMAND_haye_entekhabi_baraye_NN_NODE(:,1)=matrice_Demand2(1,:);
disp('B')

 %###########################################
% agar be maximumha khord paresh mikonad
    if im(1,i)>0
f_jame_sotun_haye_ff1(t,i)=f_jame_sotun_haye_ff1(t-1,i);
DEMAND_haye_entekhabi_baraye_NN_NODE(i,1)=matrice_Demand2(1,i);     
        continue
    end
 %###########################################
 % ek halghe dige lazem darim ke ta enteha bere
 % ekhtesase debi badi ba dar nazar gereftan nodhaye hazf shode
 
 DEMAND_haye_entekhabi_baraye_NN_NODE(i,1)=matrice_Demand2(t+1,i);    
 
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

f_jame_sotun_haye_ff1(t,i)=sum(ff_moraba_e_ekhtelaf_Ho_va_HS);
   
disp('G')

%###########################################

% pida kardan minimm ff chon ke meghdarhaye hazf shode dekhalat pida
% nakonad in kar ra kardim
% ama ma maximmha ra hazf mikonim pas nbayad in halghe ra tashkil bedim???
MINIMUM_f_dar_LOOP=f_jame_sotun_haye_ff1(1,1);
disp('H')

    for j=1:t
        if im(1,i)>0
            continue
        end
if MINIMUM_f_dar_LOOP>f_jame_sotun_haye_ff1(j,i)
   MINIMUM_f_dar_LOOP=f_jame_sotun_haye_ff1(j,i); 
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
 jam_moraba_ekhtelafha=f_jame_sotun_haye_ff1;
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
[~,n]=find(f_jame_sotun_haye_ff1==ffmin);
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
  ffsumco(t,i)=sum(f_jame_sotun_haye_ff1(:,n)-f_jame_sotun_haye_ff1(:,i));
  
 for j=1:t % baraye mohasebe darsad ekhtelaf
    percentffco(j,i)=(f_jame_sotun_haye_ff1(j,n)-f_jame_sotun_haye_ff1(j,i))./f_jame_sotun_haye_ff1(j,n);
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



%****************************fromon avalea*****************************%
NN=30;
T=ones(AN,NN); % foromone avaliye barea hame debiha
alfa=1.0; %for pheromone
beta=1; % for matlubiat


%**************************** kavosh avalea *****************************%
ccc=ones(AN,NN);
matlubiat=10*AN;
for j=1:NN
for i=1:AN
if DEMAND_BASE(j,1)<=1.0001*matrice_Demand(i,j) &&...
DEMAND_BASE(j,1)>=0.9999*matrice_Demand(i,j)%dar inja gofte agar matrice demandhaee ke dar zareab zarb shode dar in bazeh bod oon drae ra barabar matlubiat gharar bede
ccc(i,j)=matlubiat;
if im(1,j)>0
ccc(1,j)=90;    
end
end
end
end
MINIMUM_f_dar_LOOP=1;
t=1;
%%
%%%%%%%%%%%%%%%%%%%%loop asli %%%%%%%%%%%%%%%%%%%%%%%%%

while MINIMUM_f_dar_LOOP > fitness_e_morede_nazar(1,t) % fitness>3e-5
for j=1:NN
  if im(1,j)>0
  DEMAND_haye_entekhabi_baraye_NN_NODE(j,1)=matrice_Demand2(1,j);    
     continue 
  end
    
    
for i=1:AN
TT(i,1)=T(i,j);
end
for i=1:AN
cc(i,1)=ccc(i,j);
end
cca=cc.^beta;
TTb=TT.^alfa;
sig01=times(cca,TTb);
sig02=sum(sig01);
for i=1:AN
p(i,j)=(ccc(i,j)^beta*T(i,j)^alfa)/(sig02);
end
%####################################
pc=cumsum(p);
x(j)=rand;
for i=1:AN-1
if 0<= x(j) && x(j) <pc(1,j)
DEMAND(j,t)=matrice_Demand(1,j);
elseif pc(i,j)<= x(j) && x(j) <pc(i+1,j)
DEMAND(j,t)=matrice_Demand(i+1,j) ;
end
end
DEMAND_haye_entekhabi_baraye_NN_NODE(j,1)=DEMAND(j,t);

end % end for j=1:NN

%*************************************************************************
% dar in ghesmat hadaf in ast ke debihaye tekrari ra barrasi nakonad va
% hengami ke tekrari=0 shavad eyani debiha tekrari nist
tekrari=0;
if t>=tekrari_ke_moghayese_demande_jadid_shoroo_mishavad % mige agar t>=TEK_MIN(=25)+1 bood mrahel zyr ra shroo kon
for kk=moghayese_demandha_az_in_shomare_demand_be_bad_shoroo_shavad:t-1 % TEK_MIN : t-1
DEMAN_t_om(:,1)=DEMAND(:,kk); % tamam demandhay soton kk ra briz dakhel ston yk oon matric dige
if DEMAND(:,t)==DEMAN_t_om(:,1) % baad myad oon demandhye (DEMAN_T_OM )ra ba akharin demand be dast amade moghayese mikone va in braye hame demandha anjam mishe va ba akharin demand moghayese mishavad
 tekrari=1; %bad agar debihaye akhari ba har kodom az ghabliha eki bood tekrari=1 eani fromon ra kahesh nade va agar barabar nabod haman sefr mimanad va fromon dar marhale badi kam mishavad
end
end
end

%tekrari=tekrari;
if tekrari==0
T=TABKHIR*T ; % tabkhire foromon
%*************************************************************************

HSS=get_pressure(DEMAND_haye_entekhabi_baraye_NN_NODE);
for j=1:NN
if HSS(1,j)<0
HSS(1,j)=0;
end
end
for j=1:NN
HS(j,1)=HSS(1,j);
kole_fesharhaye_tolid_shode(j,t)=HS(j,1);
end
for j=1:NN
if Ho(j,1)==0
HHo(j,1)=HS(j,1);
end
HHo2=HHo+Ho;
ff(j,1)=(HS(j,1)-HHo2(j,1))^2;
end
ff_moraba_e_ekhtelaf_Ho_va_HS=ff;
f_jame_sotun_haye_ff(1,t)=sum(ff_moraba_e_ekhtelaf_Ho_va_HS);
MINIMUM_f_dar_LOOP=min(f_jame_sotun_haye_ff);
fitness(1,t)=MINIMUM_f_dar_LOOP;



%############################# BEHTARIN DEMAND HA ########################

if t==shorue_yek_kardane_pheromon
for j=1:NN
if BEHTARIN_DEMAND_HA(j,1) >=1.001*DEMAND_BASE(j,1) ||...
BEHTARIN_DEMAND_HA(j,1) <=0.999*DEMAND_BASE(j,1)
BEHTARIN_DEMAND_HA222=BEHTARIN_DEMAND_HA;
T(1:AN,j)=1;
ccc(1:AN,j)=1;
end
end
disp('shorue_yek_kardane_pheromon');
end % end for:t>shorue_yek_kardane_pheromon
if t==shorue_yek_kardane_pheromon_final
for ii=1:NN
for jj=1:AN
if ccc_demand_haye_paye(ii,1)==DEMAND_BASE(jj,1)
ccc(jj,:)=matlubiat;
end
end
end
for ii2=1:NN
if BEHTARIN_DEMAND_HA(ii2,1)~=ccc_demand_haye_paye(ii2,1)
T(:,ii2)=1;
ccc(:,ii2)=1;
ccc(afzayeshe_matlubiat:AN,ii2)=meghdar_afzayeshe_matlubiat;
end
end
end % end for:t>shorue_yek_kardane_pheromon


%%%%%%% mohemtarin beroz resani fromon %%%%% %%%%%%%%%%%%%

%************ entekhab behtarin demand va feshar ************
if t==1
BEHTARIN_DEMAND_HA=DEMAND(:,t);
BEST_TEKRAR=t;
feshar_haye_gereyi_shabih_sazi_shode=HS;
else
if f_jame_sotun_haye_ff(1,t)<=MINIMUM_f_dar_LOOP ...
&& f_jame_sotun_haye_ff(1,t)~=f_jame_sotun_haye_ff(1,t-1)
BEHTARIN_DEMAND_HA=DEMAND(:,t);
BEST_TEKRAR=t;
feshar_haye_gereyi_shabih_sazi_shode=HS;
end
end
%###########################################################
%************ berozresani feromone behtarin demandha ************
for m=1:NN
for n=1:AN
if DEMAND_haye_entekhabi_baraye_NN_NODE(m,1)==matrice_Demand(n,1)
T(n,m)=T(n,m)+zaribe_update*update_foromon_loope_asli_fit1;
end
end
end

%##################################################################
disp(['tekrare= ',num2str(t)]);
disp(['fitnes= ',num2str(fitness(1,t))]);
t=t+1;
else continue % baraye if _e_ aval(tekrari=0 or 1)
end
 %end % baraye if _e_ aval(tekrari=0 or 1)
%##################################################################
if t==TEK || fitness(1,t-1)<=0.001
break
end
%##################################################################

disp('akhare while');
end %end % baraye while
