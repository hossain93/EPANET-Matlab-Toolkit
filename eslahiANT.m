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


%%

%#######################################################
minv=1;
maxv=1.5;   
%baze afzaesh debiha
increment=.1;
AN=(maxv-minv)/increment+1; % AN bayad zoj bashad-- in bakhsh tedad bazeha afzaesh debi ra midahad
NN=30;%tedad gerehha ya tedad debiha


%%

%############################################
%******************** maghadir marboot be tekrarha *******************
TEK=2000; % tedad tekrare koli
TEK_MIN=25; % TEK_MIN<TEK


%%

%%********************************** pheromone update*******************
update_foromon_loope_asli_fit1=2*AN; %TEKRAR_e_update_shodan;% mohemtarin %%%
zaribe_update=1;
TABKHIR=.99;%.997;
fitness_e_morede_nazar=.00003*ones(1,TEK);
shorue_yek_kardane_pheromon=3200;
shorue_yek_kardane_pheromon_final=44900;
tekrari_ke_moghayese_demande_jadid_shoroo_mishavad=TEK_MIN+1;
moghayese_demandha_az_in_shomare_demand_be_bad_shoroo_shavad=TEK_MIN;


%%

%######################################################
% masaref bardashti az kountourha %
%$$####################################################
DEMAND_BASE=zeros(1,NN);
first_DEMAND=50; % first_DEMAND= avalin meghdare DEMAND
DEMAND_BASE(1,1)=first_DEMAND; % avalin meghdare DEMAND
step_DEMAND=.0; % step_c= mizan afzayeshe c
for j=2:NN
DEMAND_BASE(1,j)=DEMAND_BASE(1,j-1)+ step_DEMAND ;
end
ccc_demand_haye_paye=DEMAND_BASE;


%%

%################# baze taghirat demand ################
Bazeye_zaribe_demand_ha=zeros(AN,1);
Bazeye_zaribe_demand_ha(1,1)=minv;%zarib avali eak ast
for i=2:AN
Bazeye_zaribe_demand_ha(i,1)=Bazeye_zaribe_demand_ha(i-1,1)+ increment ;
end
%zaribi ke gharare dar demandha zarb bshe


%%

%############### kol entekhabhae mojod barae demand ################
matrice_Demand=zeros(AN,NN);
for j=1:NN
for i=1:AN
if DEMAND_BASE(1,j)~=0%mige agar demand sefr nabod bro zaribsho pida kon va zarb in debi kon
matrice_Demand(i,j)=DEMAND_BASE(1,j)*Bazeye_zaribe_demand_ha(i,1);
else matrice_Demand(i,j)=...% va dar ghir in sorat hamon demand ba zaribesh menhae kamtrin zarib kon
DEMAND_BASE(1,j)+Bazeye_zaribe_demand_ha(i,1)-minv;
end
end
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% moshahedat %%%%%%%%%%%%%%%%%%%%%%%%
% dar inja fesharhae moshahedati grefte mishe ke mitavanad be sorat autumat
% az web GIS grefte beshe
Ho=zeros(1,NN);
HHo=zeros(1,NN);
Ho(1,24)=54.579;% or J-20
Ho(1,17)=40.15;% or J-25
Ho(1,30)=52.27;% or J-26


%%
%****************************fromon avalea*****************************%
T=ones(AN,NN); % foromone avaliye barea hame debiha
alfa=1; %for pheromone
beta=1; % for matlubiat


%%

%**************************** kavosh avalea *****************************%
ccc=ones(AN,NN);
matlubiat=10*AN;
for j=1:NN
for i=1:AN
if DEMAND_BASE(1,j)<=1.0001*matrice_Demand(i,j) &&...
DEMAND_BASE(1,j)>=0.9999*matrice_Demand(i,j)%dar inja gofte agar matrice demandhaee ke dar zareab zarb shode dar in bazeh bod oon drae ra barabar matlubiat gharar bede
ccc(i,j)=matlubiat;
end
end
end


%%
%##################################################################
feshar_haye_gerehyi_shabih_sazi_shode=0;
BEHTARIN_DEMAND_HA=0;
BEST_TEKRAR=1;
disp('maghadir avalea')
t=1;
BASE=first_DEMAND-4*.05*first_DEMAND;% barae shroe kar debihae avalea dade bedon mabna khasi
%meghdar demande tamame gere ha dar avalin mohasebe
DEMAND(t,:)=BASE.*ones(t,NN);
for j=1:NN
DEMAND_haye_entekhabi_baraye_NN_NODE(1,j)=DEMAND(t,j);
end

%HSS=get_pressure( DEMAND_haye_entekhabi_baraye_NN_NODE);% matris demandha be sorat ofoghi hast in ja va dar tabe get_pressur bayad amoodi bashad


%%
%%%%%%%%%%%%%%%%%%%%   loop pida kardan minimm morabae fesharha %%%%%%%%%%

while t<TEK_MIN
for j=1:NN
    
 %################################## 
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
%##################################
pc=cumsum(p);
x(j)=rand;
for i=1:AN-1
if 0<= x(j) && x(j) <pc(1,j)
DEMAND(t,j)=matrice_Demand(1,j); %chon darsad kavoshi ra zeyad greftim
elseif pc(i,j)<= x(j) && x(j) <pc(i+1,j)
DEMAND(t,j)=matrice_Demand(i+1,j) ;
end
end
DEMAND_haye_entekhabi_baraye_NN_NODE(1,j)=DEMAND(t,j);
%##################################
tekrarhaemin(1,j)=t;
end % end for j=1:NN

HSS=get_pressure( DEMAND_haye_entekhabi_baraye_NN_NODE);
for j=1:NN
if HSS(1,j)<0
HSS(1,j)=0;
end
end

for j=1:NN
HS(1,j)=HSS(1,j);
kole_fesharhaye_tolid_shode(t,j)=HS(1,j);
end
for j=1:NN
if Ho(1,j)==0
HHo(1,j)=HS(1,j);
end
HHo2=HHo+Ho;
ff(1,j)=(HS(1,j)-HHo2(1,j))^2;
end
ff_moraba_e_ekhtelaf_Ho_va_HS=ff;
f_jame_sotun_haye_ff(1,t)=sum(ff_moraba_e_ekhtelaf_Ho_va_HS);
MINIMUM_f_dar_LOOP=min(f_jame_sotun_haye_ff);
fitness(1,t)=MINIMUM_f_dar_LOOP;
disp(['tekrare= ',num2str(t)])% dar disp bayad hame reshte bashad va tabe num2str hamin kar ra mikonad
t=t+1;
disp('akhare while for min')
end
fitness2(1,:)=fitness(1,t-1);

%%
%%%%%%%%%%%%%%%%%%%%loop asli %%%%%%%%%%%%%%%%%%%%%%%%%

while t<90
    %fitness(1,t-1) > fitness_e_morede_nazar(1,t-1) % fitness>3e-5
for j=1:NN
 
    
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
DEMAND(t,j)=matrice_Demand(1,j);
elseif pc(i,j)<= x(j) && x(j) <pc(i+1,j)
DEMAND(t,j)=matrice_Demand(i+1,j) ;
end
end
DEMAND_haye_entekhabi_baraye_NN_NODE(1,j)=DEMAND(t,j);

end % end for j=1:NN

%*************************************************************************
% dar in ghesmat hadaf in ast ke debihaye tekrari ra barrasi nakonad va
% hengami ke tekrari=0 shavad eyani debiha tekrari nist
tekrari=0;
if t>=tekrari_ke_moghayese_demande_jadid_shoroo_mishavad % mige agar t>=TEK_MIN(=25)+1 bood mrahel zyr ra shroo kon
for kk=moghayese_demandha_az_in_shomare_demand_be_bad_shoroo_shavad:t-1 % TEK_MIN : t-1
DEMAN_t_om(1,:)=DEMAND(kk,:); % tamam demandhay soton kk ra briz dakhel ston yk oon matric dige
if DEMAND(t,:)==DEMAN_t_om(1,:) % baad myad oon demandhye (DEMAN_T_OM )ra ba akharin demand be dast amade moghayese mikone va in braye hame demandha anjam mishe va ba akharin demand moghayese mishavad
 tekrari=1; %bad agar debihaye akhari ba har kodom az ghabliha eki bood tekrari=1 eani fromon ra kahesh nade va agar barabar nabod haman sefr mimanad va fromon dar marhale badi kam mishavad
end
end
end
T=TABKHIR*T ; % tabkhire foromon
%tekrari=tekrari;
if tekrari==0

%*************************************************************************

HSS=get_pressure(DEMAND_haye_entekhabi_baraye_NN_NODE);
for j=1:NN
if HSS(1,j)<0
HSS(1,j)=0;
end
end

for j=1:NN
HS(1,j)=HSS(1,j);
kole_fesharhaye_tolid_shode(t,j)=HS(1,j);
end

for j=1:NN
if Ho(1,j)==0
HHo(1,j)=HS(1,j);
end
HHo2=HHo+Ho;
ff(1,j)=(HS(1,j)-HHo2(1,j))^2;
end

ff_moraba_e_ekhtelaf_Ho_va_HS=ff;
f_jame_sotun_haye_ff(1,t)=sum(ff_moraba_e_ekhtelaf_Ho_va_HS);
MINIMUM_f_dar_LOOP=min(f_jame_sotun_haye_ff);
fitness(1,t)=MINIMUM_f_dar_LOOP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% mohasebat fitness ha %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

sorte_fitness_ha=sort(f_jame_sotun_haye_ff);
fitness1(1,1)=sorte_fitness_ha(1,1);
fitness2(1,1)=sorte_fitness_ha(1,2);
fitness3(1,1)=sorte_fitness_ha(1,3);
fitness4(1,1)=sorte_fitness_ha(1,4);
fitness5(1,1)=sorte_fitness_ha(1,5);
fitness6(1,1)=sorte_fitness_ha(1,6);
fitness7(1,1)=sorte_fitness_ha(1,7);


%############################# BEHTARIN DEMAND HA ########################

if f_jame_sotun_haye_ff(1,t)==fitness2(1,1)
BEHTARIN_DEMAND_HA_2=DEMAND(t,:);
end
if f_jame_sotun_haye_ff(1,t)==fitness3(1,1)
BEHTARIN_DEMAND_HA_3=DEMAND(t,:);
end
if f_jame_sotun_haye_ff(1,t)==fitness4(1,1)
BEHTARIN_DEMAND_HA_4=DEMAND(t,:);
end
if f_jame_sotun_haye_ff(1,t)==fitness5(1,1)
BEHTARIN_DEMAND_HA_5=DEMAND(t,:);
end
if f_jame_sotun_haye_ff(1,t)==fitness6(1,1)
BEHTARIN_DEMAND_HA_6=DEMAND(t,:);
end

%########################################
% chera ek mikone?????????????????
if t==shorue_yek_kardane_pheromon%=3200
for j=1:NN
if BEHTARIN_DEMAND_HA(1,j) >=1.001*DEMAND_BASE(1,j) ||...
BEHTARIN_DEMAND_HA(1,j) <=0.999*DEMAND_BASE(1,j)
BEHTARIN_DEMAND_HA222=BEHTARIN_DEMAND_HA;
T(1:AN,j)=1;
ccc(1:AN,j)=1;
end
end
disp('shorue_yek_kardane_pheromon');
end % end for:t>shorue_yek_kardane_pheromon

%########################################
% chera ek mikone?????????????????
if t==shorue_yek_kardane_pheromon_final%=44900
for ii=1:NN
for jj=1:AN
if ccc_demand_haye_paye(1,ii)==DEMAND_BASE(jj,1)
ccc(jj,:)=matlubiat;
end
end
end

for ii2=1:NN
if BEHTARIN_DEMAND_HA(1,ii2)~=ccc_demand_haye_paye(1,ii2)
T(:,ii2)=1;
ccc(:,ii2)=1;
ccc(afzayeshe_matlubiat:AN,ii2)=meghdar_afzayeshe_matlubiat;
end
end
end % end for:t>shorue_yek_kardane_pheromon


%%%%%%% mohemtarin beroz resani fromon %%%%% %%%%%%%%%%%%%

%************ entekhab behtarin demand va feshar ************
if f_jame_sotun_haye_ff(1,t)<=MINIMUM_f_dar_LOOP ...
&& f_jame_sotun_haye_ff(1,t)~=f_jame_sotun_haye_ff(1,t-1)
BEHTARIN_DEMAND_HA=DEMAND(t,:);
BEST_TEKRAR=t;
feshar_haye_gereyi_shabih_sazi_shode=HS;
%###########################################################
%************ berozresani feromone behtarin demandha ************
%inja mige kodom zarib demand va baraye kodom node boode ke monjar be
%entekhab khobi shode va bro foromon hamon zarib demand ra braye node
%moshakhas afzayesh bede
for m=1:NN
for n=1:AN
if DEMAND_haye_entekhabi_baraye_NN_NODE(1,m)==matrice_Demand(n,m)
T(n,m)=T(n,m)+zaribe_update*update_foromon_loope_asli_fit1;
end
end
end
end
%##################################################################
disp(['tekrare= ',num2str(t)]);
disp(['fitnes= ',num2str(fitness(1,t))]);
t=t+1;
else continue % baraye if _e_ aval(tekrari=0 or 1)
end %end % baraye if _e_ aval(tekrari=0 or 1)
%##################################################################
if t==TEK || fitness(1,t-1)<=0.001
break
end
%##################################################################

disp('akhare while');
end %end % baraye while
