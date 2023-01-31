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
TABKHIR=.9998;%.997;
fitness_e_morede_nazar=.00003*ones(1,TEK);
shorue_yek_kardane_pheromon=3200;
shorue_yek_kardane_pheromon_final=44900;
tekrari_ke_moghayese_demande_jadid_shoroo_mishavad=TEK_MIN+1;
moghayese_demandha_az_in_shomare_demand_be_bad_shoroo_shavad=TEK_MIN;


%%

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


%%

%################# baze taghirat demand ################
Bazeye_zaribe_demand_ha=zeros(AN,1);
Bazeye_zaribe_demand_ha(1,1)=minv;%zarib avali eak ast
for i=2:AN
Bazeye_zaribe_demand_ha(i,1)=Bazeye_zaribe_demand_ha(i-1,1)+ increment ;
end
Bazeye_zaribe_demand_ha=Bazeye_zaribe_demand_ha;%zaribi ke gharare dar demandha zarb bshe


%%

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


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% moshahedat %%%%%%%%%%%%%%%%%%%%%%%%
% dar inja fesharhae moshahedati grefte mishe ke mitavanad be sorat autumat
% az web GIS grefte beshe
Ho=zeros(NN,1);
HHo=zeros(NN,1);
Ho(24,1)=26.40;
Ho(17,1)=30.75;
Ho(30,1)=26.34;


%%
%****************************fromon avalea*****************************%
T=ones(AN,NN); % foromone avaliye barea hame debiha
alfa=1.0; %for pheromone
beta=1; % for matlubiat


%%

%**************************** kavosh avalea *****************************%
ccc=ones(AN,NN);
matlubiat=10*AN;
for j=1:NN
for i=1:AN
if DEMAND_BASE(j,1)<=1.0001*matrice_Demand(i,j) &&...
DEMAND_BASE(j,1)>=0.9999*matrice_Demand(i,j)%dar inja gofte agar matrice demandhaee ke dar zareab zarb shode dar in bazeh bod oon drae ra barabar matlubiat gharar bede
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
DEMAND(:,t)=BASE.*ones(NN,t);
for j=1:NN
DEMAND_haye_entekhabi_baraye_NN_NODE(j,1)=DEMAND(j,t);
end
DEMAND_haye_entekhabi_baraye_NN_NODE=DEMAND_haye_entekhabi_baraye_NN_NODE;
%HSS=DEpanet1(DEMAND_haye_entekhabi_baraye_NN_NODE,NN)% ehtemalan gofte bea fsharha ro ba in debiha bro to barname mohasebe kon va pas bear
%felan kode bala ra khareg kardim va debiha ra dasti vared mikonim
HSS=zeros(1,NN);
for j=1:NN
    HSS(1,j)=35*rand(1,1);
end
%dar inja adad tasadofi tolid kardim barea HSS


%%
HSS_AVALIYE=HSS;
for j=1:NN
if HSS(1,j)<0
HSS(1,j)=0;% dar in halghe gofte oon fesharhae ke manfi ast ra sefr gharar bede
end
end
for j=1:NN
HS(j,1)=HSS(1,j);
kole_fesharhaye_tolid_shode(j,t)=HS(j,1);% dar in halghe oomade fesharha ra az halat radifi be sotoni tabdil karde
end
HS_Mohasebati=HS;
for j=1:NN
if Ho(j,1)==0% mige agar matris Ho ke ghablan fesharhae moshahedati ra dar oon gharar dadim agar deraee az an sefr bod bead feshar mohasebati ra dar HHo gharar bede ke ghablan sefr bod hamash
HHo(j,1)=HS(j,1);
end
HHo2=HHo+Ho;% dar inja HHo2 ham feshar mohasebati ra darad va ham moshahedati ra
ff(j,1)=(HS(j,1)-HHo2(j,1))^2;% moraba ekhtelafat ra mohasebe mikone ke bebine aya fesharhae mohasebati va moshahedati eki shode ya na
end
ff_moraba_e_ekhtelaf_Ho_va_HS=ff;
f_jame_sotun_haye_ff(1,1)=sum(ff_moraba_e_ekhtelaf_Ho_va_HS);% moraba ekhtlafhae mohasebe shode barae hame fesharha ra ba ham jam mikonad
fitness(1,1)=f_jame_sotun_haye_ff(1,1);


%%

disp(' tekrar dovom : t= 2 ')
t=2;
for j=1:1
for i=1:AN
TT(i,1)=T(i,j);%fromon
cc(i,1)=ccc(i,j);% in ke har debi cheghadr arzesh entekhab darad
end
cca=cc.^beta;
TTb=TT.^alfa;
sig01=times(cca,TTb);
sig02=sum(sig01);
for i=1:AN
p(i,j)=(ccc(i,j)^beta*T(i,j)^alfa)/(sig02);%chon ke arzesh va fromon ra brabar grftim brae bar aval ehtemal entekhab hame barabar ast
end
pc=cumsum(p);% be soorat tajamoee jam mikonad
x(j)=rand;
for i=1:AN-1
    r(j)=rand;
    w=cumsum(p);
   k=find(r(j)<=w,1,'first');
   DEMAND(j,t)=matrice_Demand(k,j);
end
%%
for i=1:AN-1
if 0<= x(j) && x(j)<pc(1,j)% mige agar x(j) bin in baze bod az matrise zarieb zarb shode dar debiha avalin satr ra entekhab kon
DEMAND(j,t)=matrice_Demand(1,j);
elseif pc(i,j)<= x(j) && x(j)<pc(i+1,j) % in mige bro satrhae dige ra entekhab kon
DEMAND(j,t)=matrice_Demand(i+1,j);
end
end
DEMAND_haye_entekhabi_baraye_NN_NODE(j,1)=DEMAND(j,t);
end % end for j=1:NN
% barae bar aval myad oon debihaee ke bedon mabna entekhab karde bood ra
% eslah mikond va az dakhel matris zaraeb dar debiha entekhab mikonad
%va dar marahel bad ba estefade az arzesh va fromon behbood mibakhshim
%debiha ra ke fekr konam inja anjam nadade, barae in kar vaghti debihae bala
%bedast amad fesharha va moraba va majmoe morabe anha ra hesab karde va
%sepas makoos karde ta arzesh har kodam bedast byad


%%
%va dobare omade bar asas debehae ghabli ke bedast omade fesharha ra
%mohasebe karde
HSS=DEpanet1(DEMAND_haye_entekhabi_baraye_NN_NODE,NN);
for j=1:NN
if HSS(1,j)<0
HSS(1,j)=0
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
fitness(1,2)=f_jame_sotun_haye_ff(1,t);
T_avaliye=T;%%feromon
t=3


%%
%%%%%%%%%%%%%%%%%%%%%%%%%   loop pida kardan minimm  %%%%%%%%%%%%%%%%%%%%%%%%%

while t<TEK_MIN
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
pc=cumsum(p);
x(j)=rand;
for i=1:AN-1
if 0<= x(j) && x(j) <pc(1,j)
DEMAND(j,t)=matrice_Demand(1,j);
elseif pc(i,j)<= x(j) && x(j) <pc(i+1,j)
DEMAND(j,t)=matrice_Demand(i+1,j) ;
end
end
DEMAND_haye_entekhabi_baraye_NN_NODE(j,1)=DEMAND(j,t) ;
end % end for j=1:NN
HSS=DEpanet1(DEMAND_haye_entekhabi_baraye_NN_NODE,NN);
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
t=t+1
disp('akhare while for min')
end
fitness2(1,:)=fitness(1,t-1)


%%
%%%%%%%%%%%%%%%%%%%%loop asli %%%%%%%%%%%%%%%%%%%%%%%%%

while fitness(1,t-1) > fitness_e_morede_nazar(1,t-1)
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
pc=cumsum(p);
x(j)=rand;
for i=1:AN-1
if 0<= x(j) && x(j) <pc(1,j)
DEMAND(j,t)=matrice_Demand(1,j);
elseif pc(i,j)<= x(j) && x(j) <pc(i+1,j)
DEMAND(j,t)=matrice_Demand(i+1,j) ;
end
end
DEMAND_haye_entekhabi_baraye_NN_NODE(j,1)=DEMAND(j,t) ;
end % end for j=1:NN
tekrari=0;
if t>=tekrari_ke_moghayese_demande_jadid_shoroo_mishavad
for kk=moghayese_demandha_az_in_shomare_demand_be_bad_shoroo_shavad:t-1
DEMAN_t_om(:,1)=DEMAND(:,kk);
if DEMAND(:,t)==DEMAN_t_om(:,1)
tekrari=1;
end
end
end
tekrari=tekrari;
%%
%*************************************************************************
%*************************************************************************
if tekrari==0
T=TABKHIR*T ; % tabkhire foromon
HSS=DEpanet1(DEMAND_haye_entekhabi_baraye_NN_NODE,NN);
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% mohasebat fitness ha
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
sorte_fitness_ha=sort(f_jame_sotun_haye_ff);
fitness1(1,1)=sorte_fitness_ha(1,1)
fitness2(1,1)=sorte_fitness_ha(1,2)
fitness3(1,1)=sorte_fitness_ha(1,3)
fitness4(1,1)=sorte_fitness_ha(1,4)
fitness5(1,1)=sorte_fitness_ha(1,5);
fitness6(1,1)=sorte_fitness_ha(1,6);
fitness7(1,1)=sorte_fitness_ha(1,7);
%%
%############################# BEHTARIN DEMAND HA ########################
if f_jame_sotun_haye_ff(1,t)==fitness2(1,1)
BEHTARIN_DEMAND_HA_2=DEMAND(:,t);
end
if f_jame_sotun_haye_ff(1,t)==fitness3(1,1)
BEHTARIN_DEMAND_HA_3=DEMAND(:,t);
end
if f_jame_sotun_haye_ff(1,t)==fitness4(1,1)
BEHTARIN_DEMAND_HA_4=DEMAND(:,t);
end
if f_jame_sotun_haye_ff(1,t)==fitness5(1,1)
BEHTARIN_DEMAND_HA_5=DEMAND(:,t);
end
if f_jame_sotun_haye_ff(1,t)==fitness6(1,1)
BEHTARIN_DEMAND_HA_6=DEMAND(:,t);
end
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
%%
%%%%%%% mohemtarin beroz resani fromon %%%%% %%%%%%%%%%%%%
if f_jame_sotun_haye_ff(1,t)<=MINIMUM_f_dar_LOOP ...
&& f_jame_sotun_haye_ff(1,t)~=f_jame_sotun_haye_ff(1,t-1)
BEHTARIN_DEMAND_HA=DEMAND(:,t);BEST_TEKRAR=t;
feshar_haye_gereyi_shabih_sazi_shode=HS;
for m=1:NN
for n=1:AN
if DEMAND_haye_entekhabi_baraye_NN_NODE(m,1)==matrice_Demand(n,1)
T(n,m)=T(n,m)+zaribe_update*update_foromon_loope_asli_fit1;
end
end
end
end
%%
%##################################################################
if t==TEK
break
end
%%
%##################################################################
t=t+1
else continue % baraye if _e_ aval(tekrari=0 or 1)
end %end % baraye if _e_ aval(tekrari=0 or 1)
disp('akhare while');
end %end % baraye while
%%
%%%%%%%%%%%%%%%%%%%%%%%%entehae mohasebat %%%%%%%%%%%%%%%%%%%%%%%%%

kole_DEMAND_ha=DEMAND;
kole_fesharhaye_tolid_shode=kole_fesharhaye_tolid_shode;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
kolle_fitness_haye_ijad_shode=f_jame_sotun_haye_ff
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp ('fromon nahaea')
T_en=T
disp('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
disp(' (( javabe nahayi )) ')
disp('-----------------------------------------------------------------')
for j=1:NN
Ho_moshahedati(1,j)=Ho(j,1);
end
Ho_moshahedati=Ho_moshahedati
for j=1:NN
HS_Simulated(1,j)=feshar_haye_gereyi_shabih_sazi_shode(j,1);
end
HS_Simulated=HS_Simulated
disp('----------------------------------------------------------------')
disp('----------------------------------------------------------------')
MINIMUM_fitness_FINAL=min(f_jame_sotun_haye_ff);
Best_Fitness=MINIMUM_fitness_FINAL;
disp( [' behtarin fittnes ' ,num2str(Best_Fitness),'mibashad ke:'])
shomare_TEKRAR_e=BEST_TEKRAR;
disp([' dar shomare tekrar ',num2str(BEST_TEKRAR),'rokh dade ast.'])
Kole_tekrar_ha=t
disp('----------------------------------------------------------------')
if BEHTARIN_DEMAND_HA==0
disp('ERROR')
disp('tedad tekrar kafi nist,TEK ra ziad konid,')
disp('behtarin DEMAND ha dar tekrare aval entekhab shode ast')
disp('ke be soorate zir mibashad:')
for j=1:NN
DEMAND_haye_behine_dar_tekrare_aval(1,j)=DEMAND(j,1);
end
DEMAND_haye_behine_dar_tekrare_aval=DEMAND_haye_behine_dar_tekrare_aval
elseif BEST_TEKRAR>=2 && BEST_TEKRAR<=TEK_MIN
for j=1:NN
BEHTARIN_DEMAND_HA_nahayi(1,j)=DEMAND(j,BEST_TEKRAR);
end
BEHTARIN_DEMAND_HA_nahayi=BEHTARIN_DEMAND_HA_nahayi
else
for j=1:NN
BEHTARIN_DEMAND_HA_Nahayi(1,j)=BEHTARIN_DEMAND_HA(j,1);
end
BEHTARIN_DEMAND_HA_Nahayi=BEHTARIN_DEMAND_HA_Nahayi
end
disp('END****END****END****END****END****END****END****END****END')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% zaman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time2=clock;
hr2= time2(1,4);
min2=time2(1,5);
sec2=time2(1,6);
zaman_be_saniye=etime(time2,time1)
min00=(zaman_be_saniye)/60;
daghighe=fix(min00);
saniye=round((min00-daghighe)*60);
disp ('zaman ejra')
disp([,num2str(daghighe),'daghighe va'])
disp([,num2str(saniye),'sanea'])
disp('mibashad')
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
  