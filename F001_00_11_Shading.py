'''
From CFM-Phyto

'''

from pylab import *
from af001_energy_calculation import *
from Savetxt import *
from Figsetting import *
from sf import *
from numpy import *
from matplotlib.pyplot import *

What_is_limiting=1    #0: P-limiting  1:N-limiting

DPI = 450

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#Function beging here
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
def kkI(What_is_limiting,m,Pmax,OT,Ynphoto_chl,Cnbiosynth,Nconst_protein,Nstore_max,Cnrna_variable,Ypthylakoid_chl,Pconst_other,Qp_max,Cessential,CNprotein):   #this function calculate for the same irradiance
    I=20

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Parameters
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if What_is_limiting==0:
    #For P-limiting case
        Pin=0.002  #(mol/m3) Phosphorus concentration in the incoming medium (Healey 1985)
        Nin=0.2    #(mol/m3) Nitrate concentration in the incoming medium (Healey 1985)
        Qc=1.00*10**(-12)/12      #(molC/cell) biomass C per cell (196-18)(from Healey 1985)
    
    elif What_is_limiting==1:
    #For N-limiting case
        Pin=0.02  #(mol/m3) Phosphorus concentration in the incoming medium (Healey 1985)
        Nin=0.05    #(mol/m3) Nitrate concentration in the incoming medium (Healey 1985)
        Qc=10**(-12)/12      #(molC/cell) biomass C per cell (196-18)(from Healey 1985)

    E3=evalue()
    E=E3.E
    E=2.534 #From "Fig1 MIT.xlsx"
    
    Ddmax=1.6
    Dstep=0.001
    Dd=arange(Dstep,Ddmax+Dstep,Dstep)       #(h-1) growth rate
    U=arange(0,Ddmax/Dstep,1).astype(int)
    D=Dd/(3600*24)
    
    #==============================
    #New parameters
    #==============================
    Pchl=Pmax*(1-exp(-OT*I)) #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
    print(Pchl)
    ls=D*Qc                    #(molC s-1) Biomass synthesis rate (193-25)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Key parameters: parameterization ideas -> Kei 193-28
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #================================
    
    Molar_mass_DNA_AT_average=307.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_DNA_CG_average=307.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    Molar_mass_RNA_AT_average=316.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_RNA_CG_average=323.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    #================================
    #E coli
    #================================
    CG_Ecoli=0.506          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
    AT_Ecoli=1-CG_Ecoli     #(dimensionless) 
    
    Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of DNA unit
    Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of RNA unit
    
    RNA_DNA_mass_ratio=17.844/6.5239  #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
    
    RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli    #(mol mol-1)
    
    #================================
    #Stoichiometric parameters
    #================================
    YcyanoC_N=2                             #(molC molN) C/N molar ratio of cyanophycin
    YpgC_P=40                           #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
    
    GC=0.315                   #GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [https://www.ncbi.nlm.nih.gov/genome/13712 (accessed 04/8/2020)]
    YnucacidP_N=1/(3.5*(1-GC)+4*GC)               #(molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"
    
    YdnaC_N=3.5*(1-GC)+2.5*GC       #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    YrnaC_N=3.25*(1-GC)+2.5*GC      #(molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)

    DNAmb=1.07631                   #(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs [https://www.ncbi.nlm.nih.gov/genome/13712 (accessed 04/8/2020)]
    Avogadro=6.022*10**23           #(molecules mol-1) Avogadro constant
    
    #--------------------------
    # DNA adjusting factor
    #--------------------------
    QcProchloro = (49.2e-15)/12 #(mol cell-1) See the draft "05 K table s1.docx"
    Pdna_const=DNAmb*2*10**6/Avogadro*Qc/QcProchloro                #(molP cell-1) Constant part of DNA in phosphorus 
    Prna_const=Pdna_const*RNA_DNA_molar_ratio       #(molP cell-1) Constant part of RNA in phosphorus
    Ndna_const=Pdna_const/YnucacidP_N      #(molN cell-1) Constant part of DNA in nitrogen
    Nrna_const=Ndna_const*RNA_DNA_molar_ratio   #(molN cell-1) Constatn part of RNA in phosphorus

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Calculation
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Chl=((1+E)*ls+m)/Pchl       #(molC chl cell-1) Chlrophyll concentration (193-25) 
    Nphoto=Chl*Ynphoto_chl  #(molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth=D*Cnbiosynth             #(molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein=Nphoto+Nconst_protein+Nbiosynth    #(molN cell-1) All the proteins in N (193-26)
    Nrna_variable=Nprotein*D*Cnrna_variable        #(molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Ndna_variable=0*Dd                     #(molN cell-1) While Bremer and Dennis shows increasing trend, Parrott 1980 shows decreasing trend. 
   
    Nchl=Chl*4/55                           #(molN cell-1) Chlorophyll nitrogen (actually almost negligiable)
    Pthylakoid=Chl*Ypthylakoid_chl          #(molP cel-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
    Prna_variable=Nrna_variable*YnucacidP_N     #(molP cell-1) variable part of phosphorus in RNA (193-26)
    Pdna_variable=Ndna_variable*YnucacidP_N     #(molP cell-1) variable part of phosphorus in DNA (193-26)
    
    #=================================
    #Total calculation
    #=================================
    Qn_max=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl+Nstore_max      #(molN cell-1)  nitrogen content in the cell (193-26)                                 #(molN cell-1) total phosphorus content in the cell (193-26)                                                                             #(molP cell-1) total phosphorus content in the cell (193-26)
    Qn_min=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl            #(molN cell-1) total nitrogen in the cell without storage
    Qp_min=Pconst_other+Pthylakoid+Prna_variable+Prna_const+Pdna_variable+Pdna_const      #(molP cell-1) total phosphorus in the cell without storage
    
    #=================================
    #Vector preparation
    #=================================
    Nstore=zeros(size(Dd))
    X=zeros(size(Dd))
    Qn_test=zeros(size(Dd))
    Qp_test=copy(X)
    Qp=copy(X)
    Qn=copy(X)
    Pstore=copy(X)
    Limitation=copy(X)
    #=================================
    #Population calculation
    #=================================
    Xn_max=Nin/Qn_min
    Xp_max=Pin/Qp_min
    for i in U:
        if Xn_max[i]>Xp_max[i]:
            X[i]=Xp_max[i]
            Qp[i]=Qp_min[i]
            Qn_test[i]=Nin/X[i]
            Limitation[i]=0
            if Qn_test[i]<Qn_max[i]:
                Qn[i]=Qn_test[i]
                Nstore[i]=Qn_test[i]-Nprotein[i]-Nrna_variable[i]-Nrna_const-Ndna_variable[i]-Ndna_const-Nchl[i]  #(molN cell-1) Nitrogen storage in the cell
            else:
                Qn[i]=Qn_max[i]
                Nstore[i]=Nstore_max
        else:
            X[i]=Xn_max[i]
            Qn[i]=Qn_min[i]
            Qp_test[i]=Pin/X[i]
            Limitation[i]=1
            if Qp_test[i]<Qp_max:
                Qp[i]=Qp_test[i]
            else:
                Qp[i]=Qp_max
            Pstore[i]=Qp[i]-Pconst_other-Pthylakoid[i]-Prna_variable[i]-Prna_const-Pdna_variable[i]-Pdna_const   #(molP cell-1) Stored phosphorus in the cell
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 3 (unit adjustment)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #=======================================
    #Calculation of carbon usage (195-16)
    #=======================================
    #CNprotein=4.49   #(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
    
    #-------------------
    #C protein
    #-------------------
    Cphoto=Nphoto*CNprotein     #(molC cell-1) carbon in photosystem protein (195-16)
    Cbiosynth=Nbiosynth*CNprotein   #(molC cell-1) carbon in biosynthesis protein (195-16)
    Cconst_protein=Nconst_protein*CNprotein  #(molC cell-1) carbon in other protein assumed constant (195-16)
 
    #----------------------
    #C chlorophyll
    #----------------------
    Cchl=Chl                    #(molC cell-1) carbon in chlorophyll (195-16)
    
    #----------------------
    #C DNA RNA
    #----------------------
    Crna_const=Nrna_const*YrnaC_N       #(molC cell-1) carbon in variable part of RNA (195-16)
    Crna_variable=Nrna_variable*YrnaC_N     #(molC cell-1) carbon in variable part of RNA (195-16)
    
    Cdna_const=Ndna_const*YdnaC_N       #(molC cell-1) carbon in constant part of DNA (195-16)
    Cdna_variable=Ndna_variable*YdnaC_N     #(molC cell-1) carbon in variable part of DNA (195-16)
    
    Cnstore=Nstore*YcyanoC_N        #(molC cell-1) carbon in nitrogen storage (cyanophycin)
    CthylakoidPG=Pthylakoid*YpgC_P           #(molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes
    CthylakoidPG=0
    
    #---------------------------------------------------
    #C other: Here revised to include Nstore reduction
    #---------------------------------------------------

    #print(Qc)
    Cother_without_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-CthylakoidPG
    
    Cother_with_full_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-Cnstore-CthylakoidPG
            
    Cother=Cother_with_full_Nstore            
    
    Nstore_reduce=logical_and(Cother_without_Nstore>0, Cother_with_full_Nstore<0)
    
    Cnstore[Nstore_reduce]=Cother_without_Nstore[Nstore_reduce]
    Nstore0=copy(Nstore)
    Nstore[Nstore_reduce]=Cother_without_Nstore[Nstore_reduce]/YcyanoC_N
    Cother[Nstore_reduce]=0
    Qn[Nstore_reduce]=Qn[Nstore_reduce]+Nstore[Nstore_reduce]-Nstore0[Nstore_reduce]
    
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 1
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    BiomassC=12*X*Qc   #(mg C L-1) Biomass concentration
    NtoCplot=Qn/Qc#*14*10**6/(12*10**3)    #(mol N / mol C) biomass N to C ratio (164-20)
    
    PtoCplot=Qp/Qc#*30.97*10**6/(12*10**3)    #(mol P / mol C) biomass P to C ratio (164-20)
    NtoPplot=Qn/Qp#*14*10**6/(30.97*10**6)        #(mol N /mol P) biomass N to P ratio (164-20)
    ChltoC0=Chl/Qc         #(mol C chl mol C -1) Chlorophyll to carbon ratio
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll
    ChltoCplot=ChltoC0/12*Mchl/55     #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
    
    Dd[Cother<0]=nan

    ChltoC0=Chl/Qc         #(mol C chl mol C -1) Chlorophyll to carbon ratio
  #  ChltoCplot=ChltoC0/12/1000*Mchl/55*10**6     #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
    
    #=======================================
    #C for plot in
    #=======================================
    Numbertoarray=ones(size(Dd))            #(dimensionless) Number to array converter
    percentorratio=100       #100: percent, 1:ratio
    Cphoto_plot=Cphoto/Qc*percentorratio           
    Cbiosynth_plot=Cbiosynth/Qc*percentorratio
    Cconst_protein_plot=Cconst_protein/Qc*percentorratio*Numbertoarray
    Cchl_plot=Cchl/Qc*percentorratio
    Crna_const_plot=Crna_const/Qc*percentorratio*Numbertoarray
    Crna_variable_plot=Crna_variable/Qc*percentorratio
    Cdna_const_plot=Cdna_const/Qc*percentorratio*Numbertoarray
    Cdna_variable_plot=Cdna_variable/Qc*percentorratio
    Cother_plot=Cother/Qc*percentorratio
    Cessential_plot=Cessential/Qc*percentorratio*Numbertoarray
    Cnstore_plot=Cnstore/Qc*percentorratio
    CthylakoidPG_plot=CthylakoidPG/Qc*percentorratio
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #5.Plot   Here updated version; copy and past to the other ones and use $\mu$ for mu
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
    #==================================
    #Plot control * 1=on other=off
    #==================================
    
    a = genfromtxt('..\\Data\\Data_Felcmanova2.csv',delimiter=',').T
    
    Mu = a[13]
    Lip = a[14]
    Pro = a[15]
    Carb = a[16]
    NC = 1/a[17]/14*12
    ChlC = a[18]
    
    SDlip = a[21]
    SDpro = a[22]
    SDcarb = a[23]
    SDnc = a[24]/a[17]/a[17]/14*12
    SDChlC = a[25]

    def ebar(x,y,e):
        errorbar(x,y,(e,e),fmt='o',color='red',\
             elinewidth=1,markeredgecolor='k',ecolor='k',capthick=1,capsize=5)

    def ebar2(x,y,e,l):
        errorbar(x,y,(e,e),fmt='o',color='red',\
             elinewidth=1,markeredgecolor='k',ecolor='k',capthick=1,capsize=5,label=str(l))
    
    def ti(name):
        title(name,y=1.02)

    rcParams.update({'font.size': 15})
    
    def sup(n):
        subplot(2,3,n)
        xlabel('$\mathit{\mu}$ (d$^{-1}$)')
        ylabel('C allocation ($\%$)')
        xlim(0,0.4)
    
    #++++++++++++++++
    # For Fig 1
    #++++++++++++++++
    figure(1,figsize=(6,4.8))
 #   sup(1)
    #figure(1)
    Names = ["Oth","Pho","Bio","Sto"]
    Colors = ['#1F77B4','#FF7F0E','#2CA02C','#FFFF66']
    stackplot([],[],[],[],[],colors=Colors[::-1],labels=Names[::-1])
    stackplot(Dd,Cessential_plot+Cconst_protein_plot+Cdna_const_plot+Cdna_variable_plot,\
              Cphoto_plot+CthylakoidPG_plot+Cchl_plot,Crna_const_plot+Crna_variable_plot+Cbiosynth_plot,Cnstore_plot+Cother_plot,colors=Colors)
    #legend(loc=4)
    xlabel('$\mu$ (d$^{-1}$)')                       #copied from 73
    ylabel('C allocation ($\%$)')  #copied from 73
    ylim(top=percentorratio+1e-20)
    title('MIT 9313, Light: '+str(I)+' $\mu$mol m$^{-2}$ s$^{-1}$', y=1.05)
    xlim(0,0.4)
    DD1 = array([0,0])
    DD3 = array([100,100])
    fill_between(array([0,0.15]), DD1, DD3, where=DD3 >= DD1, facecolor='#ffffff',edgecolor = "none",alpha=0.5,zorder=10)
    
#     sup(2)
#     plot(Dd, NtoCplot, label=str(I)+'$\mu$E m$^{-2}$s$^{-1}$',color='r',zorder = -1)
#     ebar(Mu,NC,SDnc)
#     xlabel('$\mu$ (d$^{-1}$)')
#     ylabel('mol N mol C$^{-1}$')
#     ylim(bottom=0)
# 
#     sup(3)
#     plot(Dd, ChltoCplot,color='r',zorder = -1)
#     ebar(Mu,ChlC,SDChlC)
#     xlabel('$\mu$ (d$^{-1}$)')
#     ylabel('Chl:C (g g$^{-1}$)')
#     ylim(bottom=0) 
      
    sf('MIT','01',DPI)
    
    rcParams.update({'font.size': 20,
                     'lines.markersize':15,
                     'lines.markeredgewidth':1})
    rcParams['figure.figsize']=6,4.8
    
    figure(5)
    Names = ["Oth","Pho","Bio","Sto"]
    Colors = ['#1F77B4','#FF7F0E','#2CA02C','#FFFF66']
    stackplot([],[],[],[],[],colors=Colors[::-1],labels=Names[::-1])
    stackplot(Dd,Cessential_plot+Cconst_protein_plot+Cdna_const_plot+Cdna_variable_plot,\
              Cphoto_plot+CthylakoidPG_plot+Cchl_plot,Crna_const_plot+Crna_variable_plot+Cbiosynth_plot,Cnstore_plot+Cother_plot,colors=Colors)
   # legend(loc=4)
    xlabel('$\mathit{\mu}$ (d$^{-1}$)')                       #copied from 73
    ylabel('C ($\%$)')  #copied from 73
    ylim(top=percentorratio+1e-20)
  #  title('Functional allocation',y=1.02)
    xlim(0,0.4)
    DD1 = array([0,0])
    DD3 = array([100,100])
    fill_between(array([0,0.15]), DD1, DD3, where=DD3 >= DD1, facecolor='#ffffff',edgecolor = "none",alpha=0.5,zorder=10)
    
    sf('MIT','05',DPI)
    
    figure(6)
    Names = ["Oth","Pho","Bio","Sto"]
    Colors = ['#1F77B4','#FF7F0E','#2CA02C','#FFFF66']
    stackplot([],[],[],[],[],colors=Colors[::-1],labels=Names[::-1])
    legend()
    sf('MIT','06',DPI)
    
    rcParams.update({'font.size': 15,
                     'lines.markersize':15,
                     'lines.markeredgewidth':1})

    from decimal import Decimal
    def sci(Numb):
        Numb='%.2E' % Decimal(Numb)
        return Numb
    
    print("m",sci(m/Qc*86400))
    print("Pmax",sci(Pmax*86400))
    print("Apho",sci(OT))
    print("Ynphoto_chl",sci(Ynphoto_chl*CNprotein))
    print("Abio",sci(Cnbiosynth/Qc*CNprotein/86400))
    print("Cother_protein",sci(Nconst_protein/Qc*CNprotein))
    print("Arna",sci(Cnrna_variable/CNprotein/86400*YnucacidP_N))
    print("Ythylakoid_chl_P",sci(Ypthylakoid_chl))
    print("Pconst_other",sci(Pconst_other/Qc))
    print("Nstore_max",sci(Nstore_max/Qc))
    print("Cessential",sci(Cessential/Qc))
     
    #Other parameters
    print("E",sci(E))
    print("Cdna",sci(Cdna_const/Qc))
    print("Prnamin",sci(Prna_const/Qc))
    print("QpMax",sci(Qp_max/Qc))
     
    print('NrnaConst',Nrna_const)
    print('NdnaConst',Ndna_const)
    #+++++++++++++++++++++++
    # For fig 3
    #+++++++++++++++++++++++
    fc0 = 0.71
    fc1 = 1.0
    
    al = 1/100*(Cconst_protein_plot + Cphoto_plot + Cbiosynth_plot + Cessential_plot + Cother_plot + CthylakoidPG_plot)
    Cconst_protein_plot = Cconst_protein_plot/al
    Cphoto_plot = Cphoto_plot/al
    Cbiosynth_plot = Cbiosynth_plot/al
    Cessential_plot = Cessential_plot/al
    Cother_plot = Cother_plot/al
    CthylakoidPG_plot = CthylakoidPG_plot/al
    
    pro = (Cphoto_plot + Cbiosynth_plot + Cconst_protein_plot)
    carb = (Cessential_plot*fc0 + Cother_plot*fc1)
    lip = (CthylakoidPG_plot + Cessential_plot*(1-fc0) + Cother_plot*(1-fc1))
    
    figure(3,figsize=(15,7))
    sup(1)
    plot(Dd,pro,'red',label='Model',zorder = -1)
    ebar2(Mu,Pro,SDpro,'Data')
    title('Protein')
    ylim(0,70)
    legend(loc=2)
    
    sup(2)
    plot(Dd,carb,'red',zorder=-1)
    ebar(Mu,Carb,SDcarb)
    title('Carbohydrate')
    ylim(0,100)
    ylim(bottom=0)
    
    sup(3)
    plot(Dd,lip,'red',zorder = -1)
    ebar(Mu,Lip,SDlip)
    title('Lipid')
    ylim(bottom=0)
    ylim(0,30)
    
    sup(4)
    Names = ["Oth","Pho","Bio"]
    Colors = ['#9BC2E6','#FFD966',"#A9D08E"]
    stackplot([],[],[],[],colors=Colors[::-1],labels=Names[::-1])
    stackplot(Dd,Cconst_protein_plot,Cphoto_plot,Cbiosynth_plot,colors=Colors)
    legend(loc=2)
    ylim(0,70)
    
    sup(5)
    Names = ["Oth","Sto"]
    Colors = ['#C65911','#F4B084']
    stackplot([],[],[],[],colors=Colors[::-1],labels=Names[::-1])
    stackplot(Dd,Cessential_plot*fc0,Cother_plot*fc1,colors=Colors)
    legend(loc=1)
    ylim(0,100)


    sup(6)
    Names = ["Oth","Sto"]
    Colors = ['#C65911','#F4B084']
    stackplot([],[],[],[],colors=Colors[::-1],labels=Names[::-1])
    stackplot(Dd,Cessential_plot*(1-fc0),Cother_plot*(1-fc1),CthylakoidPG_plot,colors=Colors)
    legend(loc=1)
    ylim(0,30)
    
    sf('MIT','03',DPI)

    rcParams.update({'font.size': 20,
                     'lines.markersize':15,
                     'lines.markeredgewidth':1})
    rcParams.update({'xtick.major.pad': 15})
    rcParams.update({'xtick.major.pad': 15})
    rcParams.update({'figure.autolayout': True})
    rcParams['figure.figsize']=6,4.8
    rcParams.update({'figure.facecolor':'w'})
    rcParams.update({'lines.linewidth':3})   
    
    
    rcParams.update({'axes.linewidth':1.5})
    rcParams.update({'xtick.major.width':1})
    rcParams.update({'ytick.major.width':1})
    rcParams.update({'mathtext.default': 'regular' })
    
    Dd0 = copy(Dd)
    Dd1 = copy(Dd)
    Dd0[Dd0>0.15] = nan
    Dd1[Dd1<0.15] = nan
    
    figure(4)
    plot(Dd0,carb,':',color = '#0000fF')
    plot(Dd1,carb,'blue',label='Carb')
    errorbar(Mu,Carb,SDcarb,fmt='o',color='blue',markeredgecolor='k',elinewidth=2,ecolor='blue',capthick=2,capsize=10)
    plot(Dd0,pro,':',color ='red')
    plot(Dd1,pro,'red',label='Protein')
    errorbar(Mu,Pro,SDpro,fmt='o',color='red',markeredgecolor='k',elinewidth=2,ecolor='red',capthick=2,capsize=10)
    plot(Dd0,lip,':',color='green')
    plot(Dd1,lip,'green',label='Lipid')
    errorbar(Mu,Lip,SDlip,fmt='o',color='green',markeredgecolor='k',elinewidth=2,ecolor='green',capthick=2,capsize=10)
    xlabel('$\mathit{\mu}$ (d$^{-1}$)')
    ylabel('C ($\%$)')
    print(lip)
    #title('Molecular allocation',y=1.02)
    title('MIT 9313, Light = 20 $\mu$mol m$^{-2}$ s$^{-1}$',fontsize=20, y = 1.05) 
    xlim(0,0.4)
    ylim(0,100)
  #  legend(fontsize=18,loc=1)

    sf('MIT','04',DPI)
    
    figure(45)
    plot(Dd0,carb,':',color = '#0000fF')
    plot(Dd1,carb,'blue',label='Carb')
    errorbar(Mu,Carb,SDcarb,fmt='o',color='blue',markeredgecolor='k',elinewidth=2,ecolor='blue',capthick=2,capsize=10)
    plot(Dd0,pro,':',color ='red')
    plot(Dd1,pro,'red',label='Protein')
    errorbar(Mu,Pro,SDpro,fmt='o',color='red',markeredgecolor='k',elinewidth=2,ecolor='red',capthick=2,capsize=10)
    plot(Dd0,lip,':',color='green')
    plot(Dd1,lip,'green',label='Lipid')
    errorbar(Mu,Lip,SDlip,fmt='o',color='green',markeredgecolor='k',elinewidth=2,ecolor='green',capthick=2,capsize=10)
    xlabel('$\mathit{\mu}$ (d$^{-1}$)')
    ylabel('C ($\%$)')
    print(lip)
    #title('Molecular allocation',y=1.02)
    #title('MIT 9313, Light = 20 $\mu$mol m$^{-2}$, s$^{-1}$',fontsize=20, y = 1.02) 
    xlim(0,0.4)
    ylim(0,100)
    legend(fontsize=18,loc=1)

    sf('MIT','04.5',DPI)
    
    return

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAA
# Main part
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#=============================
# Parameter sets
#=============================
CNprotein = 4.49#3.82     #(molC molN) From Geider (2002) 
pf = 4.49/CNprotein   #parameter adjustment due to change in CN ratios in proteins: 4.49 is original

m=3.79146798299876E-19         #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
m=0
Pmax=0.003411  #Value from "Fig1 MIT.xlsx"

OT= 0.02439   #Value from "Fig1 MIT.xlsx"

Cnbiosynth= 4.34728279914354E-10*pf      #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
Nconst_protein= 0.0E-15*pf     #(molN cell-1) Constant protein pool in nitrogen (193-25)
Nstore_max= 2.91679384515998E-15         #(molN cell-1) Constant protein pool in nitrogen (193-25)
Cnrna_variable= 6212.59249917364/pf        #(s) Constant for Variable part of RNA (193-26)
Ypthylakoid_chl= 0.0281633095303638        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
Pconst_other= 5.44534485638617E-17               #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
Qp_max= 25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)


#parameters to tune:
Ynphoto_chl= 8.9*pf          #((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
Cessential= 1.03*4E-14          #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
#==============================

kkI(What_is_limiting,m,Pmax,OT,Ynphoto_chl,Cnbiosynth,Nconst_protein,Nstore_max,Cnrna_variable,Ypthylakoid_chl,Pconst_other,Qp_max,Cessential,CNprotein)

show()
    
    