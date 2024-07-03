def initialize_additive_DF_array(database_suffix):
    global DF_V,DF_Q,DU_V,DU_Q
    AA_type=['E','D','H','K','Y','8','9','R','C','A','F','G','L','I','M','N','P','Q','T','V','S','W']
    max_neigh=6
    DF_V=np.array([[[float('NaN') for l in range(max_neigh)] for k in range(len(AA_type))] for i in range(0,7)])
    DU_V=np.array([[[float('NaN') for l in range(max_neigh)] for k in range(len(AA_type))] for i in range(0,7)])
    for i in range(max_neigh):
        for a0 in range(0,7):
            try :
                if AA_type[a0]!='R':
                    temp_AA=AA_type[a0]
                    data_U=read_file(MODULE_PATH+'/MEAN_FIELD'+database_suffix+'/VOLUME/Prediction_parameters/U/'+temp_AA+'_'+str(i+1)+'.txt',silent=True)
                    data_F=read_file(MODULE_PATH+'/MEAN_FIELD'+database_suffix+'/VOLUME/Prediction_parameters/F/'+temp_AA+'_'+str(i+1)+'.txt',silent=True)
            except : 
                continue
            for a1 in range(len(data_F)):
                for a2 in range(len(AA_type)):
                    if data_F[a1][0].upper()==AA_type[a2] :
                        DF_V[a0,a2,i]=float(data_F[a1][1])
                        DU_V[a0,a2,i]=float(data_U[a1][1])
    
    
    DF_Q=np.array([[[float('NaN') for l in range(0,max_neigh)] for k in range(len(AA_type))] for i in range(0,7)])
    DU_Q=np.array([[[float('NaN') for l in range(0,max_neigh)] for k in range(len(AA_type))] for i in range(0,7)])
    for i in range(max_neigh):
        for a0 in range(0,7):
            try :
                if AA_type[a0]!='R':
                    temp_AA=AA_type[a0]
                    data_F=read_file(MODULE_PATH+'/MEAN_FIELD'+database_suffix+'/CHARGE/Prediction_parameters/F/'+temp_AA+'_'+str(i+1)+'.txt',silent=True)
                    data_U=read_file(MODULE_PATH+'/MEAN_FIELD'+database_suffix+'/CHARGE/Prediction_parameters/U/'+temp_AA+'_'+str(i+1)+'.txt',silent=True)
            
            except :
                continue
            for a1 in range(len(data_F)):
                for a2 in range(len(AA_type)):
                    if data_F[a1][0].upper()==AA_type[a2] :
                        DF_Q[a0,a2,i]=float(data_F[a1][1])
                        DU_Q[a0,a2,i]=float(data_U[a1][1])
def read_file(file_name,silent=False,split='') :
    data=load_file(file_name,silent=silent)
    len1=len(data)
    arr=[[] for i in range(len1)]
    for i in range(len1):
        if split!='':
            temp=data[i].split(split)
        else : 
            temp=data[i].split()
        arr[i]+=temp
    return arr

def load_file(file_name,silent=False):
    try :
        with open(file_name,'r') as f :
            data=f.readlines()
            f.close()
    except :
        if silent==False:
            print("Could not open "+file_name)
    return data

def convert_AA_1_to_3(seq,mode=1) :
    dic={   "A" : "ALA",
            "C" : "CYS",
            "D" : "ASP",
            "d" : "ASH",
            "E" : "GLU",
            "e" : "GLH",
            "F" : "PHE",
            "G" : "GLY",
            "H" : "HIP",
            "h" : "HID",
            "8" : "HEX",
            "9" : "HDX",
            "I" : "ILE",
            "K" : "LYS",
            "k" : "LYD",
            "L" : "LEU",
            "M" : "MET",
            "N" : "ASN",
            "P" : "PRO",
            "Q" : "GLN",
            "R" : "ARG",
            "S" : "SER",
            "T" : "THR",
            "V" : "VAL",
            "W" : "TRP",
            "y" : "TYR",
            "Y" : "TYO",
            "z" : "NME",
            "_" : ""}
    if mode==1:
        W=''
        for i in range(len(seq)):
            try : 
                W+=dic[seq[i]]+"\n"
            except : 
                print("Could not find "+seq[i]+", exiting.")
                exit(1)
    if mode==2:
        W=[]
        for i in range(len(seq)):
            try :
                W+=[dic[seq[i]]]
            except :
                print("Could not find "+seq[i]+", exiting.")
                exit(1)

    return W


def convert_AA_3_to_1_letter(input_arr,silent=False):
    dico = {"ALA" : "A",
            "CYS" : "C",
            "ASP" : "D",
            "ASH" : "d",
            "ASX" : "D",
            "GLU" : "E",
            "GLH" : "e",
            "GLX" : "E",
            "PHE" : "F",
            "GLY" : "G",
            "HIS" : "H",
            "HID" : "h",
            "HIE" : "h",
            "HIP" : "H",
            "HIX" : "H",
            "HEX" : "8",
            "HDX" : "9",
            "ILE" : "I",
            "LYS" : "K",
            "LYD" : "k",
            "LYX" : "K",
            "LEU" : "L",
            "MET" : "M",
            "ASN" : "N",
            "PRO" : "P",
            "GLN" : "Q", 
            "SER" : "S",
            "TYR" : "y",
            "TYO" : "Y",
            "TYX" : "Y",
            "THR" : "T",
            "VAL" : "V",
            "TRP" : "W",
            "ACE" : "z",
            "ARG" : "R", 
            "NME" : "z"}
    output_arr=[]
    for i in range(len(input_arr)):
        try :         
            output_arr+=[dico[input_arr[i]]]
        except : 
            if not silent:
                print("Did not convert "+str(input_arr[i]))
    return output_arr
          
def read_sequence(seq_3,list_res):
    
    pos_res,neg_res,base_charge,arg_res,raw_seq,seq_data_q,seq_data_id,seq_id_reduced,map,new_W=HSQ_internal_sequence_create(
        seq_3,list_res)
    titrable_residue_indexes=[]
    if titrable_residue_indexes==[]:
        list_res=np.array([i for i in range(len(seq_3))])
    else:
        list_res=np.array([titrable_residue_indexes[i] for i in range(len(titrable_residue_indexes))])

    if titrable_residue_indexes==[]:
        ind=0
        sites_num=neg_res+pos_res-arg_res
        titrable_residue_indexes=np.zeros((sites_num),dtype=int)
        for i in range(len(seq_data_q)):
            if seq_data_q[i]!=0:
                titrable_residue_indexes[ind]=i
                ind+=1
    else:
        sites_num=len(titrable_residue_indexes)
        titrable_residue_indexes=np.array(titrable_residue_indexes)

    pos_res,neg_res,base_charge,arg_res,raw_seq,seq_data_q,seq_data_id,seq_id_reduced,map,new_W=HSQ_internal_sequence_create(seq_3,list_res)

    return sites_num,titrable_residue_indexes,pos_res,neg_res,base_charge,arg_res,raw_seq,seq_data_q,seq_data_id,seq_id_reduced,map,new_W

def HSQ_internal_sequence_create(seq,list_res,no_write=False):
    
    data=seq
    new_W=''
    pos_res=0
    neg_res=0
    base_charge=0 
    arg_res=0

    seq_data_q=[]
    seq_data_id=[]
    raw_seq=[]
    raw_seq_2=[]
    map=[]

    for i in range(len(data)) :

        line=data[i]
        do_it=find_ind(i,list_res)
        raw_seq_2+=[line]
        
        if (line=="GLU" or line=="GLH"  or line=="GLX") and do_it==True :
            new_W+='GLX\n'
            raw_seq+=['GLX']
            neg_res+=1
            seq_data_id+=[0]
            seq_data_q+=[-1]
        elif (line=="ASP" or line=="ASH"  or line=="ASX") and do_it==True :
            new_W+='ASX\n'
            raw_seq+=['ASX']
            neg_res+=1
            seq_data_id+=[2]
            seq_data_q+=[-1]
        elif (line=="LYD" or line=="LYS" or line=="LYX") and do_it==True :
            new_W+='LYX\n'
            raw_seq+=['LYX']
            pos_res+=1
            seq_data_q+=[1]
            seq_data_id+=[1]
        elif line=="ARG" :
            new_W+='ARG\n'
            raw_seq+=['ARG']
            pos_res+=1
            arg_res+=1
            seq_data_q+=[0] 
            seq_data_id+=[float('NaN')]
        elif (line=="TYX" or line=="TYR" or line=="TYO") and do_it==True :
            new_W+='TYX\n'
            raw_seq+=['TYX']
            neg_res+=1
            seq_data_q+=[-1]
            seq_data_id+=[3]
        elif (line=="SXP" or line=="S1P" or line=="S2P") and do_it==True :
            new_W+='SXP\n'
            raw_seq+=['SXP']
            neg_res+=1
            base_charge+=-1
            seq_data_q+=[-1]
            seq_data_id+=[4]
        elif (line=="TXP" or line=="T1P" or line=="T2P") and do_it==True :
            new_W+='TXP\n'
            raw_seq+=['TXP']
            neg_res+=1
            base_charge+=-1
            seq_data_q+=[-1]
            seq_data_id+=[5]
        elif (line=="YXP" or line=="Y1P" or line=="Y2P") and do_it==True :
            new_W+='YXP\n'
            raw_seq+=['YXP']
            neg_res+=1
            base_charge+=-1
            seq_data_q+=[-1]
            seq_data_id+=[6]
        elif (line=="HIE" or line=="HID" or line=="HIS" or line=="HEX" or line=="HDX" or line=="HIP" or line=="HIX") and do_it==True :
            new_W+='HIX\n'
            raw_seq+=['HIX']
            pos_res+=1
            seq_data_q+=[1]
            seq_data_id+=[9]
        elif do_it==False and (line=="HIS" or line=="LYS"):
   
            new_W+=line+'\n'
            seq_data_q+=[0]
            seq_data_id+=[float('NaN')]
            raw_seq+=[line] 
            base_charge+=1
        elif do_it==False and (line=="ASP" or line=="GLU" or line=="TYO"):
            new_W+=line+'\n'
            seq_data_q+=[0]
            seq_data_id+=[float('NaN')]
            base_charge+=-1
            raw_seq+=[line]
        elif line!="END" :
            if line!="NA+" and line!="CL-" and line!="CLX" and line!="NAX" :
                raw_seq+=[line]
            new_W+=line+'\n'
            seq_data_q+=[0]
            seq_data_id+=[float('NaN')] 
        if seq_data_q[-1]!=0 and line!="END":
            map+=[i]  

    if no_write==False :
        write_file('seq_hsq.in',new_W)
    seq_id_reduced=[]

    for i in range(len(seq_data_id)):
        if math.isnan(seq_data_id[i])==False :
            seq_id_reduced+=[seq_data_id[i]]
    
    return pos_res,neg_res,base_charge,arg_res,raw_seq,seq_data_q,seq_data_id,seq_id_reduced,map,new_W


def find_ind(ind,arr) :
    for i in range(len(arr)):
        if arr[i]==ind :
            return True
    return False


def get_ref_pkas(mode,FF='OPLS',DF_only=True,version=-1,silent=True,suffix='') :
    if version==-1 :
        try :
            version=read_file('/project/fava/work/martinfossat/pKa_Calc/ALL_MDCP/'+str(FF)+'/ACE-GXG-NME/Latest.inf')[0][1]
        except :
            print("Could not read info on the last model compound free energy")
            version=str(-1)

    else :
        version=str(version)
    seq_ref_arr=['GLX','LYX','ASX','TYX','SXP','TXP','YXP','HDX','HEX','HIX']
    ref_notation=[['E','e'],['K','k'],['D','d'],['Y','y'],['U','u'],['X','x'],['Z','z'],['H','9'],['H','8'],['H','h']]
    factor=[1,-1,1,1,1,1,1,-1,-1,-1]
     
    ref=[]
    base='/project/fava/work/martinfossat/pKa_Calc/ALL_MDCP/'+str(FF)+'/ACE-GXG-NME/298K_50mMNaCl_'+version+suffix+'/'
    for i in range(len(seq_ref_arr)):
        try : 
            ref+=[float(read_key_file(base+str(seq_ref_arr[i])+'/0.05/_Analysis/Parameters.txt',['offset'],silent=silent)[0])] 
        except :
            if not silent : 
                print("Could not get the offset from the model compounds for "+seq_ref_arr[i])
            ref+=[0.]
    if mode==1 :     
        pKa_ref_arr=[4.34,10.34,3.86,9.76,5.96,6.3,5.96,7.15,6.55,6.45] 
   
        folders=[base+'GLX/0.05/_Analysis/Output.txt',\
        base+'LYX/0.05/_Analysis/Output.txt',\
        base+'ASX/0.05/_Analysis/Output.txt',
        base+'TYX/0.05/_Analysis/Output.txt',
        base+'SXP/0.05/_Analysis/Output.txt',
        base+'TXP/0.05/_Analysis/Output.txt',
        base+'YXP/0.05/_Analysis/Output.txt',
        base+'HDX/0.05/_Analysis/Output.txt',
        base+'HEX/0.05/_Analysis/Output.txt']
        
        DF_ref=[]
        DF_ref_err=[] 
        DS_ref=[]
        DS_ref_err=[]
        DU_ref=[]
        DU_ref_err=[]
        for i in range(len(folders)):
            try : 
                data=read_file(folders[i],silent=True)        
                DF_ref+=[float(data[-2][-1])]  
                DF_ref_err+=[float(data[-1][-1])]
                DS_ref+=[float(data[-6][-1])]  
                DS_ref_err+=[float(data[-5][-1])]
                DU_ref+=[float(data[-4][-1])]  
                DU_ref_err+=[float(data[-3][-1])]
            except : 
                DF_ref+=[float('NaN')]
                DF_ref_err+=[float('NaN')]
                DS_ref+=[float('NaN')]
                DS_ref_err+=[float('NaN')]
                DU_ref+=[float('NaN')]
                DU_ref_err+=[float('NaN')]
                if not silent : 
                    print("!!!!!!!!!!!!!!!!!!!!!!!!!!! No model compound pka value !!!!!!!!!!!!!!!!!!!!!!!!")
                    print("File is not found or is invalid : ",folders[i])

    elif mode==2 :
        print("! You should proably use mode 1 instead of mode 2 !")
        pKa_ref_arr=[4.25,10.28,3.65]
        base='/work/mfossat/pKa_Calc/ALL_MDCP/'
        folders=[base+'GLX/0.05/_Analysis/Output.txt',\
        base+'LYX/0.05/_Analysis/Output.txt',\
        base+'ASX/0.05/_Analysis/Output.txt']

        DF_ref=[]
        DF_ref_err=[]
        for i in range(len(folders)):
            data=read_file(folders[i])
            DF_ref+=[float(data[-2][-1])]
            DF_ref_err+=[float(data[-1][-1])]
    else :
        print("Only modes 1 and 2 are available")
    
    if DF_only==True:
        out=[seq_ref_arr,ref_notation,DF_ref,DF_ref_err,folders,factor,pKa_ref_arr,ref]
    else :
        out=[seq_ref_arr,ref_notation,DF_ref,DF_ref_err,folders,factor,pKa_ref_arr,ref,DS_ref,DS_ref_err,DU_ref,DU_ref_err]    
    return out

def read_key_file(key_file,pattern,silent=True):
    found_it=[0 for i in range(len(pattern))]
    data=load_file(key_file,silent=True)
    out=[]
    if type(pattern)==type([]) :
        W=''
        val=''
        for i in range(len(data)) :
            if len(data[i].split())>0:
                found=False
              
                for g in range(len(pattern)) :
                    if data[i].split()[0]==str(pattern[g]) :
                        found=True
                        found_it[g]=1
                        break

                if found==True:
                    if silent==False : 
                        print(pattern[g]," is ",data[i].split()[1])
                    out+=[data[i].split()[1]]                    
    return out             
                    
def get_base_contexts(map,neigh,seq_1,seq_data_q):
    
    contexts=[]
    states=[]

    
    for i in range(len(map)):
        tmp=''
        tmp2=[]
        for j in range(-neigh,neigh+1):
            
            if map[i]+j>=1 and map[i]+j<len(seq_1)-1:
                tmp+=seq_1[map[i]+j]
                tmp2+=[seq_data_q[map[i]+j]]
            
            else :
                tmp+='_'
                tmp2+=[float('nan')]
        contexts+=[tmp]
        states+=[tmp2]
    return contexts,states

def get_all_contexts(contexts,states,neigh,refs,unshifted,reverse,penta,additive,T,base_rep,reduced):
    import itertools
    missing_pat=False
    p=0
    pat_seq=[[] for i in range(len(contexts))]
    pat_ste=[[] for i in range(len(contexts))]
    pat_E=[[] for i in range(len(contexts))]

    
    
    
    for i in range(len(contexts)):
        count=0
        new_seq=''
        titratable_pre=[]
        for j in range(neigh):
            if contexts[i][j]=='E' or contexts[i][j]=='K' or contexts[i][j]=='H' or contexts[i][j]=='D' or contexts[i][j].upper()=='Y':
                titratable_pre+=[1]
            else :
                titratable_pre+=[0]
        titratable_post=[]
        for j in range(neigh+1,neigh*2+1):
            if contexts[i][j]=='E' or contexts[i][j]=='K' or contexts[i][j]=='H' or contexts[i][j]=='D' or contexts[i][j].upper()=='Y':
                titratable_post+=[1]
            else :
                titratable_post+=[0]
        for j in range(neigh):
            if states[i][j]==-1 or states[i][j]==1 :
                count+=1
        
        for curr in itertools.product([1,0],repeat=count):
            tmp=''
            new_seq=''
            offset=0
            txt_ste=''
            h=0

            for k in range(len(titratable_pre)):
                if titratable_pre[k]==1 :
                    if contexts[i][k]=='E' or contexts[i][k]=='D' or contexts[i][k].upper()=='Y':
                        if curr[h]==0:
                            tmp+=contexts[i][k].upper()
                        elif curr[h]==1:
                            tmp+=contexts[i][k].lower()
                    elif contexts[i][k]=='K' or contexts[i][k]=='H' :
                        if curr[h]==1:
                            tmp+=contexts[i][k].upper()
                        elif curr[h]==0:
                            tmp+=contexts[i][k].lower()
                    txt_ste+=str(curr[h])
                    h+=1

                else :
                    txt_ste+='1'
                    tmp+=contexts[i][k]

            new_seq+=tmp
            
            new_seq+=contexts[i][neigh].upper()

            for k in range(neigh+1,neigh*2+1):
                
                if contexts[i][k]=='E' or contexts[i][k]=='D' or contexts[i][k].upper()=='Y':
                    new_seq+=contexts[i][k].lower()
                elif contexts[i][k]=='K' or contexts[i][k]=='H' :
                    new_seq+=contexts[i][k].upper()
                else :
                    new_seq+=contexts[i][k]

            pat_ste[i]+=[txt_ste]
            pat_seq[i]+=[new_seq]

            
            if unshifted==True :
                for j in range(len(refs[1])):
                    if refs[1][j][1].upper()==new_seq[neigh].upper() :
                        break

                F_mdcp_expt=R*T*np.log(10**(refs[6][j]))
                if reverse==False:
                    pat_E[i]+=[-(F_mdcp_expt)]
                else :
                    pat_E[i]+=[F_mdcp_expt]
                continue
            else :
                
                if new_seq[neigh]=='H':
                    He_seq=''
                    Hd_seq=''
                    for l in range(2*neigh+1):
                        if l==neigh:
                            He_seq+='8'
                            Hd_seq+='9'
                        else :
                            He_seq+=new_seq[l]
                            Hd_seq+=new_seq[l]
                    out_e=get_pattern_energy_from_database(He_seq,refs,reverse,penta,additive,T,base_rep,neigh,silent=True,reduced=reduced)
                    out_d=get_pattern_energy_from_database(Hd_seq,refs,reverse,penta,additive,T,base_rep,neigh,silent=True,reduced=reduced)

                    if out_e==None or out_d==None :
                        missing_pat=True
                        continue
                    else :
                        F_pat_e,F_pat_err,F_mdcp_e,S_pat,S_mdcp,U_pat,U_mdcp,offset_pat_e,mdcp_offset_e,U_mdcp_expt,S_mdcp_expt,sign,F_mdcp_expt_e,T_pat=out_e
                        F_pat_d,F_pat_err,F_mdcp_d,S_pat,S_mdcp,U_pat,U_mdcp,offset_pat_d,mdcp_offset_d,U_mdcp_expt,S_mdcp_expt,sign,F_mdcp_expt_d,T_pat=out_d

                else :
                    out=get_pattern_energy_from_database(new_seq,refs,reverse,penta,additive,T,base_rep,neigh,silent=True,reduced=reduced)

                    if out==None :
                        missing_pat=True
                        continue
                    else :
                        F_pat,F_pat_err,F_mdcp,S_pat,S_mdcp,U_pat,U_mdcp,offset_pat,mdcp_offset,U_mdcp_expt,S_mdcp_expt,sign,F_mdcp_expt,T_pat=out

            if reverse==False:
                sign2=-1
            else :
                sign2=1
            if new_seq[neigh]=='H':
                
                pat_E[i]+=[get_G_sum([sign2*(F_mdcp_expt_e-sign*((F_pat_e-offset_pat_e)-(F_mdcp_e-mdcp_offset_e))),sign2*(F_mdcp_expt_d-sign*((F_pat_d-offset_pat_d)-(F_mdcp_e-mdcp_offset_e)))],T)]
            else :
                pat_E[i]+=[sign2*(F_mdcp_expt-sign*((F_pat-offset_pat)-(F_mdcp-mdcp_offset)))]
                

        pat_ste[i]=np.array(pat_ste[i])

    if missing_pat==False:
        for i in range(len(pat_E)):
            if math.isnan(np.sum(pat_E[i])):
                missing_pat=True
                break

    if (missing_pat):
        print("Some pattern free energy are missing from the database.")
        quit()
    return pat_E,pat_ste,pat_seq

def get_pattern_energy_from_database(pattern,refs,reverse,penta,additive,T,data_dir,neigh,silent=True,reduced=True,unshifted=True,ent_corr=False):
    import os
    found=False
    for j in range(len(refs[1])):
        if refs[1][j][1].upper()==pattern[neigh].upper()  :
            found=True
            break

    if found==False :
        print(pattern[neigh])
        input("This should not happen")

    sign,sign2=check_sign(pattern[neigh],refs,reverse)
    
    F_mdcp_expt=R*T*np.log(10**(refs[6][j]))
    good=0
    tot=0
    
    if penta==True :

        
        
        mdcp_offset=refs[7][j] 


        
        F_mdcp=refs[2][j]
        F_mdcp_err=refs[3][j]
        S_mdcp=refs[8][j]
        U_mdcp=refs[10][j]
        
        S_mdcp_expt=(S_mdcp/F_mdcp)*F_mdcp_expt
        U_mdcp_expt=(U_mdcp/F_mdcp)*F_mdcp_expt
    
        temp=convert_pattern_to_simplified(pattern,neigh)
        if additive==False :
            try :  
                dirs=os.listdir(data_dir+'/'+temp)
                if reduced==False : 
                    try : 
                        dirs_dum=os.listdir(data_dir+'/'+temp+'/'+pattern)
                    except : 
                        print("Did not find exact peptide ", pattern)
                        not_found=True  
                        setup_motif_database_sequence(temp,pattern,data_dir)
                        return None 
            except : 
                print("Did not find peptide ", temp)
                not_found=True  
                setup_motif_database_sequence(temp,pattern,data_dir)
                return None 
        good=0
         
        
        temp_pat=''
        for l in range(len(pattern)):
            if l==neigh :
                temp_pat+=pattern[l].upper()
            else :
                temp_pat+=pattern[l]

        for l in range(len(dirs)): 
            if dirs[l]==temp_pat : 
                dups=os.listdir(data_dir+'/'+temp+'/'+dirs[l])
                
                F_pat_err=0
                F_pat=0
                S_pat_err=0
                S_pat=0
                U_pat=0
                U_pat_err=0
                temp_vals=[]
                tot=0
                for n in range(len(dups)):
                    if os.path.isdir(data_dir+'/'+temp+'/'+dirs[l]+'/'+dups[n]):
                        try :
                            tmp=read_file(data_dir+'/'+temp+'/'+dirs[l]+'/'+dups[n]+'/_Analysis/Output.txt',silent=True)
                            
                            F_pat_err+=float(tmp[-1][-1])**2
                            F_pat+=float(tmp[-2][-1])

                            S_pat_err+=float(tmp[-5][-1])**2
                            S_pat+=float(tmp[-6][-1])
                            
                            U_pat_err+=float(tmp[-3][-1])**2
                            U_pat+=float(tmp[-4][-1])
                            
                            good+=1
                            tot+=1
                            
                            
                            
                            temp_vals+=[float(tmp[-2][-1])]

                        except :
                            continue

                if tot!=0 :
                    if tot==1 :
                        print('!!!!!!!!Warning!!!!!!!')
                        print(temp+'/'+dirs[l]+' does not have any duplicates. This could result in undetected bad free energy being used')
                    for k in range(len(temp_vals)): 
                        if abs(np.mean(temp_vals)-temp_vals[k])>0.05:
                            print(temp+'/'+dirs[l]+' Has at least one outlier') 
                     
                        F_pat=F_pat/tot
                        F_pat_err=np.sqrt(F_pat_err)
                
                        S_pat=S_pat/tot
                        S_pat_err=np.sqrt(S_pat_err)
                
                        U_pat=U_pat/tot
                        U_pat_err=np.sqrt(U_pat_err)
                    break

        
        if good==0 and reduced:
            tot=0
            F_pat_err=0
            F_pat=0
            S_pat_err=0
            S_pat=0
            U_pat_err=0
            U_pat=0
            for l in range(len(dirs)):
                dups=os.listdir(data_dir+'/'+temp+'/'+dirs[l])
                for n in range(len(dups)):
                    if os.path.isdir(data_dir+'/'+temp+'/'+dirs[l]+'/'+dups[n]): 
                        try :
                            tmp=read_file(data_dir+'/'+temp+'/'+dirs[l]+'/'+dups[n]+'/_Analysis/Output.txt',silent=True)
                      
                            F_pat_err+=float(tmp[-1][-1])**2
                            F_pat+=float(tmp[-2][-1])
        
                            S_pat_err+=float(tmp[-5][-1])**2
                            S_pat+=float(tmp[-6][-1])
        
                            U_pat_err+=float(tmp[-3][-1])**2
                            U_pat+=float(tmp[-4][-1])
        
                            tot+=1
                            good+=1
                            
                        except :     
                            continue
        
            if good==0:
                print("Could not find free energy files for ",temp)
                
                return None
             
            else : 
                F_pat_err=np.sqrt(F_pat_err)
                F_pat=F_pat/tot       
                S_pat=S_pat/tot
                S_pat_err=np.sqrt(S_pat_err)
                U_pat=U_pat/tot
                U_pat_err=np.sqrt(U_pat_err)

        elif good==0 and not reduced :
            print("For real ?")
            print(reduced)
            print("No free energy file for "+temp+"/"+pattern)
            return None
        try : 
            
            tmp=read_key_file(data_dir+'/'+temp+'/'+dirs[0]+'/1/run.key',['FMCSC_PKA_REP_POT','FMCSC_TEMP'],silent=silent)
            offset_pat=float(tmp[0])
            T_pat=float(tmp[1])
        except : 
            T_pat=T
            offset_pat=0.
            print(data_dir+'/'+temp+'/'+dirs[0]+'/1/run.key is missing')
            input("This is a problem. This file is necessary to know the offset ")

    elif additive:
        T_pat=T
        F_pat=get_pattern_additive_F(pattern,T,neigh,ent_corr=ent_corr) 
        F_pat_err=0.
        F_mdcp=0.
        S_pat=0.
        S_mdcp=0.
        U_pat=F_pat
        U_mdcp=0.
        offset_pat=0.
        mdcp_offset=0.
        U_mdcp_expt=0.
        S_mdcp_expt=0.

    else : 
        
        T_pat=T
        F_pat=0.
        F_pat_err=0.
        F_mdcp=0.
        S_pat=0.
        S_mdcp=0.
        U_pat=0.
        U_mdcp=0.
        offset_pat=0.
        mdcp_offset=0.
        U_mdcp_expt=0.
        S_mdcp_expt=0.

    if F_pat_err>0.1:
        print("pentpeptide "+pattern+" has a large error ("+str(F_pat_err)+")")
    return F_pat,F_pat_err,F_mdcp,S_pat,S_mdcp,U_pat,U_mdcp,offset_pat,mdcp_offset,U_mdcp_expt,S_mdcp_expt,sign,F_mdcp_expt,T_pat

def get_G_sum(Garr,T):
    
    Garr=np.array(Garr)
    out=-R*T*np.log(np.sum(np.exp(-(Garr)/(R*T)),axis=0))
    return out 

def check_sign(AA,refs,reverse):
    for j in range(len(refs[1])):
        if refs[1][j][1].upper()==AA.upper()  :
            break

    if j==len(refs[2]):
        return 0,0
    if reverse==True :
        sign2=-refs[5][j]
    else :
        sign2=refs[5][j]
    sign=refs[5][j]

    return sign,sign2


def convert_pattern_to_simplified(pattern,neigh,central_only=True,ion_simplify=True,ion_sides_only=True,separate=[]):
    ion='EDY-HK+89'
    if ion_simplify:
        ion=ion+'R'
        ion_simp='----++++++'

    if central_only==True :
        
        found=False
        for j in range(len(ion)):
            if ion[j]==pattern[neigh]:
                found=True
                break
        if found==False :
            return
    pattern_out=''
    for i in range(len(pattern)):
        found_one=False
        for j in range(len(ion)):
            if ion[j]==pattern[i]:
                found_one=True
                if ion_simplify==True and (i<neigh-1 or i>neigh+1):
                    if ion[j].islower():
                        pattern_out+='0'
                    else : 
                        pattern_out+=ion_simp[j]  
                else : 
                    pattern_out+=pattern[i]
                break   

        if found_one==False:
            if ion_simplify and (i<neigh-1 or i>neigh+1) and pattern[i]!='_':
                if separate!=[]: 
                    found_sep=False
                    for m in range(len(separate)):
                        if separate[m]==pattern[i]:
                            found_sep=True
                            break
                    if found_sep: 
                        pattern_out+=pattern[i]
                    else :
                        pattern_out+='0'
                else :
                    pattern_out+='0'
            elif pattern[i]=='_' :
                pattern_out+='_'
            else :
                pattern_out+=pattern[i] 

    return pattern_out
    
def setup_motif_database_sequence(reduced,pattern,base_dir):

    check_and_create_rep(base_dir)


    seq_3=convert_AA_1_to_3(pattern,mode=2)
    neigh=int(len(pattern)/2)   

    if pattern[neigh]=='E':
        seq_3[neigh]='GLX'
    elif pattern[neigh]=='D':
        seq_3[neigh]='ASX'
    elif pattern[neigh]=='K':
        seq_3[neigh]='LYX'
    elif pattern[neigh]=='H':
        seq_3[neigh]='HEX'
    elif pattern[neigh]=='Y':
        seq_3[neigh]='TYX'
    W='ACE\n'
    for o in range(len(seq_3)):
        if len(seq_3[o])>0:
            W+=seq_3[o]+'\n'
    W+='NME\nEND'

    print("Creating pentapeptide database folder to be submited :")    
    print(base_dir)
    check_and_create_rep(base_dir+'/'+reduced+'/'+pattern) 
    write_file(base_dir+'/'+reduced+'/'+pattern+'/seq.in',W)

def get_pattern_additive_F(pattern,T,neigh,ent_corr=False): 
    global refs
    global AA_type
    global groups
    global R
    global DF_V,DF_QDU_V,DU_Q 
    
    signs=refs[5] 
    F_pat=0.

    for a0 in range(len(refs[0])):
        if pattern[neigh]==AA_type[a0]: 
            break

    
    ions='EDH89KY'
    
    seq_F=[]
    seq_F_q=[]

    if ent_corr :
        DX_V=np.copy(DU_V)
        DX_Q=np.copy(DU_Q)
    else : 
        DX_V=np.copy(DF_V)
        DX_Q=np.copy(DF_Q)

    
    for a1 in range(len(AA_type)):
        if 'G'==AA_type[a1]:
            break

    F_pat+=-DX_V[a0,a1,0]*2   

    for i in range(len(pattern)): 
        for a1 in range(len(AA_type)):
            if AA_type[a1]==pattern[i].upper() and pattern[i]!='_' and i!=neigh :
                F_pat+=DX_V[a0,a1,abs(i-neigh)-1] 
                seq_F+=[DX_V[a0,a1,abs(i-neigh)-1]]
                
                if a1<=5 and  pattern[i].isupper():
                    F_pat+=DX_Q[a0,a1,abs(i-neigh)-1]     
                    seq_F_q+=[DX_Q[a0,a1,abs(i-neigh)-1]]
     
    if math.isnan(F_pat):
        print(pattern," has a problem")
        input(F_pat) 
    return F_pat

def main_prediction(T,neigh,contexts,pat_E,pat_ste,pat_seq,map,prun_during,per_res,base_E,test_time,max_diff,silent=False):
    import itertools
    blocks=[]
    tmp=[]
    for i in range(len(contexts)):
        tmp+=[i]
        if i%2==1:
            blocks+=[np.array(tmp)]
            tmp=[]

    
    dim=(len(contexts)+1,)
    for i in range(neigh):
        dim=dim+(2,)
    if not per_res:
        lvl_disc_E=np.zeros(dim,dtype=np.double)
        lvl_disc_E[:]=np.double('+inf')
    
    
    
    
    
    
    dim=(len(contexts)+1,len(contexts))
    for i in range(neigh):
        dim=dim+(2,)

    if per_res:
        lvl_context0=np.zeros(dim,dtype=np.double)
        lvl_context0[:]=float('+inf')
        lvl_context1=np.zeros(dim,dtype=np.double)
        lvl_context1[:]=float('+inf')
    else :
        lvl_ste=[np.array(['']) for i in range(len(contexts)+1)]
        
        lvl_E=[np.array([base_E],dtype=np.double) for i in range(len(contexts)+1)]
    tot=0
    
    for i in range(len(pat_E)):
        if not silent :
            print(i+1," out of ",len(pat_E))
        temp_E0=[]
        for t in range(i+1):
            if not per_res:
                temp=np.zeros((len(lvl_ste[t])),dtype=np.double)
                temp[:]=float('nan')
                for h in range(len(lvl_ste[t])):
                    break_it=False
                    
                    for l in range(len(pat_ste[i])):
                        ind=[]
                        offset=0
                        for o in range(1,neigh+1):
                            if i-o>=0:
                                while map[i]-o-offset!=map[i-o] and neigh-o-offset>=0:
                                    offset+=1
                                if neigh-o-offset>=0:
                                    ind=[neigh-o-offset]+ind

                        min_ind=max(0,i-len(ind))
                        count_good=0
                        for o in range(len(ind)) :
                            if lvl_ste[t][h][min_ind+o]==str(pat_ste[i][l][ind[o]]):
                                count_good+=1

                        if count_good==len(ind):
                            break_it=True
                            break

                    if break_it==False:
                        print(pat_seq[i])
                        print(pat_ste[i])
                        print(lvl_ste[t])
                        input("Problem")

                    temp[h]=lvl_E[t][h]+pat_E[i][l]
                temp_E0+=[np.array(temp)]
        if not per_res :
            temp_ste1=[np.char.add(lvl_ste[t],np.array(['1'])) for t in range(i+1)]
            temp_ste0=[np.char.add(lvl_ste[t],np.array(['0'])) for t in range(i+1)]

        dim=(2,)
        for l in range(neigh-1):
            dim=dim+(2,)

        temp_pat_E=np.zeros(dim,dtype=np.double)
        temp_pat_E[:]=float('+inf')

        
        for l in range(len(pat_ste[i])):
            all_ind=[()]
            ind=[]
            offset=0
            for o in range(1,neigh+1):
                if i-o>=0:
                    while neigh-o-offset>=0:
                        if map[i]-o-offset==map[i-o]:
                            ind=[neigh-o-offset]+ind
                            for h in range(len(all_ind)):
                                all_ind[h]=(int(pat_ste[i][l][neigh-o-offset]),)+all_ind[h]
                            break
                        else:
                            offset+=1
            true_ind=[]
            for curr in  itertools.product([0,1],repeat=neigh-len(ind)):
                ind_temp=()
                for h in range(len(ind)):
                    ind_temp+=(int(pat_ste[i][l][ind[h]]),)
                true_ind+=[curr+ind_temp]
            for h in range(len(true_ind)):
                temp_pat_E[true_ind[h]]=pat_E[i][l]

        if prun_during and not per_res :
            disc_temp_E=np.zeros(np.shape(lvl_disc_E),dtype=np.double)
            disc_temp_E[:,:]=float('+inf')

        if per_res:
            lvl_temp_context0= np.zeros(np.shape(lvl_context0),dtype=np.double)
            lvl_temp_context0[:,:]=float('+inf')

            lvl_temp_context1= np.zeros(np.shape(lvl_context1),dtype=np.double)
            lvl_temp_context1[:,:]=float('+inf')

        
        for t in range(i+1):
            
            if per_res==True:
                if t==0:
                    curr1=()
                    for k in range(neigh):
                        curr1+=(1,)
                    curr0=curr1[:-1]+(0,)
                    
                    
                    lvl_context0[(t+1,i,)+curr0]=base_E+temp_pat_E[curr1]
                    lvl_context1[(t,i,)+curr1]=base_E
                for k in range(i):
                    
                    base1=(t,k)
                    base0=(t+1,k)
                    for curr in itertools.product([0,1],repeat=neigh):
                        cont=base1+curr
                        pat_cont=curr

                        new_cont0=base0+curr[1:]+(0,)
                        new_cont1=base1+curr[1:]+(1,)

                        lvl_temp_context0[new_cont0]=get_G_sum([lvl_temp_context0[new_cont0],lvl_context0[cont]+temp_pat_E[pat_cont]],T)
                        lvl_temp_context1[new_cont0]=get_G_sum([lvl_temp_context1[new_cont0],lvl_context1[cont]+temp_pat_E[pat_cont]],T)
                        lvl_temp_context0[new_cont1]=get_G_sum([lvl_temp_context0[new_cont1],lvl_context0[cont]],T)
                        lvl_temp_context1[new_cont1]=get_G_sum([lvl_temp_context1[new_cont1],lvl_context1[cont]],T)

                k=i
                
                if t==0 :
                    k1=k
                    k2=k
                if t!=0 :
                    k1=k-1
                    k2=k
                base1k1=(t,k1)
                base1k2=(t,k2)
                base0k1=(t+1,k1)
                base0k2=(t+1,k2)
                
                for curr in itertools.product([0,1],repeat=neigh):
                    cont=base1k1+curr
                    pat_cont=curr
                    new_cont0=base0k2+curr[1:]+(0,)
                    new_cont1=base1k2+curr[1:]+(1,)
                    
                    lvl_temp_context0[new_cont0]=get_G_sum([lvl_temp_context0[new_cont0],lvl_context0[cont]+temp_pat_E[pat_cont],lvl_context1[cont]+temp_pat_E[pat_cont]],T)
                    lvl_temp_context1[new_cont1]=get_G_sum([lvl_temp_context1[new_cont1],lvl_context1[cont],lvl_context0[cont]],T)

            if prun_during and not per_res:
                base0=(t+1,)
                base1=(t,)
                for curr in itertools.product([0,1],repeat=neigh):
                    cont=base1+curr
                    pat_cont=curr
                    new_cont0=base0+curr[1:]+(0,)
                    new_cont1=base1+curr[1:]+(1,)

                    disc_temp_E[new_cont0]=get_G_sum([disc_temp_E[new_cont0],lvl_disc_E[cont]+temp_pat_E[pat_cont]],T)
                    disc_temp_E[new_cont1]=get_G_sum([disc_temp_E[new_cont1],lvl_disc_E[cont]],T) 

        if per_res :
            lvl_context0=lvl_temp_context0
            lvl_context1=lvl_temp_context1

        if prun_during and not per_res:
            lvl_disc_E=disc_temp_E
        if per_res:
            continue
        
        
        for t in range(i+1):
            lvl_ste[t]=temp_ste1[t]

        
        for t in range(i+1):
            if len(lvl_ste[t+1][0])!=0:
                lvl_ste[t+1]=np.append(lvl_ste[t+1],temp_ste0[t])
                lvl_E[t+1]=np.append(lvl_E[t+1],temp_E0[t])
            else :
                lvl_ste[t+1]=temp_ste0[t]
                lvl_E[t+1]=temp_E0[t]

        for t in range(i+1):
            
            if prun_during and not (math.isinf(max_diff) and max_diff>0):
                if len(lvl_ste[t])<2 :
                    continue
                
                if lvl_ste[t][0]=='':
                    break
                
                W=np.exp(-(lvl_E[t])/(R*T))
                ind=np.argsort(lvl_E[t])
                
                if not (math.isinf(max_diff)) and (max_diff<0):
                    Z_up=np.cumsum(W[ind[::-1]])[::-1]
                    Z_down=np.cumsum(W[ind[::1]])

                    excess=np.sum(np.exp(-lvl_disc_E[t]/(R*T)))
                    G_down=-R*T*np.log(Z_down)
                    G_up=-R*T*np.log(Z_up)

                    G_up_all =-R*T*np.log(np.exp(-G_up /(R*T)) +excess)
                    G_down_all =-R*T*np.log(np.exp(-G_down /(R*T)) +excess)
                    temp=np.where(G_up_all>G_down_all+max_diff)

                    if len(temp[0])<1:  
                        continue
                    elif len(temp[0])>=len(ind):  
                        temp=[[temp[0][-1]]]
                else :
                    temp=[[1]]

                ste_temp=lvl_ste[t][ind[temp[0][0]:]]
                lvl_ste[t]=lvl_ste[t][ind[:temp[0][0]]]
                lvl_E[t]=lvl_E[t][ind[:temp[0][0]]]
                W_temp=W[ind[temp[0][0]:]]

                
                
                
                for curr in itertools.product([1,0],repeat=neigh):
                    temp_char=''
                    for c in range(-min(0,len(ste_temp[0])-neigh),len(curr)):
                        temp_char+=str(curr[c])
                    bool_temp=np.char.endswith(ste_temp,temp_char,len(ste_temp[0])-min(len(ste_temp[0]),neigh))
                    if bool_temp.any()==True:
                        lvl_disc_E[(t,)+curr]=-R*T*np.log(np.exp(-lvl_disc_E[(t,)+curr]/(R*T))+np.sum(W_temp[bool_temp==True]))
                        tot+=len(W_temp[bool_temp==True])

        if prun_during and t!=0 and math.isinf(get_G_sum(lvl_disc_E[t].flatten(),T)) and not (math.isinf(max_diff) and max_diff>0):
            print("MEDOC failed due to lack of machine precision. Try a 64-bit python environment.")
            quit(1)

    if per_res:
        return lvl_context0,lvl_context1
    else :
        return lvl_ste,lvl_E,lvl_disc_E

def write_file(file_name,W,silent=True):
    try :
        with open(file_name,'w') as f :
            f.writelines(W)
            f.close()
        if not silent :
            print(file_name+" written.")
    except :
        print("Could not open "+file_name)




def get_microstates_population_mesostates(DFs,DFs_err,T):
    R=1.987*10**(-3)
    W=[]
    W_err=[]    
    for i in range(len(DFs)):
       
        W+=[np.exp(-DFs[i]/(R*T))]
        W_err+=[(np.exp(-DFs[i]/(R*T))*abs(-DFs_err[i]/(R*T)))]
    p=[]
    p_err=[]
    W_tot_err=np.sqrt(sum([W_err[i]**2 for i in range(len(W_err))]))

    for i in range(len(DFs)):
        p+=[W[i]/sum(W)]
        p_err+=[p[-1]*np.sqrt((W_err[i]/W[i])**2+(W_tot_err/sum(W))**2)]

    return p,p_err

def check_and_create_rep(directory):
    if len(directory.split('/')[-1].split('.'))>1:
        new_dir='./'
        for i in range(len(directory.split('/'))-1):
            new_dir+=directory.split('/')[i]+'/'
        directory=new_dir
    if not os.path.exists(directory):
        os.makedirs(directory)    

def get_mesostate_site_spe_G(lvl_ste,lvl_E,T):
    out0=[]
    out1=[]
    
    
    for i in range(len(lvl_ste)):
        temp0=[[] for j in range(len(lvl_ste[0][0]))]
        temp1=[[] for j in range(len(lvl_ste[0][0]))]
        for j in range(len(lvl_ste[i])):
            for k in range(len(lvl_ste[i][j])):
                if lvl_ste[i][j][k]=='0':
                    temp0[k]+=[lvl_E[i][j]]
                elif lvl_ste[i][j][k]=='1':
                    temp1[k]+=[lvl_E[i][j]]
        tmp1=[]
        tmp0=[]
        for j in range(len(temp1)):
            tmp1+=[get_G_sum(temp1[j],T)]
            tmp0+=[get_G_sum(temp0[j],T)]
        out1+=[tmp1]
        out0+=[tmp0]
    out1=np.array(out1)
    out0=np.array(out0)
    return out0,out1


def plot_proba_F_per_res(save0,save1,pH,T):
    proba0=np.zeros((len(save0[0]),len(pH)),dtype=np.double)
    proba1=np.zeros((len(save0[0]),len(pH)),dtype=np.double)
    
    for k in range(len(save0[0])):
        W0=np.zeros((len(save0),len(pH)),dtype=np.double)
        W1=np.zeros((len(save0),len(pH)),dtype=np.double)
        for i in range(len(save0)):
            W0[i]=np.exp(-(save0[i,k]-i*np.log(10.)*R*T*pH[:])/(R*T))
            W1[i]=np.exp(-(save1[i,k]-i*np.log(10.)*R*T*pH[:])/(R*T))
            
            




        proba1[k]=np.sum(W1,axis=0)/(np.sum(W1,axis=0)+np.sum(W0,axis=0))
        proba0[k]=np.sum(W0,axis=0)/(np.sum(W1,axis=0)+np.sum(W0,axis=0))

    return proba1,proba0

def compute_charge_density(frac,charge,pH,seq,indices,suffix,win_sz=3,get_back=False):
    import matplotlib.pyplot as plt
    for i in range(len(seq)):
        if seq[i]=='R':
            charge[i]=1
    indices=np.array(indices)
    out=np.zeros((len(seq),len(pH)))
    for p in range(len(pH)):
        for r in range(len(seq)):
            temp=0.
            eff_len=0
            for n in range(max(0,r-win_sz),min(len(seq)-1,r+win_sz)+1):
                if np.any(n==indices): 
                    ind=np.argwhere(n==indices)[0]
                    if charge[n]>0 :
                        temp+=frac[ind,p]*charge[n] 
                    else : 
                        temp+=(1-frac[ind,p])*charge[n]
                else : 
                    temp+=charge[n]
                eff_len+=1.
            
            out[r,p]=temp/eff_len 
    
    seq_label=make_charge_color_string_for_plt(seq,[1 for i in range(len(seq))])
    plt.figure(figsize=(len(seq)/6,6))
    plt.imshow(out.transpose(),vmin=-1,vmax=1,cmap='bwr_r',aspect='auto',origin='lower')
    temp=[i+1 for i in range(len(seq)-2)]
    plt.xticks(temp,seq_label)
    temp=[i for i in range(0,len(pH),int(len(pH)/7))]+[len(pH)-1]
    temp2=[int(pH[i]) for i in range(0,len(pH),int(len(pH)/7))]+[14]
    plt.ylabel('pH') 
    plt.yticks(temp,temp2)
    plt.xlim(0.5,len(seq)-1.5) 

    cbar=plt.colorbar()
    cbar.set_label('Charge density',rotation=270)
    check_and_create_rep('./Results/Plots/')
    plt.tight_layout()
    plt.savefig('./Results/Plots/Charge_density_vs_pH_'+suffix+'.pdf')
    plt.close()
    if get_back :
        return out.transpose()

def make_charge_color_string_for_plt(seq_1,seq_id,include_caps=False,asString=False):
    if include_caps :
        st=0
        end=len(seq_id)
    else :
        st=1
        end=len(seq_id)-1
        
    seq_id_out=[]
    
    for k in range(st,end):

        if str(seq_id[k])=='1' :
            tmp=seq_1[k].upper()
            if seq_1[k]=='E' or seq_1[k]=='D':
                seq_id_out+=[r'{\textcolor{red}{'+tmp+'}}']
            elif seq_1[k]=='K' or seq_1[k]=='R' or seq_1[k]=='H':
                seq_id_out+=[r'{\textcolor{blue}{'+tmp+'}}']
            else : 
                seq_id_out+=[r'{\textcolor{black}{'+tmp+'}}']
        elif str(seq_id[k])=='2' :
           tmp=seq_1[k].lower()
           seq_id_out+=[r'{\textcolor{black}{\textbf{'+tmp+'}}}']
        elif seq_1[k]=='R' :

           tmp=seq_1[k].upper()
           seq_id_out+=[r'{\textcolor{blue}{'+tmp+'}}']
        else : 
            tmp=seq_1[k].upper()
            seq_id_out+=[r'{\textcolor{black}{\textbf{'+tmp+'}}}']
    
    if asString :
        new=''
        for i in range(len(seq_id_out)):
            new+=seq_id_out[i]
        seq_id_out=new
    return seq_id_out 
        
def plot_probas(Proba,Proba_err,st_count,pH,title='Proba',subrep='./',separate=0,print_micro=0,
                labely='Mesostate probability',color='rainbow',labels=[],publi_figure=False):
    
    import matplotlib
    from matplotlib import cm
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"]="Times New Roman"
    from matplotlib.backends.backend_pgf import FigureCanvasPgf
    matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
    
    
    pgf_with_latex = {
        "text.usetex": True,
        "pgf.preamble":
            r'\usepackage{color}'
        ,
        "font.family": "Times New Roman"
    }

    import matplotlib
    matplotlib.rcParams.update(pgf_with_latex)
    import matplotlib.pyplot as plt
    from matplotlib import cm
    if color=='rainbow':
        colors=[cm.jet_r(float(i)/len(Proba)) for i in range(len(Proba))]
    elif type(color)==type([]):
        colors=[color[i] for i in range(len(Proba))]
    else : 
        colors=[color for i in range(len(Proba))]
    if len(labels)==0:
        labels=np.array(['' for i in range(len(Proba))])

    
    
    

    plt.close()
    linewidth=1.0
    if publi_figure:

        fig=plt.figure(figsize=(4.5,2.5))
    else :
        fig=plt.figure(figsize=(7.5,3.5))

    plt.xlim(0,14)
    plt.gca().set_ylim(top=1.)    
    plt.xlabel("pH")
    plt.ylabel(labely)
    for i in range(len(Proba)):
        low_err=[]
        high_err=[]
        if len(Proba_err[i])!=0:
            for p in range(len(Proba[i])):

                if math.isinf(Proba_err[i][p]):
                    low_err+=[Proba[i][p]-Proba_err[i][p]]
                    high_err+=[Proba[i][p]+Proba_err[i][p]]
                else :
                    low_err+=[Proba[i][p]]
                    high_err+=[Proba[i][p]]
        else :
            low_err+=Proba[i]
            high_err+=Proba[i]
        if len(Proba[i])==0 :
            continue
      
        plt.plot(pH,Proba[i],color=colors[i],linestyle='-',label=labels[i],linewidth=linewidth)
        
        if print_micro==1:
            for j in range(len(st_count[i])):
                plt.plot(np.array(pH),np.array(Proba[i])*st_count[i][j],color='k',linestyle='--',linewidth=0.3)        





    check_and_create_rep('Results/Plots/')
    check_and_create_rep('Results/Plots/'+subrep)  
    plt.legend()
    plt.tight_layout()
    plt.savefig('Results/Plots/'+subrep+'/'+title+'.pdf')
    plt.close()

def plot_frac_and_deriv(pH_a,addi,name,points=[[],[]],points_syn=[[],[]],loc='./',prefix='',hill_param=[],fit=True,plot=True):
    import matplotlib.pyplot as plt
    import scipy 
    
    pH_ind_pKa=np.argwhere(abs(addi[:]-0.5)==np.amin(abs(addi[:]-0.5)))[0][0]
    pKa_a=pH_a[pH_ind_pKa]
    
    xd,d_addi=deriv(pH_a,addi)
    unsh=(10**(-pKa_a+pH_a))/(1+10**(-pKa_a+pH_a))
    xd,d_unsh=deriv(pH_a,unsh)
    if plot :
        plt.close()
        fig=plt.figure(1)
        fig.add_subplot(2,1,1)

    if hill_param!=[] and plot:
        y=hill_equation_from_hill_param(pH_a,hill_param)
        inds_lin=get_pH_equivalence(points[0],pH_a)
        corr=np.corrcoef(y[inds_lin[:]],points[1])
        plt.plot(pH_a,y,color='r',alpha=0.5,label='Hill fit (Payliss et al.) (r$^2$='+str(np.round(corr[0][1]**2,5))+')')
    height_u=np.amax(d_unsh)
    height_a=np.amax(d_addi)
    cross_height=0.0
    conf1=0.01
    conf2=1-conf1
    cross_m_ind=np.argwhere(abs(addi[:]-conf1)==np.amin(abs(addi[:]-conf1)))[0][0]
    cross_p_ind=np.argwhere(abs(addi[:]-conf2)==np.amin(abs(addi[:]-conf2)))[0][0]
   
    median_trans_ind=int(np.round((cross_p_ind+cross_m_ind)/2.,0))

    
    w_d=2. 
    
    for p in range(len(xd)):
        if xd[p]<pKa_a-w_d:
            min_ind=p+1 
            continue
        
        elif xd[p]>pKa_a+w_d:
            max_ind=p     
            break
    try :
        mean,std,asy=get_skewness(xd[min_ind:max_ind],d_addi[min_ind:max_ind])
    except :
        print("Could not compute skewdness. Try incresing the pH range.")
        exit(1)
    coop=height_a/height_u
    W_title=name+'\nCooperativity factor : '+str(np.round(coop,2))+'\nAsymmetry : '+str(np.round(asy,3))
    arg_addi=np.argwhere(d_addi==np.amax(abs(d_addi)))[0][0]

    if len(points[0])>0:
        inds_lin=get_pH_equivalence(points[0],pH_a)  
        corr=np.corrcoef(addi[inds_lin[:]],points[1])
    else : 
        corr=np.zeros((2,2))*np.double('NaN')
    if plot :
        plt.title(W_title)
        plt.subplot(211)
        if fit:
            plt.plot(pH_a,addi,color='g',label='Fit (r$^2$='+str(np.round(corr[0][1]**2,5))+')')
        else :
            plt.plot(pH_a,addi,color='g',label='Prediction')
        if len(points_syn[0])!=0:
            plt.scatter(points_syn[0],points_syn[1],color='None',edgecolors='grey',label='Synthetic points')
        if len(points[0])!=0:
            plt.scatter(points[0],points[1],color='k',label='Real points')
        
        plt.xlim(0,14)
        plt.ylabel('Fraction deprotonated')
        plt.legend()
        plt.subplot(212)
        plt.xlim(0,14)
        plt.xlabel('pH')
        plt.ylabel('Fraction deprotonated\nfirst derivative')

        if len(points[0])>0:
            d_points_pH,d_points=deriv(points[0],points[1])

            plt.scatter(d_points_pH,d_points,color='k')
        if hill_param!=[]:
            xd,dy=deriv(pH_a,y)
            plt.plot(xd,dy,color='r',alpha=0.5)
        plt.vlines(pKa_a,0,d_addi[pH_ind_pKa],color='orange',linestyle='--',label='pK$_a$')

        
        
        plt.plot(xd,d_addi,color='g')
        plt.legend()
        plt.savefig(loc+prefix+'Fraction_protonated_derivative_'+name+'.pdf')
        plt.close()

    
    
    
    
    

    

    
    
    return pKa_a,coop,asy

def get_sym(pH,proba,pKa_ind):
    sym1=proba.copy()
    sym2=proba.copy()
    for p in range(len(proba)):   
        if 2*pKa_ind-p<0:
            sym1[p]=0.
            sym2[p]=0.
            continue
        elif  2*pKa_ind-p>=len(proba):
            sym1[p]=1.
            sym2[p]=1.
            continue
        if p>pKa_ind:  
            sym1[p]=1-proba[2*pKa_ind-p]
        elif p<pKa_ind:
            sym2[p]=1-proba[2*pKa_ind-p]










def deriv(x,y):                                                                                                                                                          
    return (x[1:]+x[0:-1])/2.,(y[1:]-y[0:-1])/(x[1:]-x[0:-1])

def hill_equation(pH,pKa,n,plat1,plat2): 
    y=(10**(n*(pH-pKa)))/(1.+10**(n*(pH-pKa)))*(plat2-plat1)+plat1
    return y

def hill_equation_from_hill_param(pH,hill_param):
    return hill_equation(pH,hill_param[0],hill_param[1],hill_param[2],hill_param[3])
    
def get_pH_equivalence(pH_dis,pH_lin):
    
    inds=np.zeros((len(pH_dis)),dtype=np.int)
    for i in range(len(pH_dis)):
        inds[i]=np.argwhere(np.amin(abs(pH_dis[i]-pH_lin[:]))==abs(pH_dis[i]-pH_lin[:]))[0][0]
    return inds

def get_skewness(x,y):
    
    
    N=np.sum(y)
    mean=np.sum(x*y)/N
    p=y/N
    std=np.sqrt((np.sum(p*(x-mean)**2)))
    mu3=np.sum(p*((x-mean)/(std))**3)    
    return mean,std,mu3

def mix_probas_fitting(Proba_meso,Proba,pH):
    proba_out=np.zeros((len(Proba[0]),len(pH)))
    for r in range(len(Proba[0])):
        for i in range(len(Proba_meso)):
            proba_out[r]+=Proba_meso[i]*Proba[i,r]
    proba_out=1.-proba_out                                                      
    return proba_out 

def get_mesostate_Fs(Fs,T):
    F_meso=[]
    for i in range(len(Fs)):
        offset=0
        temp=0 
            
        if np.mean(Fs[i])<-306. or np.mean(Fs[i])>306.:
            offset=np.mean(Fs[i])
        for j in range(len(Fs[i])):
            temp+=np.exp(-(Fs[i][j]-offset)/(R*T))
        F_meso+=[-R*T*np.log(temp)+offset]
    return F_meso


def get_probas_2(Fs,Fs_err,pH,T,ign_norm_err=True,unsafe=0):

    Fs=np.array(Fs)
    Fs_err=np.array(Fs_err)
    Weights_count=np.array([ len(Fs)-k-1 for k in range(len(Fs))])
    Weights=np.zeros((len(Fs),len(pH)),dtype=np.double)
    Weights_err=np.zeros((len(Fs),len(pH)),dtype=np.double)
    Weights_norm=np.zeros((len(pH)),dtype=np.double)
    Weights_norm_err=np.zeros((len(pH)),dtype=np.double)
    
    MIN=np.amin(Fs)
    MAX=np.amax(Fs)

    Exponents=np.zeros((len(Fs),len(pH)),dtype=np.double)
    Exponents2=np.zeros((len(Fs),len(pH)),dtype=np.double)
    Exponents3=np.zeros((len(Fs),len(pH)),dtype=np.double)
    Fs=np.asarray(Fs)
    Proba=np.zeros((len(Fs),len(pH)))
    Proba_err=np.zeros((len(Fs),len(pH)))
    
    for i in range(len(Fs)):
        if math.isnan(Fs[i]):
            continue
        Exponents[i,:]=(-Fs[i]+np.log(10.)*R*T*pH[:]*Weights_count[i])/(R*T)
        Exponents2[i,:]=(np.log(10.)*R*T*pH[:]*Weights_count[i])/(R*T)
        Exponents3[i,:]=-Fs[i]


    for i in range(len(Fs)):
        Weights[i,:]=np.exp(Exponents[i,:])
    
    Weights_norm[:]=np.sum(Weights[:,:],axis=0)

    for i in range(len(Fs)): 
        if math.isnan(Fs[i]):
            if unsafe==1 :
                Proba[i]+=np.zeros((len(pH)))
            continue
        for p in range(len(pH)):
            Proba[i,p]=Weights[i,p]/Weights_norm[p]        
            if ign_norm_err==False :
                if p>0 and p<len(pH)-1 :
                    cor=np.corrcoef([Weights[i,p-1:p+1],Weights_norm[p-1:p+1]])
                else :
                    cor=np.zeros(((2,2)),dtype=float)
                A=Weights_norm[p]
                S_A=Weights_norm_err[p]
                B=Weights[i,p]
                S_B=Weights_err[i,p]
                if math.isnan(cor[0][1])==True :
                    cor[0][1]=0.0

                if S_A==float('inf') or S_B==float('inf') or A==0.0 or B==0.0:
                    Proba_err[i]+=[0.0]

                else:
                    Proba_err[i]+=[abs(Proba[i,p])*\
                    np.sqrt((S_A/A)**2+(S_B/B)**2-2*cor[0][1]*S_A*S_B/(A*B))]
            else :
                Proba_err[i]+=[Weights_err[i,p]/Weights_norm[p]]

    
    for i in range(len(Fs)):
        if math.isnan(Fs[i]):
            continue
 
        Weights_err[i,:]=abs(((Fs_err[i])/(R*T)))*np.exp((-Fs[i]+np.log(10.)*R*T*pH[:]*Weights_count[i])/(R*T))
        
        Weights_norm[:]+=np.exp((-Fs[i]+np.log(10.)*R*T*pH[:]*Weights_count[i])/(R*T))
        if ign_norm_err==False :
            for p in range(len(pH)):
                if p>0 and p<len(pH)-1 :
                    cor=np.corrcoef([Weights_norm[p-1:p+1],Weights[i,p-1:p+1]])
                else :
                    cor=np.zeros(((2,2)),dtype=float)
                if math.isnan(cor[0][1])==True :
                    cor[0][1]=0.0
                if math.isnan(Fs_err[i])==False :
                    if Weights_norm_err[p]!=0:
                        Weights_norm_err[p]+=np.sqrt((Weights_norm_err[p])**2+(Weights_err[i,p])**2)
                    else :
                        Weights_norm_err[p]+=Weights_err[i,p]

    return Weights_norm,Weights_norm_err,Proba,Proba_err

def plot_probas2(Proba,Proba_err,st_count,pH,title='Proba',subrep='./',legend=[]):
    
    import matplotlib
    
    from matplotlib import cm
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"]="Times New Roman"
    from matplotlib.backends.backend_pgf import FigureCanvasPgf
    matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
    
    pgf_with_latex = {
        "text.usetex": True,
        "pgf.preamble":
            r'\usepackage{color}'
        ,
        "font.family": "Times New Roman"}
    import matplotlib
    matplotlib.rcParams.update(pgf_with_latex)
    import matplotlib.pyplot as plt
    from matplotlib import cm

    colors=[cm.gnuplot(float(i)/len(Proba)) for i in range(len(Proba))] 
    fig=plt.figure(figsize=(7,3))

    ax1=fig.add_axes([0.08,0.15,0.8,0.80])

    plt.xlim(np.round(np.amin(pH)),np.round(np.amax(pH)))
    plt.gca().set_ylim(top=1.)    
    pH=np.array(pH)

    if len(legend)==0:
        legend=['' for i in range(len(Proba))]
    plt.xlabel("pH")
    plt.ylabel("Mesostate probability")

    for i in range(len(Proba)):
        low_err=[]
        high_err=[]
        if len(Proba_err[i])!=0:
            for p in range(len(Proba[i])):
                if math.isinf(Proba_err[i][p]):
                    low_err+=[Proba[i][p]-Proba_err[i][p]]
                    high_err+=[Proba[i][p]+Proba_err[i][p]]
                else :
                    low_err+=[Proba[i][p]]
                    high_err+=[Proba[i][p]]
        else :
            low_err+=Proba[i]
            high_err+=Proba[i]
        if len(Proba[i])==0 :
            continue

        plt.plot(pH,Proba[i],color=colors[i],linestyle='-',label=str(legend[i]))
        plt.fill_between(pH, low_err, high_err, color=colors[i],alpha=0.5)
    plt.legend()
    check_and_create_rep('Results/Plots/')
    check_and_create_rep('Results/Plots/'+subrep)
    plt.savefig('Results/Plots/'+subrep+'/Mesostates_'+title+'.pdf')

def write_per_res_F(save0,save1,suffix):
    check_and_create_rep('Results/Fs/')
    W=''                                                                                                                                                                  
    for i in range(len(save0)):
        for j in range(len(save0[i])):
            W+=str(save0[i][j])+'\t'
        W+='\n'
    write_file('Results/Fs/Per_res_0_'+suffix+'.txt',W)
    W=''
    for i in range(len(save1)):
        for j in range(len(save1[i])):
            W+=str(save1[i][j])+'\t'
        W+='\n'
    write_file('Results/Fs/Per_res_1_'+suffix+'.txt',W)
 
def plot_quantity_vs_pH2(vals,Proba_1,pH,seq_1,st_count,name,publi_figure=False,labely=''):
    import matplotlib 
    matplotlib.use("pgf")
    from matplotlib.backends.backend_pgf import FigureCanvasPgf 
    matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
    
    pgf_with_latex = {
        "text.usetex": True,
        "pgf.preamble":
            r'\usepackage{color}'
        ,
        "font.family": "Times New Roman"
    }
    import matplotlib
    matplotlib.rcParams.update(pgf_with_latex)
    import matplotlib.pyplot as plt    
    from matplotlib import cm
    plt.close()
    if publi_figure:
        plt.figure(figsize=(4,3))
    else :
        plt.figure(figsize=(8,6))
    Values=np.zeros((len(pH)),dtype=float)      
    Values=np.transpose(Values)
    SSP=np.asarray(vals) 
    for i in range(len(st_count)):     
        for j in range(len(st_count[i])):
            if math.isnan(st_count[i][j]*SSP[i][j])==False:           
                Values[:]=Values[:]+Proba_1[i]*st_count[i][j]*SSP[i][j]
    plt.xlim(0,14)
    plt.ylim(np.min(Values)-0.5,np.max(Values)+0.5)
    plt.ylabel(labely)
    plt.xlabel('pH')
    plt.plot(pH,Values)
    plt.tight_layout()
    plt.savefig('Results/Plots/'+name+'_vs_pH.pdf')
    plt.close()
    save_pH_dependent_curve(Values,pH,name)   

def save_pH_dependent_curve(vals,pH,name):
    W='pH\t'+name+'\n'
    for p in range(len(pH)):
        W+=str(pH[p])+'\t'+str(vals[p])+'\n'
    check_and_create_rep('./Results/Plots/pH_plots_raw_data')
    write_file('./Results/Plots/pH_plots_raw_data/'+name+'.txt',W)

import matplotlib
import time
matplotlib.use("pgf")                                                                    
from matplotlib.backends.backend_pgf import FigureCanvasPgf 
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
pgf_with_latex = {
    "text.usetex": True,
    "pgf.preamble":
        r'\usepackage{color}',
    "font.family": "Times New Roman" }
import matplotlib
matplotlib.rcParams.update(pgf_with_latex)
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import scipy.misc
import math
import os
global AA_type
MODULE_PATH = os.path.dirname(__file__)
AA_type=['E','D','H','K','Y','8','9','R','C','A','F','G','L','I','M','N','P','Q','T','V','S','W']
default_res=0.01
R=np.double(1.98720425864083*10**(-3))

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--site_specific","-ss",
                        help="Whether the scope of the prediction is global (0) or site-specific (1)",
                        default=0)
    parser.add_argument("--seq_file","-s",
                        help="Sequence file containing a single line of uninterrupted amino acid sequence",
                        default='seq.fasta')
    parser.add_argument("--pH_range","-pH",nargs=2,default=[0,14])
    parser.add_argument("--max_frac", "-mf",
                        help="Max fraction represented by discarded states",default=9999)
    parser.add_argument("--nneigh", "-nn",
                        help="Number of immediate neighbors (on one side) for peptides",default=2)
    parser.add_argument("--predict_type", "-pt",
                        help="Prediction type ([I]mplicit, [E]xplicit, [U]nshifted)",
                        default="I")
    parser.add_argument("--plot_res", "-pr",
                        help="Plot resolution",
                        default=default_res)
    parser.add_argument("--det_level", "-dl",
                        help="Level of details for output ",
                        default=2)
    parser.add_argument("--temp", "-t",
                        help="Temperature ",
                        default=298.)
    args = parser.parse_args()
    if args.pH_range :
        if len(args.pH_range)!=2:
            print("Must provide 2 values for the pH range")
            exit()
        else :
            pH_range=[float(args.pH_range[0]),float(args.pH_range[1])]
    else :
        pH_range=[0,14]
        pH_range=[0,14]

    if args.seq_file:
        seq_file=args.seq_file

    if args.site_specific:
        try:
            if int(args.site_specific)==1:
                per_res=True
            elif int(args.site_specific)==0:
                per_res=False
            else :
                per_res=False
                print("Invalid option for site_specific, defaulting to ",per_res)
        except :
            per_res=False
            print("Invalid option for site_specific, defaulting to ",per_res)
    else :
        per_res=False

    if args.temp:
        try: 
            T=np.double(args.temp)
        except :
            T=np.double(298.)
            print("Invalid option for det_level, defaulting to ",T)

    if args.det_level:
        try: 
            det_lvl=int(args.det_level)
        except :
            det_lvl=2
            print("Invalid option for det_level, defaulting to ",det_lvl)

    if args.plot_res:
        try: 
            res=float(args.plot_res)
        except : 
            res=default_res
            print("Invalid option for plot res, defaulting to ",res)
    else : 
        res=default_res

    if args.predict_type :
        if args.predict_type.upper()=='I' :
            additive=True
            penta=False 
            unshifted=False
            suffix='Implicit'
        elif args.predict_type.upper()=='U' :
            additive=False
            penta=False
            unshifted=True
            suffix='Unshifted'
        else :
            print("Invalid option for predict type, aborting")
            quit()
    else : 
        additive=False
        penta=True
        unshifted=False

    prun_during=True
    prun_after=False
    if args.nneigh :
        try : 
            if int(args.nneigh)>=1:
                neigh=int(args.nneigh)
            else :
                neigh=2
                print("Invalid value for nneigh, defaulting to ",neigh)
        except : 
            neigh=2
            print("Invalid value for nneigh, defaulting to ",neigh)
    fraction_kept=0,
    max_frac=float('+inf')
    test_time=True
    publi_figure=True
    reduced=False
    test=''
    recal=''
    HT=True
    debug=False
    reverse=True
    base_E=-709.0*(R*T)
    if HT==True :
        database_suffix='_HT'     
    else : 
        database_suffix=''

    base_rep=''
    initialize_additive_DF_array(database_suffix)
    pH=np.arange(pH_range[0],pH_range[1],res,dtype=np.double)
    
    
    seq_1=read_file(seq_file)[0][0]
    if seq_1[0]!='z':
        seq_1='z'+seq_1
    if seq_1[-1]!='z':
        seq_1=seq_1+'z'

    seq_3=convert_AA_1_to_3(seq_1,mode=2)
    list_res=[i for i in range(len(seq_3))]
    

    T_expt=298.
    seq_1=convert_AA_3_to_1_letter(seq_3)
    max_diff=-R*T*np.log(max_frac)

    t_start=time.time()
    sites_num,titrable_residue_indexes,pos_res,neg_res,base_charge,arg_res,raw_seq,seq_data_q,seq_data_id,seq_id_reduced,map,new_W=read_sequence(seq_3,list_res)
    Nmes=pos_res+base_charge+1-(-neg_res+arg_res+base_charge)
    layers_q=np.array([i for i in range(-neg_res+arg_res+base_charge,pos_res+base_charge+1)])
    refs=get_ref_pkas(1,DF_only=False,suffix=recal+test)

    contexts,states=get_base_contexts(map,neigh,seq_1,seq_data_q)
    pat_E,pat_ste,pat_seq=get_all_contexts(contexts,states,neigh,refs,unshifted,reverse,penta,additive,T,base_rep,reduced)
    if per_res :
        lvl_context0,lvl_context1=main_prediction(T,neigh,contexts,pat_E,pat_ste,pat_seq,map,prun_during,per_res,base_E,test_time,max_diff)
    else :
        lvl_ste,lvl_E,lvl_disc_E=main_prediction(T,neigh,contexts,pat_E,pat_ste,
                                                                                         pat_seq,map,prun_during,
                                                                                         per_res,base_E,test_time,
                                                                                         max_diff)

    disc_E=[]
    print("There are ",pos_res-arg_res+neg_res," ionizable residue")
    if test_time :
        print("MEDOC ran in "+str(time.time()-t_start)+" seconds")
        write_file('./Compute_time.txt',str(time.time()-t_start))
    if not per_res:
        SUM=[]
        for i in range(Nmes):
            SUM+=[len(lvl_ste[i])]
    if (prun_during or prun_after) and not per_res:
        if fraction_kept==1. :
            SW=''
            EW=''
            popW=''
            for i in range(len(lvl_ste)):
                for j in range(len(lvl_ste[i])):
                    SW+=lvl_ste[i][j]+'\t'
                    EW+=str(lvl_E[i][j])+'\t'
                SW+='\n'
                EW+='\n'
                temp=get_microstates_population_mesostates(lvl_E[i],lvl_E[i]*0,T)
                for j in range(len(temp[0])):
                    popW+=str(temp[0][j])+'\t'
                popW+='\n'
            print('./Results/Fs/Fs_'+suffix+'.txt')
            write_file('./Results/Fs/Populations_'+suffix+'.txt',popW)
            check_and_create_rep('./Results/States_details/')
            write_file('./Results/States_details/States_'+suffix+'.txt',SW)
            write_file('./Results/Fs/Fs_'+suffix+'.txt',EW)

    Meso_G_all=np.zeros((Nmes),dtype=np.double)

    if prun_during and not per_res:
        if det_lvl>=2: 
            print("q\tTotal energy\tKept energy\tDiscarded energy")
        for i in range(len(lvl_disc_E)):
            temp_1=get_G_sum(np.concatenate((lvl_disc_E[i].flatten(),lvl_E[i])),T)
            temp_2=get_G_sum(lvl_E[i],T)
            temp_3=get_G_sum(lvl_disc_E[i].flatten(),T)
            if det_lvl>=2:
                print(layers_q[len(layers_q)-1-i],'\t',"%0.4f" % temp_1,'\t',"%0.4f" % temp_2,'\t',"%0.4f" % temp_3)
    elif not per_res :
        if det_lvl>=2:
            print("Total energy\tKept energy\tDiscarded energy")
        for i in range(Nmes): 
            temp_1=get_G_sum(lvl_E[i],T)
            temp_2=disc_E[i] 
            temp_3=get_G_sum([temp_1,temp_2],T)
            if det_lvl>=2:
                print(temp_3,temp_1,temp_2)

    check_and_create_rep('Results/Fs')
    if per_res:
        if debug==True:
            reconstructed=np.zeros((2,)+np.shape(lvl_context0))
            reconstructed[:,:,:]=float('+inf')
            for i in range(Nmes):
                for j in range(len(lvl_ste[i])):
                    ind=()        
                    for k in range(len(lvl_ste[i][j])):
                        if lvl_ste[i][j][k]=='0':
                            base=(0,i,k)
                        else :
                            base=(1,i,k)
                        for l in range(len(lvl_ste[i][j])-neigh,len(lvl_ste[i][j])):
                            base+=(int(lvl_ste[i][j][l]),)
                        ind+=(base,)
                    for k in range(len(ind)):             
                        reconstructed[ind[k]]=get_G_sum([reconstructed[ind[k]],lvl_E[i][j]],T)
                if det_lvl>=3 :                
                    for k in range(len(reconstructed[0][i])):
                        print(i,k)
                        print("0should be ")
                        print(reconstructed[0][i][k])
                        print("0is ")
                        print(lvl_context0[i][k])
                        print()

            out0,out1=get_mesostate_site_spe_G(lvl_ste,lvl_E,T)
        save0=np.zeros((Nmes,len(lvl_context0[0])))
        save1=np.zeros((Nmes,len(lvl_context0[0])))

        for i in range(Nmes):
            for k in range(len(lvl_context0[0])):
                save0[i,k]=get_G_sum(lvl_context0[i,k].flatten(),T)
                save1[i,k]=get_G_sum(lvl_context1[i,k].flatten(),T)
        save0_tmp=save0
        save1_tmp=save1
        print("Computing residue wise probabilities")

        proba0,proba1=plot_proba_F_per_res(save0-2*base_E,save1-2*base_E,pH,T)

        if det_lvl>=1:
            
            print("Plotting charge density vs pH")
            compute_charge_density(proba0,seq_data_q,pH,seq_1,map,suffix)
            
            

        print("Plotting residue wise probabilities")
        IAAs=['D','E','Y','H','K']

        colors=['orange','red','maroon','darkturquoise','blue']
        st_count_fake_r=[[1.] for i in range(len(proba0))]
        temp_all=[]
        temp_temp=[]
        colors_all=[]
        all_colors_all=[]
        all_labels=[]
        all_labels_all=[]
        used=np.zeros(np.shape(IAAs))
        for t in range(len(IAAs)):
            temp=[]
            temp_labels=[]
            temp_colors=[]
            for i in range(len(proba0)):
                if seq_1[map[i]]==IAAs[t]:
                    temp+=[proba0[i]]
                    temp_all+=[proba0[i]]
                    colors_all+=[colors[t]]
                    temp_colors+=[colors[t]]
                    if used[t]==0:
                        used[t]=1
                        temp_labels+=[IAAs[t]]
                        all_labels+=[IAAs[t]]
                    else :
                        all_labels+=['_'+IAAs[t]]
                        temp_labels+=['_'+IAAs[t]]
            temp_temp+=[temp]
            all_labels_all+=[temp_labels]
            all_colors_all+=[temp_colors]
        if det_lvl>=2:
            try :
                print("Plotting all residue probabilities")
                plot_probas(np.array(temp_all),np.array(temp_all)*0.,st_count_fake_r,pH,
                           title='Proba_Residue_all_'+suffix,labely='Fraction protonated',
                           color=colors_all,labels=np.array(all_labels),publi_figure=publi_figure)
            except :
                print('Could not plot all probas at once')

        if det_lvl>=3:
            for t in range(len(IAAs)):
                if len(temp_temp[t])==0:
                    continue
                try :
                    print("Plotting all "+IAAs[t]+" probabilities")
                    plot_probas(np.array(temp_temp[t]),np.array(temp_temp[t])*0.,st_count_fake_r,pH,
                                   title='Proba_Residue_all_'+IAAs[t]+'_'+suffix,labely='Fraction protonated',
                                   color=all_colors_all[t],labels=all_labels_all[t],publi_figure=publi_figure)
                except :
                    print("Could not plot probas for "+IAAs[t])

    if det_lvl>=2 and per_res:
        print("Printing site specific derivative and transition parameters")
        param_all=np.zeros((len(proba0),3),dtype=np.double)
        names=[]

        for i in range (len(proba0)):
            names+=[seq_1[map[i]]+str(map[i]+1)]
            print(names[i])
            if det_lvl<4 :
                param_all[i]=plot_frac_and_deriv(pH,proba1[i,:],name=seq_1[map[i]]+str(map[i]+1),loc='./Results/Plots/',fit=False,plot=False)
            else :
                param_all[i]=plot_frac_and_deriv(pH,proba1[i,:],name=seq_1[map[i]]+str(map[i]+1),
                                                    loc='./Results/Plots/',fit=False,plot=True)

        x=[i for i in range(len(param_all))]
        plt.bar(x,param_all[:,1],color='orange')

        plt.xticks(x,names,rotation=90)
        plt.savefig('Results/Plots/Cooperativity_factor_per_res_pred.pdf')
        plt.close()
        W=''
        for i in range(len(names[i])):
            W+=names[i]+'\t'+str(param_all[i,1])+'\n'
        write_file('./Results/Plots/pH_plots_raw_data/Cooperativity_factor_per_res.txt',W)

        plt.bar(x,param_all[:,2],color='green')
        plt.xticks(x,names,rotation=90)
        plt.savefig('Results/Plots/Asymetry_factor_per_res_pred.pdf')
        plt.close()
        W=''
        for i in range(len(names[i])):
            W+=names[i]+'\t'+str(param_all[i,2])+'\n'
        write_file('./Results/Plots/pH_plots_raw_data/Asymetry_factor_per_res.txt',W)






    if per_res:
        W='pH\t'
        for i in range(len(temp_all)):
            W+='res_'+str(map[i])+'\t'
        W+='\n'

        for p in range(len(pH)):
            W+=str(pH[p])+'\t'
            for r in range(len(temp_all)):
                W+=str(temp_all[r][p])+'\t'
            W+='\n'
        check_and_create_rep('./Results/Plots/pH_plots_raw_data/')
        write_file('./Results/Plots/pH_plots_raw_data/Fraction_protonated_per_residue_'+suffix+'.txt',W)

        W=''
        pKas_out=[]
        pKas_ind=[]
        for i in range(len(proba0)):
            for p in range(1,len(proba0[i])):
                if proba0[i][p]<0.5  and proba0[i][p-1]>=0.5 :
                    pKas_out+=[pH[p]]
                    pKas_ind+=[p]
                    W+=seq_1[map[i]]+str(map[i])+'\t'+str(np.round(pH[p],len(str(res))-1))+'\n'
                    break
                elif p==len(proba0[i])-1:
                    pKas_out+=[float('NaN')]
                    pKas_ind+=[p]
                    W+=seq_1[map[i]]+str(map[i])+'\t'+str('NaN')+'\n'
        check_and_create_rep('Results/pKas/')
        write_file('Results/pKas/pKas_per_res_'+suffix+'.txt',W)

    pop_count=[[1.] for i in range(Nmes)]
    qs=[[layers_q[i]] for i in range(len(layers_q))]

    if det_lvl>=5 and max_frac==0.0:
        F_meso=get_mesostate_Fs(lvl_E,T) 
        print("Computing weights restricted")
        Weights_norm,Weights_norm_err,Proba_meso,Proba_meso_err=get_probas_2(F_meso,F_meso,pH,T)
        print("Plotting probas restricted")
        plot_probas2(Proba_meso,Proba_meso_err,pop_count,pH,subrep='',title='Proba_restricted_'+suffix)

    if per_res :
        for i in range(Nmes):
            Meso_G_all[i]=get_G_sum([get_G_sum(lvl_context0[i,0].flatten(),T),get_G_sum(lvl_context1[i,0].flatten(),T)],T)
        write_per_res_F(save0,save1,suffix)
    else :
        for i in range(Nmes):
            
            Meso_G_all[i]=get_G_sum(np.concatenate((lvl_disc_E[i].flatten(),lvl_E[i])),T)

    Meso_G_all=Meso_G_all[::-1]-2*base_E  

    W=''
    for i in range(len(Meso_G_all)):
        W+=str(Meso_G_all[i])+'\n'
    check_and_create_rep('Results/Fs')
    write_file('Results/Fs/Mesostate_F_'+suffix+'.txt',W,silent=False)

    if det_lvl>=1:
        print("Computing weights")
        Weights_norm,Weights_norm_err,Proba_meso,Proba_meso_err=get_probas_2(Meso_G_all,Meso_G_all*0.,pH,T)
    if det_lvl>=2:
        print("Plotting mesostate probas")
        plot_probas2(Proba_meso,Proba_meso_err,pop_count,pH,subrep='',title='Proba_all_'+suffix,legend=layers_q)
    if det_lvl>=1:
        print("Plotting the charge profile")
        check_and_create_rep('Results/Plots')
        plot_quantity_vs_pH2(qs,Proba_meso,pH,seq_1,pop_count,'Net_charge_'+suffix,publi_figure=publi_figure,labely='Net charge')
