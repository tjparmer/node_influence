#GTN construction
#Thomas Parmer, 2023

import networkx as nx
import sympy as sp
from sympy.logic.boolalg import to_dnf,is_dnf
import re
from cana import BooleanNetwork #only necessary for file_to_pi and alter_network functions

illegals={'/',';','-','+','*','.'}
bools={'AND','OR','NOT'}

#create expression directly from look-up table; takes a BooleanNode as input
#NOTE: this isn't necessary for GTN construction but may be helpful to understand logical rules
def expr_from_LUT(node,match=1):
    """ takes a BooleanNode as input; only consider rows that result in output=match """
    lut=node.look_up_table()
    rhs=[]
    for i,row in enumerate(lut['In:']):
        if lut['Out:'][i]!=match: continue #we only need to consider the cases that match the state we want
        inputs=[node.inputs[j] if lut['In:'][i][j]=='1' else 'not '+node.inputs[j] for j,inp in enumerate(row)]
        rhs.append('('+' and '.join(inputs)+')')
        #print row,lut['Out:'][i],' and '.join(inputs)
    return node.name+'*='+' or '.join(rhs)


#turn a BooleanNetwork into a general threshold network based on node LUTs
#threshold networks are digraphs that need variable, state, type, tau, delay, label (what displays as the node label)
def create_threshold_network_LUTs(n, single_nodes=[]):
    """ Create a digraph compatible with threshold network manipulation
    Expects a BooleanNetwork as input, returns a threshold-like digraph where every node is a literal and/or threshold node
    
    Valid nodes include these properties: label, threshold (tau), type, time delay, variable, state
    Valid edges include these properties: type 
    
    Set single_nodes to specify nodes that should NOT have an opposing state (e.g. non-binary nodes) """
    
    tn = nx.DiGraph()
    tn.add_nodes_from([node.name+'-1' for node in n.nodes],state=1) #add positive s-units as nodes
    tn.add_nodes_from([node.name+'-0' for node in n.nodes if node.name not in single_nodes],state=0) #add negative s-units as nodes
    #add node properties
    for name in tn.nodes():
        tn.node[name]['delay']=1
        tn.node[name]['tau']=0
        tn.node[name]['type']='variable'
        tn.node[name]['variable']=name[:-2]
        tn.node[name]['label']=name[:-2]
    
    #add threshold nodes/edges by adding one per LUT entry
    #NOTE: if the node doesn't have an opposing state, ignore all LUT entries that result in 0
    for node in n.nodes:
        #print node.name,node.look_up_table()
        for i,row in enumerate(node.look_up_table()['In:']): 
            if node.name in single_nodes and node.look_up_table()['Out:'][i]==0: continue #ignore the opposing state
            inputs=[inp+'-'+str(row[j]) for j,inp in enumerate(node.inputs)]
            not_in_singles = lambda x: x[-1]!='0' or x[:-2] not in single_nodes #filter out opposing states
            inputs = filter(not_in_singles, inputs) #solution from https://stackoverflow.com/questions/4915920/how-to-delete-an-item-in-a-list-if-it-exists
            length=len(inputs)
            if not length: continue #this row was filtered out
            output=node.name+'-'+str(node.look_up_table()['Out:'][i])
            name='T-'+str(i)+'_'+output
            #print row,node.look_up_table()['Out:'][i],inputs,output,name,length
            tn.add_node(name,delay=0,tau=length,type='threshold',label=length,group=node.name) #counter increases regardless of state
            for inp in inputs: tn.add_edge(inp,name) #add edges to threshold unit
            tn.add_edge(name,output) #add edge from threshold unit to output s-unit
    
    return tn


#convert Cell Collective expr file to standard Boolean text format
#additionally, add in external components with self copy functions
#NOTE: this copies directly without trying to get rid of illegal characters or other issues
def expressions_to_boolean_text(filename,write_file=None,external_file=None):
    lines=""
    with open(filename) as f:
        for line in f:
            if '=' in line:
                lines+=line.split('=')[0].strip()+'*= '+line.split('=')[1].strip()+'\n'
    if external_file:
        with open(external_file) as f:
            for line in f:
                lines+=line.strip()+'*='+line.strip()+'\n'
    #lines=lines.lower() #python only recognizes lower case as logical operators
    for b in bools:
        lines=lines.replace(' '+b+' ',' '+b.lower()+' ')
    if write_file:
        with open(write_file,'w') as f:
            f.write(lines)
    return lines


#convert Cell Collective expr file to standard Boolean text format using ASCII values for illegal characters
#additionally, add in external components with self copy functions
def expressions_to_boolean_text_ascii(filename,write_file=None,external_file=None):
    lines=""
    with open(filename) as f:
        for line in f:
            if '=' in line:
                for illegal in illegals:
                    if illegal in line: line=line.replace(illegal,'a'+str(ord(illegal))+'a')
                #print line
                lines+=line.split('=')[0].strip()+'*= '+line.split('=')[1].strip()+'\n'
    if external_file:
        with open(external_file) as f:
            for line in f:
                for illegal in illegals:
                    if illegal in line: line=line.replace(illegal,'a'+str(ord(illegal))+'a')
                lines+=line.strip()+'*='+line.strip()+'\n'
    for b in bools:
        lines=lines.replace(' '+b+' ',' '+b.lower()+' ')
    #special cases:
    if "PC12 Cell Differentiation" in filename:
        lines=lines.replace('G(ia47ao)','Ga40aia47aoa41a') #replace only these parentheses
    if "Signaling Pathway for Butanol Production in Clostridium beijerinckii NRRL B-598" in filename:
        lines=lines.replace('NAD(P)H','NADa40aPa41aH')
        lines=lines.replace('3a45ahydroxybutyryla45aCoA','a51aa45ahydroxybutyryla45aCoA')
    if "TOL Regulatory Network" in filename: #replace leading number
        lines=lines.replace('3MBz','a51aMBz')
    if "Signal Transduction in Fibroblasts" in filename:
        lines=lines.replace('Palpha_iR','PalphaiR') #causes issue with DNF variable replacement
        lines=lines.replace('Palpha_1213R','Palpha1213R')
        lines=lines.replace('Palpha_sR','PalphasR')
    if "SKBR3 Breast Cell Line Long-term ErbB Network" in filename:
        lines=lines.replace('ERBB2','ER_BB2') #causes issue with DNF variable replacement
    if "Signaling Pathway for Butanol Production in Clostridium beijerinckii NRRL B-598" in filename:
        lines=lines.replace("acetaldehyde","acetal_dehyde") #causes issue with DNF variable replacement
    if "IL-6 Signalling" in filename:
        lines=lines.replace("camk24","camk_24") #causes issue with DNF variable replacement
        lines=lines.replace("dum_gab1_kin_or_jak1_gab1_mem_p","dum_gab1_kin_jak1_gab1memp")
        lines=lines.replace("dum_il6rc_p_or_grb2_vav","dum_il6rc_p_grb2_vav")
        lines=lines.replace("dum_mtorc1_or_pkcd_stat3_ta","dum_mtorc1_pkcd_stat3_ta")
    if "Signaling in Macrophage Activation" in filename:
        lines=lines.replace("BAG4_TNFRSF1A","BAG4_TN_FRSF1A")
        lines=lines.replace("TNF_BAG4_TN_FRSF1A","TN_BAG4_TN_RSF1A")
    if "T Cell Receptor Signaling" in filename:
        lines=lines.replace("vav1","vaa118a1")
    if write_file:
        with open(write_file,'w') as f:
            f.write(lines)
    return lines


#NOTE: this function breaks if the node name has a comma or other punctuation in it, e.g. Clb1,2
#NOTE: this function breaks if a node has more than 36 inputs
#TODO: make this more robust to bad naming to avoid manual fixes
def convert_dnf(expr,nodes):
    """ convert logical expresion to disjunctive normal form; returns both the dnf and the dnf of the expression's negation 
    expects expr in a format such as PH_P1*=PTC_P1 and (not HH_P4 or not HH_P2)
    expects nodes as list such as ['PTC_P1','HH_P4','HH_P2']
    nodes is required as a list of all variable names that may occur in the expression """
    
    expr=expr.split('*=')[1]
    literals=[node for node in nodes if node in re.findall(r'\w+',expr)] #variables present in the expression; use regex to separate words
    literals=sorted(literals,reverse=True) #reverse sort to avoid replacing partial names
    variables=['v'+str(i) if i<10 else 'v'+chr(i+55) for i in range(len(literals))] #over 10 variables, must use letters; input limit=36
    #print literals,variables
    for i,node in enumerate(literals): #convert to generic variables that work with exec and eval
        expr=expr.replace(node,variables[i])
    exec(','.join(variables) + " = sp.symbols(','.join(variables))")
    expr=expr.replace('and','&') #convert to operators compatible with sympy
    expr=expr.replace('or','|')
    expr=expr.replace('not','~')
    #print expr,eval(expr)
    #eval(expr).subs({v0: True, v1: False, v2: False})
    dnf_expr=str(to_dnf(eval(expr)))
    nexpr='~('+expr+')'
    #print nexpr
    dnf_nexpr=str(to_dnf(eval(nexpr)))
    for i,node in enumerate(literals): #convert back to node names
        dnf_expr=dnf_expr.replace(variables[i],node)
        dnf_nexpr=dnf_nexpr.replace(variables[i],node)
    return dnf_expr,dnf_nexpr


#read file with expressions in any logical form and return positive/negative expressions in dnf
def file_to_dnf(filename):
    with open(filename) as f:
        lines=f.readlines()
        nodes=[line.split('*=')[0] for line in lines if line[0]!='#']
        #print len(nodes)
    expressions={}
    for line in lines:
        if line[0]!='#':
            var=line.split('*=')[0]
            #print line.strip(),var
            expr,nexpr=convert_dnf(line.strip(),nodes)
            #print '\n'.join(['Positive: '+expr,'Negative: '+nexpr,''])
            expressions[var]=(expr,nexpr)
            
    return expressions


#turn prime implicants into DNF expression
def pi_to_dnf(node): #takes a BooleanNode as input
    pi_lut=node.schemata_look_up_table(type='pi')
    inputs=node.inputs
    pos,neg="(","("
    #print pi_lut,inputs
    for k,row in enumerate(pi_lut['In:']):
        expr=""
        for i in range(len(row)):
            if row[i]=='0': expr+='~'+inputs[i]+' & '
            elif row[i]=='1': expr+=inputs[i]+' & '
        if '&' in expr: expr=expr[:-3] #cut off trailing &
        output=pi_lut['Out:'][k]
        if output==0: neg+=expr+') | ('
        else: pos+=expr+') | ('
        #print expr,pi_lut['Out:'][k]
    if '|' in pos: pos=pos[:-4]
    if '|' in neg: neg=neg[:-4]
    return pos,neg


#read file with expressions in any logical form as BooleanNetwork, reduce to PIs, and return positive/negative expressions in dnf
def file_to_pi(filename):
    
    expressions={}
    n=BooleanNetwork.from_file(filename,file_type='logical') #must convert to BooleanNetwork first to use QM algorithm
    #print f,len(n.nodes),max([len(node.inputs) for node in n.nodes]),n.name
    for node in n.nodes:
        var=node.name
        expr,nexpr=pi_to_dnf(node)
        #print var,'\n'.join(['Positive: '+expr,'Negative: '+nexpr,''])
        expressions[var]=(expr,nexpr)
    return expressions


#turn a BooleanNetwork into a general threshold network based on logical rules in disjunctive normal form (dnf)
#threshold networks are digraphs that need variable, state, type, tau, delay, label (what displays as the node label)
def create_threshold_network_dnf(n, filename, single_nodes=[],pi=False):
    """ Create a digraph compatible with threshold network manipulation
    Expects a BooleanNetwork as input, returns a threshold-like digraph where every node is a literal and/or threshold node
    
    Valid nodes include these properties: label, threshold (tau), type, time delay, variable, state
    Valid edges include these properties: type 
    
    Set single_nodes to specify nodes that should NOT have an opposing state (e.g. non-binary nodes) 
    Assumes that a filepath is passed that contains all of the logical expressions (not necessarily in dnf) """
    
    tn = nx.DiGraph()
    tn.add_nodes_from([node.name+'-1' for node in n.nodes],state=1) #add positive s-units as nodes
    tn.add_nodes_from([node.name+'-0' for node in n.nodes if node.name not in single_nodes],state=0) #add negative s-units as nodes
    #add node properties
    for name in tn.nodes():
        tn.node[name]['delay']=1
        tn.node[name]['tau']=0
        tn.node[name]['type']='variable'
        tn.node[name]['variable']=name[:-2]
        tn.node[name]['label']=name[:-2]
        
    #add threshold nodes/edges by adding one per disjunction
    #NOTE: if the node doesn't have an opposing state, ignore all LUT entries that result in 0
    if pi: expressions=file_to_pi(filename) #reduce to prime implicants first and then find dnf
    else: expressions=file_to_dnf(filename) #positive/negative expressions in dnf
    for node in n.nodes:
        #print node.name,expressions[node.name]
        for k in [1,0]: #ON/OFF state
            for i,clause in enumerate(expressions[node.name][abs(k-1)].split("|")):
                if k==0 and node.name in single_nodes: continue #ignore the opposing state
                literals=re.findall(r'~*\w+',clause.strip())
                #inputs=[lit+'-0' if '~'+lit in clause else lit+'-1' for lit in literals] #note: we assume no contradictions
                inputs=[lit[1:]+'-0' if '~' in lit else lit+'-1' for lit in literals] #note: we assume no contradictions
                not_in_singles = lambda x: x[-1]!='0' or x[:-2] not in single_nodes #filter out opposing states
                inputs = filter(not_in_singles, inputs) #solution from https://stackoverflow.com/questions/4915920/how-to-delete-an-item-in-a-list-if-it-exists
                length=len(inputs)
                if not length: continue #this row was filtered out
                output=node.name+'-'+str(k)
                name='T-'+str(i)+'_'+output
                #print clause.strip(),inputs,output,name,length
                tn.add_node(name,delay=0,tau=length,type='threshold',label=length,group=node.name)
                for inp in inputs: tn.add_edge(inp,name) #add edges to threshold unit
                tn.add_edge(name,output) #add edge from threshold unit to output s-unit
        
    return tn


#alter a Boolean network to match the original node names; takes a BooleanNetwork as input
#TODO: the modifications here are only partially complete
def alter_network(n):
    name=n.name
    nodes=n.Nnodes
    logic=n.logic
    #print name,nodes,logic
    for key in logic:
        #print key,logic[key]['name'],
        for illegal in illegals:
            #if 'a'+str(ord(illegal))+'a' in logic[key]['name']: 
            logic[key]['name']=logic[key]['name'].replace('a'+str(ord(illegal))+'a',illegal)
        #special cases:
        if name=="PC12 Cell Differentiation":
            logic[key]['name']=logic[key]['name'].replace('Ga40ai/oa41a','G(i/o)') #replace only these parentheses
        if name=="Signaling Pathway for Butanol Production in Clostridium beijerinckii NRRL B-598":
            logic[key]['name']=logic[key]['name'].replace('NADa40aPa41aH','NAD(P)H')
            logic[key]['name']=logic[key]['name'].replace('a51a-hydroxybutyryl-CoA','3-hydroxybutyryl-CoA')
        if name=="TOL Regulatory Network": #replace leading number
            logic[key]['name']=logic[key]['name'].replace('a51aMBz','3MBz')
        if name=="Signal Transduction in Fibroblasts":
            logic[key]['name']=logic[key]['name'].replace('PalphaiR','Palpha_iR') #causes issue with DNF variable replacement
            logic[key]['name']=logic[key]['name'].replace('Palpha1213R','Palpha_1213R')
            logic[key]['name']=logic[key]['name'].replace('PalphasR','Palpha_sR')
        #print logic[key]['name']
    return BooleanNetwork(name=name, logic=logic, Nnodes=nodes)


#alter a list of networks
def alter_networks(ls):
    return [alter_network(x) for x in ls]


#convert seeds to DNF ASCII versions to calculate DNF modules
def convert_seeds_ascii(seeds,name=None):
    new_seeds=set([])
    for seed in seeds:
        new_seed={s for s in seed}
        for illegal in illegals:
            new_seed={s[:-2].replace(illegal,'a'+str(ord(illegal))+'a')+s[-2:] for s in new_seed}
            
        #special cases
        if name=="PC12 Cell Differentiation":
            new_seed={s.replace('G(ia47ao)','Ga40aia47aoa41a') for s in new_seed} #replace only these parentheses
        if name=="Signaling Pathway for Butanol Production in Clostridium beijerinckii NRRL B-598":
            new_seed={s.replace('NAD(P)H','NADa40aPa41aH') for s in new_seed}
            new_seed={s.replace('3a45ahydroxybutyryla45aCoA','a51aa45ahydroxybutyryla45aCoA') for s in new_seed}
        if name=="TOL Regulatory Network": #replace leading number
            new_seed={s.replace('3MBz','a51aMBz') for s in new_seed}
        if name=="Signal Transduction in Fibroblasts":
            new_seed={s.replace('Palpha_iR','PalphaiR') for s in new_seed} #causes issue with DNF variable replacement
            new_seed={s.replace('Palpha_1213R','Palpha1213R') for s in new_seed}
            new_seed={s.replace('Palpha_sR','PalphasR') for s in new_seed}
        if name=="SKBR3 Breast Cell Line Long-term ErbB Network":
            new_seed={s.replace('ERBB2','ER_BB2') for s in new_seed} #causes issue with DNF variable replacement
        if name=="Signaling Pathway for Butanol Production in Clostridium beijerinckii NRRL B-598":
            new_seed={s.replace("acetaldehyde","acetal_dehyde") for s in new_seed} #causes issue with DNF variable replacement
        if name=="IL-6 Signalling":
            new_seed={s.replace("camk24","camk_24") for s in new_seed} #causes issue with DNF variable replacement
            new_seed={s.replace("dum_gab1_kin_or_jak1_gab1_mem_p","dum_gab1_kin_jak1_gab1memp") for s in new_seed}
            new_seed={s.replace("dum_il6rc_p_or_grb2_vav","dum_il6rc_p_grb2_vav") for s in new_seed}
            new_seed={s.replace("dum_mtorc1_or_pkcd_stat3_ta","dum_mtorc1_pkcd_stat3_ta") for s in new_seed}
        if name=="Signaling in Macrophage Activation":
            new_seed={s.replace("BAG4_TNFRSF1A","BAG4_TN_FRSF1A") for s in new_seed}
            new_seed={s.replace("TNF_BAG4_TN_FRSF1A","TN_BAG4_TN_RSF1A") for s in new_seed}
        if name=="T Cell Receptor Signaling": #causes issue with DNF variable replacement
            new_seed={s.replace("vav1","vaa118a1") for s in new_seed}
        #print seed,new_seed
        new_seeds.add(frozenset(new_seed))
    return new_seeds


#convert from DNF ascii names to standard names used in expression files
def convert_seeds_standard(modules,name=""):
    new_modules={}
    for key in modules:
        new_key={s for s in key}
        new_value={s for s in modules[key]}
        #print new_key,new_value
            
        #special cases first
        if name=="PC12 Cell Differentiation":
            new_key={s.replace('Ga40aia47aoa41a','G(ia47ao)') for s in new_key} #replace only these parentheses
            new_value={s.replace('Ga40aia47aoa41a','G(ia47ao)') for s in new_value}
        if name=="Signaling Pathway for Butanol Production in Clostridium beijerinckii NRRL B-598":
            new_key={s.replace('NADa40aPa41aH','NAD(P)H') for s in new_key}
            new_value={s.replace('NADa40aPa41aH','NAD(P)H') for s in new_value}
            new_key={s.replace('a51aa45ahydroxybutyryla45aCoA','3a45ahydroxybutyryla45aCoA') for s in new_key}
            new_value={s.replace('a51aa45ahydroxybutyryla45aCoA','3a45ahydroxybutyryla45aCoA') for s in new_value}
            new_key={s.replace("acetal_dehyde","acetaldehyde") for s in new_key} #causes issue with DNF variable replacement
            new_value={s.replace("acetal_dehyde","acetaldehyde") for s in new_value}
        if name=="TOL Regulatory Network": #replace leading number
            new_key={s.replace('a51aMBz','3MBz') for s in new_key}
            new_value={s.replace('a51aMBz','3MBz') for s in new_value}
        if name=="Signal Transduction in Fibroblasts":
            new_key={s.replace('PalphaiR','Palpha_iR') for s in new_key} #causes issue with DNF variable replacement
            new_value={s.replace('PalphaiR','Palpha_iR') for s in new_value}
            new_key={s.replace('Palpha1213R','Palpha_1213R') for s in new_key}
            new_value={s.replace('Palpha1213R','Palpha_1213R') for s in new_value}
            new_key={s.replace('PalphasR','Palpha_sR') for s in new_key}
            new_value={s.replace('PalphasR','Palpha_sR') for s in new_value}
        if name=="SKBR3 Breast Cell Line Long-term ErbB Network":
            new_key={s.replace('ER_BB2','ERBB2') for s in new_key} #causes issue with DNF variable replacement
            new_value={s.replace('ER_BB2','ERBB2') for s in new_value}
        if name=="IL-6 Signalling":
            new_key={s.replace("camk_24","camk24") for s in new_key} #causes issue with DNF variable replacement
            new_value={s.replace("camk_24","camk24") for s in new_value}
            new_key={s.replace("dum_gab1_kin_jak1_gab1memp","dum_gab1_kin_or_jak1_gab1_mem_p") for s in new_key}
            new_value={s.replace("dum_gab1_kin_jak1_gab1memp","dum_gab1_kin_or_jak1_gab1_mem_p") for s in new_value}
            new_key={s.replace("dum_il6rc_p_grb2_vav","dum_il6rc_p_or_grb2_vav") for s in new_key}
            new_value={s.replace("dum_il6rc_p_grb2_vav","dum_il6rc_p_or_grb2_vav") for s in new_value}
            new_key={s.replace("dum_mtorc1_pkcd_stat3_ta","dum_mtorc1_or_pkcd_stat3_ta") for s in new_key}
            new_value={s.replace("dum_mtorc1_pkcd_stat3_ta","dum_mtorc1_or_pkcd_stat3_ta") for s in new_value}
        if name=="Signaling in Macrophage Activation":
            new_key={s.replace("TN_BAG4_TN_RSF1A","TNF_BAG4_TN_FRSF1A") for s in new_key}
            new_value={s.replace("TN_BAG4_TN_RSF1A","TNF_BAG4_TN_FRSF1A") for s in new_value}
            new_key={s.replace("BAG4_TN_FRSF1A","BAG4_TNFRSF1A") for s in new_key}
            new_value={s.replace("BAG4_TN_FRSF1A","BAG4_TNFRSF1A") for s in new_value}            
        if name=="T Cell Receptor Signaling": #causes issue with DNF variable replacement 'TNF_BAG4_TNFRSF1A-1'
            new_key={s.replace("vaa118a1","vav1") for s in new_key}
            new_value={s.replace("vaa118a1","vav1") for s in new_value}
            
        for illegal in illegals:
            new_key={s[:-2].replace('a'+str(ord(illegal))+'a',illegal)+s[-2:] for s in new_key}
            new_value={s[:-2].replace('a'+str(ord(illegal))+'a',illegal)+s[-2:] for s in new_value}

        #if key!=new_key: print key,new_key,'\n'
        #if new_value!=modules[key]: print modules[key],new_value,'\n'
        new_modules[frozenset(new_key)]=new_value
    return new_modules

