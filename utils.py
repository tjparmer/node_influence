#Utility functions for influence calculations
#Thomas Parmer, 2022

from itertools import combinations
import random
from brute_force_computations import reverse_sunit_map
from modules import sunit_to_var

def compute_jaccard(set1,set2):
    """ return the jaccard index of two sets: intersection / union """
    return float(len(set1.intersection(set2))) / len(set1.union(set2))


#list method that ensures that s-unit combinations are in the correct order
#order is based on node number and state, starting with 0
def to_list(seed,smap,translator=None):
    """ determine the sunits involved in the seed, based on the sunit_map 
    this make take either a tuple of s-units, or a string of s-units (requires the translator argument)
    the smap maps individual node numbers to names, the translator maps seed strings to node numbers """
    
    ls = []
    if isinstance(seed,str):
        seed = translator[seed]
        
    for node in seed:
        ls.append(smap[node])
    return ls


#define s-units and modules; map any node set to numbers starting with 0
def get_sunits(N):
    sunits,sunit_map=set([]),{}
    num=0
    for node in N.nodes:
        for state in ['0','1']:
            sunits.add(num)
            sunit_map[num]=str(node.name)+'-'+state
            num+=1
        
    return sunits,sunit_map


#reduce seeds based on length (range inclusive) and contradiction
def reduce_seeds(seeds,sunit_map,translator,lrange=None,length=1):
    new_seeds=set([])
    if not lrange: #NOTE: if lrange is set, length argument is ignored
        lrange=[length,length]
    for seed in seeds: 
        sunits=to_list(str(seed),sunit_map,translator)
        if len(sunits)<lrange[0] or len(sunits)>lrange[1]: #incorrect length
            continue
        if len({sunit[:-2] for sunit in sunits})!=len(sunits): #there is a contradiction
            continue
        new_seeds.add(seed)
        
    return new_seeds


#create translator function if needing a translator but not wanting to run the actual mf-approximation
def create_translator(N,s=1,sunits=None,sunit_map=None,translator={},seeds=None,samples=None):
    """ find all pathway modules seeds for a given network N and seed size s, with the given translator,
    can iteratively add to translator for different s values; returns updated translator
    set seeds to a list of which seeds you want to find modules for (by s-unit number)
    or set samples to sample the space of possible seeds """
    
    #define s-units if they are not defined outside the function
    if not sunits or not sunit_map:
        sunits,sunit_map = get_sunits(N)
    #sunits={i for i in range(len(lsunits))}
    #define seeds
    if not seeds:
        if samples: #sampling with replacement to avoid memory issues with high combinations
            seeds=[]
            while len(seeds)<samples:
                ls=tuple(sorted([random.choice(list(sunits)) for j in range(s)])) #use a tuple to be consistent
                if len(set(ls))==s:
                    seeds.append(ls)
        else:
            seeds=list(combinations(sunits,s)) #[['en-1']] or list(combinations(sunits,s)) for example
    #print 'seeds:',len(seeds)
    for seed in seeds:
        if isinstance(seed,str): seed=eval(seed) #NOTE: eval in case these are strings
        translator[str(seed)] = seed #map between the string and the actual seed numbers
        
    return translator


##ADDED FOR 2023 PAPER ON TARGET CONTROL

#filter a seed for contradictions
def filter_seed(seed):   
    vals=[node[:-2] for node in seed]
    return len(set(vals))==len(vals) #check to see if we have multiple nodes sharing the same variable


#get all possible sunit combinations
def enumerate_sunits(nodes,num=1,samples=None):
    
    if samples: #sampling to avoid memory issues with high combinations
        seeds=[[random.choice(nodes) for j in range(num)] for i in range(samples)]
    else: #get all combinations
        seeds=combinations(nodes,num)
        
    filtered_seeds=set([])
    #filter for contradictions
    return {frozenset(seed) for seed in seeds if filter_seed(seed)}


#get all possible sunit combinations using a sunit map
def enumerate_sunit_map(sunit_map,num=1,samples=None):
    
    rmap=reverse_sunit_map(sunit_map)
    seeds=enumerate_sunits(rmap.keys(),num=num,samples=samples)
    return {frozenset(to_list(seed,rmap)) for seed in seeds}


#transform a frozenset seed into a tuple
def set_to_tuple(seed,sunit_map):
    rmap=reverse_sunit_map(sunit_map)
    return tuple([rmap[s] for s in seed])


#return sunits from a BooleanNetwork
def get_sunit_names(n):
    sunits,sunit_map=get_sunits(n) #gets sunit numbers
    return reverse_sunit_map(sunit_map).keys()


#retrieve samples from the given nodes, including no contradictions; num is how large each sample should be
#NOTE: this takes as input a list of sunits; default is selection without replacement
#TODO: add a timeout; if given more samples than possible seeds, it will run forever!
def get_sample(sunits,num=1,samples=1000,replacement=False):
    if replacement: seeds=[] #may select the same seed multiple times
    else: seeds=set([]) #will select each seed at most once
    while len(seeds)<samples:
        seed=[random.choice(sunits) for j in range(num)]
        vals=[sunit_to_var(s) for s in seed]
        #print seed,vals
        if len(set(vals))==len(vals): #no contradictions
            if replacement: seeds.append(seed) #this is a list of lists
            else: seeds.add(frozenset(seed)) #this is a set of sets
    return seeds

