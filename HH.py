# %%
import numpy as np
from itertools import combinations

# %%
#Swap two objects
def swap(a,b):
    return [b.copy(),a.copy()]

# %%
#Binary reduced-row echelon form
def brref(M):
    B=M.copy()
    result={}
    rows,cols=M.shape
    pivots,non_pivots=[],[]
    if not rows*cols:
        result['matrix']=B
        result['pivots']=[]
        result['nonpivots']=[t for t in range(cols)]
        return result
    
    piv_row=0
    for j in range(cols):
        found=False
        for i in range (piv_row,rows):
            if B[i,j]:
                if i!=piv_row:
                    B[i],B[piv_row]=swap(B[i],B[piv_row])
                found=True
                pivots.append(j)
                break
        if not found:
            non_pivots.append(j)
            continue
        for i in range(0,rows):
            if B[i,j] and i!=piv_row:
                B[i]=B[i]^B[piv_row]  
        piv_row+=1

    result['matrix']=B
    result['pivots']=pivots
    result['nonpivots']=non_pivots
    return result

# %%
#Find kernel given a matrix in rref
def ker(rref):
    kernel={}
    rows,columns=rref['matrix'].shape
    for x in rref['nonpivots']:
        vector=np.zeros(columns, np.int8)
        vector[x]=1
        kernel[x]=vector
    l=len(rref['pivots'])
    for i in range (l):
        for x in rref['nonpivots']:
            if rref['matrix'][i,x]==1:
                kernel[x][rref['pivots'][i]]=1
    kernel=[kernel[x] for x in kernel]
    return (kernel)

# %%
#Find Image given a matrix and its rref
def im(M,rref):
    image=[]
    for column in rref['pivots']:
        image.append(M[:,column])
    return image

# %%
#Read the maximal simplices and put them in a list
def read():
    simps=[]
    while True:
        string=input()
        if ';' in string:
            break
        simplex= tuple(map(int, string.split()))
        simps.append(simplex)
    return simps

# %%
#Generate simplicial complex from maximal simplices
def gen_scx(MaxSimp,shift=1):
    scx = {}  # Dictionary to store subtuples by size
    if shift:
        scx[-2*shift]=[]
    else:
        scx[0]=[()]
    for simplex in MaxSimp:
        # Iterate over all possible subtuple lengths
        for size in range(1-2*shift,len(simplex) + 1):
            # Generate all possible subtuples of the current size
            subsimplex = set(combinations(simplex, size+shift))

            # Add subtuples to the dictionary, grouped by size
            if size not in scx:
                scx[size] = set()
            scx[size].update(subsimplex)
    for x in scx:
        scx[x]=list(scx[x])

    return scx


# %%
#Simplicial Homology Matrix
def SH_Matrix(simps, dim):
    n=len(simps[dim-1])
    m=len(simps[dim])
    M=np.zeros((n,m),dtype=np.int8)
    for i in range (n):
        for j in range (m):
            if set(simps[dim-1][i]).issubset(set(simps[dim][j])):
                M[i,j]=1
    return M

# %%
#Find all RSH differentials
def delta(simps):
    maxdim=len(simps)-2
    matrices={}
    for dim in range (-1,maxdim):
        matrices[dim]=SH_Matrix(simps,dim)
    return matrices

# %%
#RSH bases computing
def RSH(simp_cx):
    maxdim=len(simp_cx)-2
    kernel={}
    image={}
    gen={}
    rref={}
    delt=delta(simp_cx)
    res={}
    for dim in range(-1,maxdim):
        rref[dim]=brref(delt[dim])

    for dim in range(-1,maxdim):
        image[dim]=im(delt[dim],rref[dim])
        kernel[dim]=ker(rref[dim])

    for dim in range(-1,maxdim-1):
        if (kernel[dim]==[]):
            continue
        gen[dim]=[]
        vectors=[x for x in image[dim+1]]
        l=len(vectors)
        for x in kernel[dim]:
            vectors.append(x)
        M=np.column_stack(vectors)
        ref=brref(M)
        for x in ref['pivots']:
            if x>=l:
                gen[dim].append(M[:,x])
        rank=len(gen[dim])
        if not rank:
            continue
        res[dim]={}
        res[dim]['rank']=rank
        res[dim]['basis']=gen[dim]
        res[dim]['modded']=image[dim+1]
    return res

# %%
#Generate full subcomplex
def gen_fscx(simp_cx,subset):
    subs=set(subset)
    new_scx={}
    used_sets=[]
    for i in range(-2,len(subset)+2):
        new_scx[i]=[]

    for l in simp_cx:
        for simp in simp_cx[l]:
            s=tuple(sorted(set(simp)&subs))
            if s in used_sets:
                continue
            else:
                new_scx[len(s)-1].append(s)
                used_sets.append(s)

    return new_scx

# %%
#Hochster Decomposition
def hochster(simp_cx,m):
    triang=[tuple(range(1,m+1))]
    all_subs=gen_scx(triang,0)
    full_scx={}
    max_deg=0
    #This finds the hochster decomposition and sorts it by size
    for size in all_subs:
        list_of_fullsubs=[]
        for subset in all_subs[size]:
            result={}
            result['subset']=subset
            result['complex']=gen_fscx(simp_cx,subset)
  

            i=len(result['complex'])-2
            #This loop ensures the complexes are only of the needed size
            for i in reversed(range(-1,i)):
                if (result['complex'][i]==[] and result['complex'][i-1]==[]):
                    del result['complex'][i]
                else:
                    break
            result['homology']=RSH(result['complex'])
            #This weeds out things of homology 0
            for x in result['homology']:
                max_deg=max(max_deg,x)
            if result['homology']=={}:
                continue
            list_of_fullsubs.append(result)
        full_scx[size]=list_of_fullsubs
    #This re-sorts it by homological degree and then size for easier use when computing HH
    hoch={}
    for i in range(-1,max_deg+1):
        hoch[i]={}
    for size in full_scx:
        for component in full_scx[size]:
            for degree in component['homology']:
                result={}
                result['subset']=component['subset']
                result['complex']=component['complex']
                result['homology']=component['homology'][degree]
                if size not in hoch[degree]:
                    hoch[degree][size]=[]
                hoch[degree][size].append(result)
    #This adds empty complexes at the endpoints of each homological degree to simplify computing HH
    for degree in hoch:
        if hoch[degree]!={}:
            hoch[degree][sorted(hoch[degree])[0]-1]=[]
            hoch[degree][sorted(hoch[degree])[len(hoch[degree])-1]+1]=[]

    #For each subset, it assigns a starting index number to get a numbered basis later
    for degree in hoch:
        for size in hoch[degree]:
            k=0
            for subset in hoch[degree][size]:
                subset['start_index']=k
                k+=subset['homology']['rank']
                
    return hoch

# %%
#Finds a presentation of a a vector in terms of the bases for homology
def image(vector,basis,modulo):
    vectors=[]
    for v in basis:
        vectors.append(v)
    for v in modulo:
        vectors.append(v)
    vectors.append(vector)

    rref=brref(np.column_stack(vectors))
    
    b=len(basis)
    solution=np.zeros(b)
    k=0
    for index in rref['pivots']:
        if index>b-1:
            break
        solution[index]=rref['matrix'][k,-1]
        k+=1
    return solution


# %%
#Computes the differential in CH
def HH_Dif(CH, l, k):
    hom_deg=k-1
    #figuring out the shape
    cols=0
    for subs in CH[hom_deg][l]:
        cols+=subs['homology']['rank']
    rows=0
    for subs in CH[hom_deg][l+1]:
        rows+=subs['homology']['rank']

    D=np.zeros([rows,cols],np.int8)
    init_column=0

    for small_set in CH[hom_deg][l]:
        for big_set in CH[hom_deg][l+1]:
            if (set(small_set['subset']).issubset(set(big_set['subset']))):
                shift=0
                for vector in small_set['homology']['basis']:
                    #Writes the vector when included in the larger simplicial complex
                    included_vector=np.zeros(len(big_set['complex'][hom_deg]),np.int8)
                    for i in range(len(vector)):
                        if vector[i]:
                            simplex=small_set['complex'][hom_deg][i]
                            index_in_big_set=big_set['complex'][hom_deg].index(simplex)
                            included_vector[index_in_big_set]=1

                    #Finds its representation in homology
                    sol=image(included_vector,big_set['homology']['basis'],big_set['homology']['modded'])

                    #Writes that representation in the matrix
                    a=big_set['start_index']
                    for i in range(len(sol)):
                        D[a+i,init_column+shift]=sol[i]
                    shift+=1
        init_column+=small_set['homology']['rank']
    return D
    

# %%
#Computes Rank of HH at a given homological degree given the hochster decomposition, it's indexed by l
def computehh(CH,hom_deg):
    k=hom_deg+1
    sizes=sorted(CH[hom_deg])
    if sizes==[]:
        return {}
    sizes.pop(-1)
    diff={}
    kernel={}
    image={}
    for l in sizes:
        diff[l]=HH_Dif(CH,l,k)
        rref=brref(diff[l])
        kernel[l]=ker(rref)
        image[l]=im(diff[l],rref)
    
    sizes.pop(0)

    rank={}
    for l in sizes:
        rank[l]=len(kernel[l])-len(image[l-1])

    return rank


#%%
#Reads the maximal simplices either from the terminal or a textfile
def read_simps(name,option):
    if option==1:
        file=open(name,"r")
        simps=[]
        while True:
            string=file.readline()
            if not string:
                break
            simplex= tuple(map(int, string.split()))
            simps.append(simplex)
        return simps
    else:
        print("Input each maximal simplex, to finish press enter again")
        simps=[]
        while True:
            string=input()
            if not string:
                break
            simplex= tuple(map(int, string.split()))
            simps.append(simplex)
        return simps

#%%
#Reading %%
option=int(input("Press 0 if you want to input the complex directly or press 1 to read from the \"input\" file\n"))
MaxSimp=read_simps("Input.txt",option)
#MaxSimp=[(1,2,3),(1,2,4),(1,3,4),(2,3,4)]
vertex_set=set()
for simplex in MaxSimp:
    vertex_set=vertex_set|set(simplex)
m=len(vertex_set)
simp_cx=gen_scx(MaxSimp)

#%%
#Printing %%
HD=hochster(simp_cx,m)
print("\t\t\t (-k,2l)\tRank\n")
for hom_deg in HD:
    print("Homological degree ",hom_deg)
    k=hom_deg+1
    ranks=computehh(HD,hom_deg)
    for l in ranks:
        if ranks[l]:
            print("\t\t\t(",k-l,",",2*l,")\t",ranks[l])
input("\n\nPress any key to exit")