import numpy as np, pandas as pd
df_file = 'space_groups.pkl'
df = pd.read_pickle(df_file)

def apply_sym_cif(cif):
    with open(cif,'r') as f:lines=f.readlines()

    #### ops
    loops = [i for i,l in enumerate(lines) if 'loop_' in l]
    li = [i for i,l in enumerate(lines) if '_symmetry_equiv_pos_as_xyz' in l][0]
    ops = [l.split("'")[1] for l in lines[li+1:loops[-1]]]
    #### positions
    li = [i for i,l in enumerate(lines) if '_atom_site_occupancy' in l][0]
    lf = [i for i,l in enumerate(lines) if '#End' in l][0]
    xyz =np.array([np.array(l.split('  ')[1:4],dtype=float) for l in lines[li+1:lf]])

    X = []
    for x0 in xyz:
        x,y,z = x0
        X += [np.array(eval(op)) for op in ops]    #;print(X)
    X = np.array(X)
    X[X<0]+=1
    X[X>1]-=1
    X = np.unique(np.array(X),axis=0)
    return X

def apply_sym(xyz,no):
    '''
    xyz : asymmetric unit cell
    no : space group
    '''

    c=df.loc[no]
    spg = c['space group']                              #;print(spg)
    ops = [l for l in c.ops if l.strip() != '']         #;print(ops)
    X = np.array([np.array(eval(op)) for op in ops])    #;print(X)
    X[X<0]+=1
    X[X>1]-=1
    dx[i] = X
    return X

if __name__=="__main__":
    cif='/home/tarik/Documents/git/github/Felix/samples/GaAs_short/felix.cif'
    with open(cif,'r') as f:lines=[l.strip() for l in f.readlines()]
    #### ops
    loops = [i for i,l in enumerate(lines) if 'loop_' in l]
    li = [i for i,l in enumerate(lines) if '_symmetry_equiv_pos_as_xyz' in l][0]
    ops = [l.split("'")[1] for l in lines[li+1:loops[-1]]]
    #### positions
    li = [i for i,l in enumerate(lines) if '_atom_site_occupancy' in l][0]
    lf = [i for i,l in enumerate(lines) if '#End' in l][0]
    xyz =np.array([np.array(l.split('  ')[1:4],dtype=float) for l in lines[li+1:lf]])

    # apply_sym_cif(cif)
