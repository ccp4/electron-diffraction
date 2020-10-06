from crystals import Crystal,Atom,symmetry_expansion,write_xyz
# import rotating_crystal as rcc; imp.reload(rcc)

cif='225.cif'  #201(SC), 212,225,229

lattice_vectors = np.diag([2,2,2])
cCIF = Crystal.from_cif('../dat/Carbone_SC/'+cif)
# cCIF = Crystal.from_database('Fe')
syms = cCIF.symmetry_operations()
atoms = [Atom(z, c) for z,c in zip(['C','C','C'], [[0.5,0.5,0],[0,0.5,0.5],[0.5,0.,0.5]]) ]
# atoms = [Atom('C', [0.5,0.5,1])]
# atoms = [Atom('C', [1.,1.,1.])]
syms4 = [ np.hstack([np.vstack([s[0],[0,0,0]]), np.vstack([s[1][:,None],1]) ]) for s in syms]
unitcell = symmetry_expansion(atoms,syms4)
C = Crystal(unitcell, lattice_vectors)
print('international_number:', cCIF.symmetry()['international_number'])
print('international_number:', C.symmetry()['international_number'])
print(C)
# write_xyz(C,'test.xyz')
