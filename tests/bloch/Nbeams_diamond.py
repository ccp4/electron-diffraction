from utils import*                  ;imp.reload(dsp)
# from blochwave import bloch         ;imp.reload(bloch)
from blochwave import bloch_pp as bl;imp.reload(bl)
from EDutils import utilities as ut ;imp.reload(ut)
plt.close('all')

path = 'dat/bloch'

simus = ['4+2w']
# simus = ['2','3','3+1w','4','4+2w']
# bloch_rock[opts] : s(solve), t(set_thicks), I(Iz), R(rocking curve) X(Xi) S(Setup) W(Sw vs theta)
opts = 'W'
opt = 'p'
bloch_args = {'cif_file':'diamond','Nmax':8,'Smax':0.01}



# exact 2 beam orientation
if '2' in simus:
    u2 = np.array([0.13846063, 0.01406432, 0.99026807])
    rock = bl.bloch_rock(tag='diamond_2beam',u=u2,u1=[0,1],
        omega=np.linspace(-0.3,0.3,128),bloch_args = bloch_args,
        thicks=(0,800,400),ts0=0,
        zs = np.arange(168,800,168), hkls=[[[0,4,0]]],
        path=path,opts=opts,opt=opt)


# exact 3 beam orientation
if '3' in simus:
    u3 = np.array([0.23274900552973615, 0.022797452273781056, 0.9722696007768337])
    bloch_args['Smax'] = 0.008
    rock = bl.bloch_rock(tag='diamond_3beam',u=u3,u1=[0,1],
        omega=np.linspace(-0.06,0.06,128),bloch_args = bloch_args,
        thicks=(0,3000,500),ts0=0,
        zs = np.arange(400,3000,600), hkls=[[[-6,-8,2]], [[5,-3,-1]] ],
        path=path,opts=opts,opt='ps')

# exact 3+1w beam orientation
if '3+1w' in simus:
    u3 = np.array([0.23274900552973615, 0.022797452273781056, 0.9722696007768337])
    bloch_args['Smax'] = 0.01
    rock = bl.bloch_rock(tag='diamond_3_1beam',u=u3,u1=[1,0],
        omega=np.linspace(-0.06,0.06,128),bloch_args = bloch_args,
        thicks=(0,3000,500),ts0=0,strong={'I':20},
        zs = np.arange(400,3000,600), hkls=[[[-6,-8,2]], [[5,-3,-1]] ],
        path=path,opts=opts,opt=opt)

# 4 beam configurations
if '4' in simus:
    u4 = np.array([0.0702957163976856, -0.00026010222764890744, 0.9975261623651619])
    bloch_args['Smax'] = 0.01
    rock = bl.bloch_rock(tag='diamond_4beam',u=u4,u1=[0,1],
        omega=np.linspace(-0.06,0.09,128),bloch_args = bloch_args,
        thicks=(0,1200,500),ts0=0.015,zs = np.arange(200,800,150),
        hkls=[ [[4,-8,0],[4,8,0]], [[2,-6,0],[2,6,0]] ],#cond='Sw>1e3',
        path=path,opts=opts,opt=opt)


# 4+2w beam configurations
if '4+2w' in simus:
    u4 = np.array([0.0702957163976856, -0.00026010222764890744, 0.9975261623651619])
    bloch_args['Smax'] = 0.025
    rock = bl.bloch_rock(tag='diamond_4_2beam',u=u4,u1=[0,1],
        omega=np.linspace(-0.06,0.09,128),bloch_args = bloch_args,
        thicks=(0,1200,500),ts0=0.015,zs = np.arange(200,800,150),
        hkls = [ [[4,-8,0],[4,8,0]], [[2,-6,0],[2,6,0]] ],cond='(Vga>1e-3)',
        path=path,opts=opts,opt=opt)
