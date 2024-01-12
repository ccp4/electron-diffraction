from utils import*                          ;imp.reload(dsp)
from EDutils import dials_utils as dials    ;imp.reload(dials)
from EDutils import utilities as ut         ;imp.reload(ut)



d = dials.Dials('dat/dials')

# print('R1:',d.R(1))
# print('R2:',d.R(2))
print(d.get_UB(0.5))
print(d.UBs[0])
# print(UB)
