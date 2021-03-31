import pickle,matplotlib
import numpy as np, pandas as pd
import scipy.fftpack as fft
from scipy.integrate import nquad,trapz,quad
import utils.displayStandards as dsp
import utils.physicsConstants as cst
import utils.glob_colors as colors

class Multislice:
  def __init__(self,
          pattern,ax,bz,
          keV=200,Nx=1,tilt=0,
          dz=1,nz=1,copt=1,eps=1,sg=-1,
          TDS=False,nTDS=8,ndeg=2**10,wobble=0.05,
          iZs=1,opts='q',iZv=1,v=1,ppopt=''):
      self.version = 1.0
      self.pattern = pattern
      self.keV     = keV
      self.Nx      = Nx
      self.dz      = dz
      self.ax      = ax
      self.bz      = bz
      #energy dependent
      self.lam  = cst.keV2lam(keV)
      self.sig  = cst.keV2sigma(keV)
      self.k0   = 1/self.lam
      self.eps  = eps
      self.tilt = tilt
      self.copt = copt
      #TDS
      self.TDS  = TDS
      self.nTDS = nTDS
      self.ndeg = ndeg
      self.wobble = self._wobble(wobble)

      #Computations
      if v:print(colors.red+'\t\t 2D multislice simulation '+colors.black)

      # Thermal diffuse scattering
      if self.TDS:
          self._TDS(eps,copt,sg,v)
      else:
          # standard
          self._set_transmission_function(eps,v)
          self._set_propagator(sg,copt)
          self.set_Psi0()
          if nz:self.propagate(nz,iZs,iZv,opts,v)

      display(ppopt)
