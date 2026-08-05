This is RBF interpolation table for CEE lambda and final radius after CE. 
The key module is Adiab_CEE.py. Try to use test.py to read the module and the .npy files.

there will be two solver in Adiab_CEE.py, determined by these parameters:
Mi: initial mass in solar scale (Msun) before CEE
Ri:  initial radius in solar scale (Rsun) before CEE
Menv: envelope mass residual in solar scale (Msun)

stage: 'HG','RGB','AGB','massive', choose one
core: 'core1' or 'core2'. where core1(XH=0.01), core2(XH=0.1) is two type of core boundaries.
use_grid: 'abs' or 'rel'. abs is for Absolute scheme(Menv=eta*1), rel is for Relative scheme(Menv=eta*Mi)

try to use ''from Adiab_CEE import RecCEE_solver'' or ''from Adiab_CEE import Rfad_solver''
After the right input of {Mi,Ri,Menv,stage,core,use_grid} in the function, 
the ''RecCEE_solver'' shall return a list of [lambda_core, F1,F2], where F1,F2 is the correction of CEE envelope residual and Adiabatic expansion.
the ''Rfad_solver'' shall just return Rf in solar scale (Rsun).