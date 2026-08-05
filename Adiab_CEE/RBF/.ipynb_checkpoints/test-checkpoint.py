#!/usr/bin/env python
# coding: utf-8

# In[1]:


import Adiab_CEE


# In[2]:


from Adiab_CEE import RecCEE_solver

stage='RGB'#HG,RGB,AGB,massive
core='core1'#core1:XH=0.01, core2:XH=0.1
use_grid='abs'#abs:Absolute(Menv=eta*1), rel: Relative(Menv=eta*Mi)

Mi=1.4
Ri=120
Menv=1e-1

#-------------------------RecCEE_solver------------------------#

result=RecCEE_solver(Mi,Ri,Menv,stage,core,use_grid)

lam0=result[0]
F1=result[1]
F2=result[2]

lam_ad=lam0*F1*F2

print(lam0,F1,F2)
print(lam_ad)


# In[4]:


from Adiab_CEE import Rfad_solver

stage='RGB'#HG,RGB,AGB,massive
core='core1'#core1:XH=0.01, core2:XH=0.1

Mi=1.4
Ri=120
Menv=1e-1 #Menv<0.1*Mi

#-------------------------Rfad_solver------------------------#

Rf=Rfad_solver(Mi,Ri,Menv,stage,core)

print(Rf)


# In[ ]:




