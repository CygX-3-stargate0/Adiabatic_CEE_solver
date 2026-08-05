#!/usr/bin/env python
# coding: utf-8

# In[9]:


#3D
from formula_test import lambda_ce_3d_fgb
from formula_test import lambda_ce_3d_agb
from formula_test import lambda_ce_3d_hg
from formula_test import lambda_ce_3d_massive
#2D
from formula_test import lambda_ce_fgb_eta_1e3
from formula_test import lambda_ce_fgb_eta_1e2
from formula_test import lambda_ce_agb_eta_1e3
from formula_test import lambda_ce_agb_eta_1e2
from formula_test import lambda_ce_hg_eta_1e3
from formula_test import lambda_ce_hg_eta_1e2
from formula_test import lambda_ce_massive_eta_1e3
from formula_test import lambda_ce_massive_eta_1e2

#3d poly formula test

path='lambda_CE_formula'
Mi=1.5
Ri=60
eta=0.0056
lam_HG_3d=lambda_ce_3d_fgb(Mi,Ri,eta,path)
print(lam_HG_3d)

#2d poly formula test

Mi=1.5
Ri=60
eta=1e-2
lam_HG_2d=lambda_ce_fgb_eta_1e2(Mi,Ri)
print(lam_HG_2d)


# In[ ]:




