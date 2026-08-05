This is Polynomial function for CEE lambda using Adiabatic Mass loss models. 
The key module is formula_test.py. Try to use test_poly.py to read the formulae in folder lambda_CE_formula.

there will be multiple functions in formula_test.py, determined by these parameters:
Mi: initial mass in solar scale (Msun) before CEE
Ri:  initial radius in solar scale (Rsun) before CEE
eta: envelope mass residual in relative solar scale (Menv/Mi)
folder_path: the path to .csv files under this folder

try to use ''from formula_test import lambda_ce_3d_fgb'' to read 3d fitting formula
or ''from formula_test import lambda_ce_fgb_eta_1e3'' to read 2d fitting formula

After the right input of {Mi,Ri,eta,folder_path} in the 3d function, or {Mi,Ri} in the 2d function
every function shall return a single value of lambda_ad, which will be close to lambda_core*F1*F2 in RBF tables.