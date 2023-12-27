from IPython.display import display
import fct

(t, date, P_el_total) = fct.GetProfile('groupe_5.csv', Shouse=1.)  
LDC_elec = fct.LoadDurationCurve(t, date, P_el_total, color='blue', name='Electricity Demand', steps=250)
display(LDC_elec)



(t, date, tau_pv) = fct.GetPVProduction('SolarData.csv', Pnom = 3000) #Pnom = puissance crête de l'installation PV en kWc
LDC_PV = fct.LoadDurationCurve(t, date, tau_pv, color='green', name='PV production', steps=250)
display(LDC_PV)

fct.CalculFacture(P_el_total, tau_pv, buyCost=150, sellCost=50) #Prix d'achat et de vente en €/MWh

