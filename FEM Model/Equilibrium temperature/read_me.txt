- CalcEqTemp_FEM_Geo is the main file of this folder
it calculates the equilibrium temperatures according to the desired conditions
and calculates the residuals as well 

- fraction_power and gebhart_3d files calculate the beam heat flux inside the cavity
and the radiation loss factors

- plot_geometry plots the geometry in a figure with the desired parameter
such as heat flux and temperature
the legend needs to be modified inside accordingly

- Teq_geometry_conditions set the governing equations for each node
rad - inner radiation only
radcond - inner radiation loss and conduction between nodes
to add outer radiation to one of them in the governinng equation
for each node it needs to be added ( - emi*Area_node*stefan_constant*(T_node^4-T_amb^4))
