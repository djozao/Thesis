- Main.m is the main file where the transient simulation runs
it is commented with respect to each section
open this first and the other files purpose is easier to understand
- cone/cylinder/radiativeloss/outer_convloss/outer_radloss and similarly named files are related to 
geometrical definition and their respective losses 
- air, ammonia, hydrogen, nitrogen and water have the fluid properties and are used for the propellant properties as well as convection losses
- "duct" named files are used to calculate the heat transfre to the propellant, the main ones are named spiral_duct and linear_duct
- T_insulation_calculation is used when insulative layer is considered and calculates the outer surface temperature
-power absorbed function calculates the absorbed power considering the geometry and surface properties of the cavity
- nistdata is used to extract fluid properties