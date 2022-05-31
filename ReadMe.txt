University of Bilkent
ME430: HEAT EXCHANGERS § DESIGN
SHELL AND TUBE HEAT EXCHANGER DESIGN FOR OLYMPIC POOLS
Instructor: Dr. Barbaros ÇETIN
Prepared by Yunus Ozkan and Bahadir Bilir
Please read the instructions before use it.
For any questions, please raise an issue on git or reach us by using e-mail, yunuszkn01@gmail.com

#All formulas taken from Heat Exchangers Selection, Rating and Thermal Design by Sadik Kakaç ,Hongtan Liu, Anchasa Pramuanjaroekenji

According to TEMA:

0.0666 < (Shell diameter/Tube length) < 0.2 
1.25 < Pitch/Outer diameter of tube < 1.5 

POOL HX DESIGNER

You must read the followings before using the program:

Tube sides are in HX: Cold Water

Shell sides are in HX: Hot Water

Fluids in flow must be in single phase. Condensation or boiling cannot be tolerated, check the values beforehand.

We are designing a pool hx, because of that, we get water inlet in certain small degree range and we want water outlet in certain small degree range.
Since we cannot change the outlet water temperature and we are expecting to get a certain outlet temperature. 
What we can do for get correct desing is that iterating heat transfer coefficient value values until both preliminary design and main part gives very close results. 
Briefly, we iterate the whole design since we cannot change the outlet temperature.


We designed main part variables as flexible because industry products are not exactly the same with offered results. 
One should keep the values and choices of design and the main part as much close as possible to each other if wants to get accurate results.
Every small change will create a deviation in the outlet temperature of the pool. Changes should be kept small as possible. 

Run the preliminary design step by step from left to right, top to bottom.
You will encounter errors and wrorng results if you change design variables after calculation in the preliminary design part.(ex: Changing temperature values after calculating fluid properties)
If you want to change preliminary design and fluid properties, we suggest you to re-start the program.
Iteration and main part can be changed freely without facing any error as long as values are logical.

This software is not the final product. However, when used properly, it can be used in pool heat exhanger designs.

Maximum working pressure on both sides are 6 bar.



Best regards,

Yunus Özkan
Bahadir Bilir
