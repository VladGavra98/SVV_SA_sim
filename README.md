# SVV_SA_sim
The simulation tool for the Structural Assignment, SVV 2020

When you write a function, PLEASE  copy the template that I used to write my first functions.

        TO DO:
Report: 
-For the structure, more elaboration on the problem in the introduction or a seperate chapter as problem analysis instead of putting everything in the numerical model part. 
- I would also adjust the position of your verification model, so that it would be analysed before your verification process. **Done**
- The structual model description needs to be improved.
- Mistake in the load case. Actuator I and II should have the same positions. **Done**
- Make FBD more clear (redo the 3D one).
- Clearly state where you put the flexural axis and how to translate the point loads and the resultant aerodynamic load.
- Clearer description about beam section.
- A clarification on how you use this VM definition to superimpose the effect of bending, shear and torsion(probably) is expected. And do take care of your reference system trasition. **Done**
- Flowchart is missing the corresponding variables on your links which can be very useful for verification.
- Change Verification Model. **Done**


So far the contained modules are:

main_sim.py :
  - shearflowCalc()  -- shear flow due to Sy and Sz  [consider moving it to a new file]
  - inertiaCalc()
  - posStCalc()      -- stringer postions
  - Aicraft class    -- constructor defines geomtry of the specified aircraft (everything in m)
  - drawSection()    -- plots the yz-section 
  - centroidCalc()   -- return the Zc 
  -circCalc()        -- length of the circumference (aslo used in shear flow)
integration.py
interpolation.py
