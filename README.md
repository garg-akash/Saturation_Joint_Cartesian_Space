# Saturation_Joint_Cartesian_Space

**Implementing a simple Cartesian tracking controller**
<br />
Result of cartesian tracking control loop. (Not tunned(Left); Tracking output with duration T = 10(Center); Tracking output with duration T = 20(Right))
![alt text](images/cartesian_tracking_circular.png)
<br />
**Vrep result** of the Original algorithm (joint limit constraints having higher priority than Cartesian limit constraints)vs Proposed algorithm (no priority ordering of joint limit constraints and Cartesian limit constraints; all constraints assigned same priority) 
<br />
* Original Algorithm
![alt text](images/Orig.png)
<br />
* Proposed Algorithm
![alt text](images/Prop.png)
<br />
**Graphical result** of the Original algorithm vs Proposed algorithm
<br />
* Original Algorithm
![alt text](images/Orig_graph.png)
<br />
* Proposed Algorithm
![alt text](images/Prop_graph.png)