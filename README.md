# Saturation_Joint_Cartesian_Space

<mark>The final project has been done for fulfilment of credits for the *Robotics II course* given by *Prof. Alessandro De Luca*, A.Y. 2019-2020. </mark>

This project deals with the implementation of the algorithm presented in the paper: [*Physical human-robot interaction under joint and Cartesian constraints*](https://ieeexplore.ieee.org/document/8981579). We report the limitations of the original algorithm. Not only this the algorithm in the paper has been modified and a new proposed algorithm has been presented which efficiently overcomes the drawbacks of the original algorithm. 
<br />
**Implementing a simple Cartesian tracking controller**
<br />
Result of cartesian tracking control loop. (Not tunned(Left); Tracking output with duration T = 10(Center); Tracking output with duration T = 20(Right))
![alt text](images/cartesian_tracking_circular.png)
<br />
**Vrep result** of the Original algorithm (joint limit constraints having higher priority than Cartesian limit constraints)vs Proposed algorithm (no priority ordering of joint limit constraints and Cartesian limit constraints; all constraints assigned same priority) 
<br />
* Original Algorithm
![alt text](images/Orig.png)

* Proposed Algorithm
![alt text](images/Prop.png)

**Graphical result** of the Original algorithm vs Proposed algorithm

* Original Algorithm
![alt text](images/Orig_graph.png)

* Proposed Algorithm
![alt text](images/Prop_graph.png)

**Contributors**

The project has been completely developed by *Akash Garg* and *Amirhossein Kazemipour*