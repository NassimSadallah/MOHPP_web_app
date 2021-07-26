This paper presents a multi-objective hybrid path planning method MOHPP for Unmanned Aerial Vehicles (UAVs) in urban dynamic environments.Two objectives are considered: the safety level and the travel time. First, we
construct two models of obstacles; static and dynamic. The static obstacles model
is based on Fast Marching Square (FM²) method to deal with the uncertainty of
the geography map, and the unexpected dynamic obstacles model is constructed
using the perception range and the safety distance of the UAV. Then, we developed
a jointly offline and online search mechanism to retrieve the optimal path. The
offline search is applied to find an optimal path vis-a-vis the static obstacles, while
the online search is applied to quickly avoid unexpected dynamic obstacles.

A script is added to control an UAV via dronekit api over the Mavlink communication.
