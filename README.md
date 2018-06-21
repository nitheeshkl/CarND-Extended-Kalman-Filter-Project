# Udacity_CarND_ExtendedKalmanFilters
Term2-Project1: Extended Kalman Filters

This project impelments the Extended Kalman filters in C++ for RADAR and LIDAR
data. LIDAR data constitutes of x & y positions and RADAR data constitutes of
radius, angle and radial velocity. These values are combined using a linear
motion model in a 2D plane.

This project works with the Udacity's simulator, which reads input data from the
dataset file and the calculated RMSE is displayed on the simulator.

To run the project, first start the simulator and then run the generated
FusionEKF binary from the project.
