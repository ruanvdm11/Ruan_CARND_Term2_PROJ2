# CARND-Term-2-Project-2-Unscented-Kalman-Filter
[//]: # (Image References)
[image1]: https://raw.githubusercontent.com/ruanvdm11/Ruan_CARND_Term2_PROJ1/master/Reference_Images/Dataset_1_Lidar_and_Radar.JPG "Dataset1"
[image2]: https://raw.githubusercontent.com/ruanvdm11/Ruan_CARND_Term2_PROJ1/master/Reference_Images/Dataset_2_Lidar_and_Radar.JPG "Dataset2"
[image3]: https://raw.githubusercontent.com/ruanvdm11/Ruan_CARND_Term2_PROJ1/master/Reference_Images/Dataset_1_Laser.JPG "Dataset1 Laser"
[image4]: https://raw.githubusercontent.com/ruanvdm11/Ruan_CARND_Term2_PROJ1/master/Reference_Images/Dataset_1_Radar.JPG "Dataset1 Radar"

This project required the implementation of the Unscented Kalman Filter on a given set of data.
#### Build Platform: Ubuntu 16.04 LTS

The following images quickly illustrates that the submitted source code is within the accuracy bounds.

## Dataset 1: Both Laser and Radar being used

![alt text][image1]

## Dataset 2: Both Laser and Radar being used

![alt text][image2]

# Accuracy Difference Between Laser Only and Radar Only
Here is a quick view at the accuracy difference between Laser Only and Radar Only implementations on dataset 1.
## Dataset 1: Only Laser

![alt text][image3]

## Dataset 1: Only Radar

![alt text][image4]

The RMSE values of the Radar only run is much higher than that of the Laser Only run. This could be due to the fact that the Radar generated data in a polar coordinate system. This means that errors can be induced during data conversions.
