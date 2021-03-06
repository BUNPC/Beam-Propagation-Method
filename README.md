# Beam-Propagation-Method
Our model based on the beam propagation method calculates the wavefront propagation in the scattering medium with the scattering mean free path and anisotropy factor characterized.

To run the code you need to:
1. Download all the Matlab files into one folder
2. Run demo.m
3. To modify the optical properties of the biological medium and the incident beam geometry, change the corresponding parameters in ModelSetting.m

After you run the demo.m provided, the input and output wavefront and xz cross section profiles should look like this:

![Figure](InputWavefront.jpg)
![Figure](OutputWavefront.jpg)
![Figure](crossSection.jpg)


For any issue reporting or suggestions, please contact Dr. Xiaojun Cheng, xcheng17@bu.edu

Notes:
1. The layer distance d should be smaller than or equal to the scattering mean free path ls to obtain a feasible beam profile in z.
2. The anisotropy factor g is not an input parameter but is related to sigma_x.
3. The model is able to calculate wave propagation for any incident wavefront. If you are interested in a special input wavefront, please contant us.

Updated 09/09/2019



