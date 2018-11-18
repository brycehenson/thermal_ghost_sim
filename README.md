# Thermal_Ghost_sim
**Bryce M. Henson**   
Numerical proof of principle for thermal ghost imaging 

The goals of this project are
* Test if thermal ghost imaging is practicaly realizable in our system.
  * QE
  * FLux
* Come up with an algorithm for constructing fake data with arbitrary correlations
  * Thermal correlations
  * Arbitrary function
* Demonstrate the many ways that correlation functions may be calculated and compare their scaling
  * Radial point by point
  * x,y,z point by point
  * x,y,z bin and fft
  * see if the [fast binning algorithms](https://github.com/brycehenson/fast_search_based_histogram) i have previously developed can be helpfull 
  


Demonstrate that the spatio-temporal correlations between the ports of a beamsplitter with a thermal cloud as the input can be used to reconstruct the 


## Contributions  
This project would not have been possible without the many open source tools that it is based on. In no particular order: 
* ***James Conder*** [gaussfilt](https://au.mathworks.com/matlabcentral/fileexchange/43182-gaussfilt-t-z-sigma)
* ***Ander Biguri*** [Perceptually uniform colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
* ***Benjamin Kraus*** [nanconv](https://au.mathworks.com/matlabcentral/fileexchange/41961-nanconv)
* ***Daniel Eaton***    [sfigure](https://au.mathworks.com/matlabcentral/fileexchange/8919-smart-silent-figure)
* ***Denis Gilbert***    [M-file Header Template](https://au.mathworks.com/matlabcentral/fileexchange/4908-m-file-header-template)
