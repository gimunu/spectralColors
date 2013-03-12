Spectral Colors
==============

A python utility suite to calculate the color of physical samples with spectra obtained from *first principles* calculations
such has produced by e.g. the [Octopus](http://www.tddft.org/programs/octopus) TDDFT code.

The suite at the moment is composed of a single tool **color.py** only. Another tool that calculates the spectral density 
from microscopic molecular polarizability is soon to come.



## color.py

Calculates the color of a light source (also reflection and transmission) from a given spectral density.
The functionality of color.py are in partial overlap with [colormath](http://code.google.com/p/python-colormath/).
For this reason is not excluded the possibility to include colormath as required dependency in the future.

External dependencies:
* nupy
* scipy
* matplotlib (optional)
* colormath (otptional)


####Usage

Calculate the color of a spectral distribution contained in ```file```:  
```color.py -f file```  
or from standard input  
```color.py < file```  

For more info and complete list of features:  
```color.py --help``` 
