spectralColors
==============

Python utility to calculate the color of physical samples with spectra obtained from *first principles* calculations
such has produced by e.g. [Octopus](http://www.tddft.org/programs/octopus) TDDFT code.

The suite at the moment is composed of a single tool **color.py**. Another tool that calculates the spectral density 
from the molecular polarizability is soon to come.



## color.py

Calculate the color from a given spectral density.
The functionality of color.py are partially overlapping with colormath (http://code.google.com/p/python-colormath/).

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
