# Astro-utilities
Plotting tools
The code inputs a file of the form 
```
 1            1183300.00   .000000  56485.0391   .691747    176.140839   35.274055   35.734264
 ```
Each line represents a D-dimensional vector.

Upon opening the file from the file dialog, a list of entry widgets is created. There are some predefined column names, but they can be changed. 
If the column name is deleted from this list, it won't be plotted. Adding a column name plots the column. The default column name output is LaTeX

The following two columns of entries represent lower and upper limits for each column, respectively. 
Thereby generated cuts in the data can be saved and loaded from the file menu.

An output folder can be chosen from the file menu as well. Not choosing the output folder selects the default folder /results.

Clicking run generates the histogram in the folder results along with 2D histograms of projections of the data on different axes. Contours show are the 1 sigma and 2 sigma contours.
