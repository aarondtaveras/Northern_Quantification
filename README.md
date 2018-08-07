# Northern_Quantification
A general purpose application that computes and quantifies useful information about Northern Blot Analyses.

# Build Instructions:
Run `g++ main.cpp` to produce your executable. 

Run the executable with the input data as follows: `./exe.out <data.csv> `. 

All given input data MUST be cleaned to achieve the proper results, and be in a `.csv` format. 

# Run-time Instructions: 

Use the following commands to compute statistics:
  `b/gapdh` - Computes band over gapdh per each trial.
  
  `tscore` - Computes t scores for band over gapdh per trial.
  
  `average < tscore || m/area || b/gapdh >` - Computes the average for the given argument. Arguments ARE case/space-sensitive.
  
  `extract <output.txt>` - Extracts current computed data to a text file named "results.txt". The file is overwritten by        default, and is written to "results.txt" by default.

# Cleaned Data:

Cleaned data refers to data that is in the proper row/column format. Individual cells are case insensitive, and ignore leading and trailing whitespace. The accepted format: 

  # Columns : Name (data type) : {set}
    Probe (string) : {"ITS1A","ITS2B","GAPDH"}
    Genotype (string) : {"KO","AB","WT"}
    Strip (string) : {"2a","1s","3a"}
    Band (string) : {"45","30","32","21","18SE","12s", empty}
    Lane (int) : {+I}
    Area (float) : {+R}
    Mean (float) : {+R}
    Min (int) : {+I}
    Max (int) : {+I} 
    Mean/Area (float) : {+R}
    
Feel free to submit bug reports or clarifications. I am by no means a biology expert - just a programmer. Happy computing! :)
