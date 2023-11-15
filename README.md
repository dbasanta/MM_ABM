# MM_ABM
Multiple Myeloma Mathematical Model
This code is written in Java so requires a recent version of the JDK which can be freely downloaded from here: https://www.oracle.com/java/technologies/downloads/

We built our Agent-Based Model using the Hybrid Automata Layer, a library that facilitates the development of computationally efficient models. The library is developed by our colleagues at the Moffitt Cancer Center and can be downloaded from here: https://halloworld.org/start_here.html

HAL has extensive documentation and the approach was previously published here:
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007635

For convenience we recommend to use the IntelliJ development environment and to follow the instructions in HAL to integrate the library with it.

The code found here should be copied into the HAL-Master folder of the HAL library in the folder named ‘samples’. From there open IntelliJ, chose the BoneGrid_2022May17.java class and execute. The default parameters for the simulation will be read from params.csv and can be easily modified from there.
