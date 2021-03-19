## Introduction

QPAD requires an input file for initializing the simulation. When QPAD starts running, it will look for the input file named "qpinput.json" located in the current working directory and read parameters from this file. If the file does not exist in the same folder where the executable file is located, the program will stop with an error message.

## The Format of the Input File

The format of the input file is JSON. You can find more information about JSON at [here](http://www.json.org/). It should be noted that JSON does not allow trailing commas. In addition, JSON does not allow comment. However, in QPAD we use [JSON-FORTRAN](https://github.com/jacobwilliams/json-fortran) which can ignore all the comments before processing the file. In "qpinput.json", you can add any comment after an exclamation mark "!", and any text following this "!" will be ignored until the end of the line is reached.

## Sessions of Input File

In "qpinput.json", there are "simulation", "beam", "species", "neutrals" and "field" sessions. All parameters are using normalized units unless it is specified.