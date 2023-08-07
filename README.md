# bs-villungerlab
Repo for work completed as part of the VBC summer school in the villunger lab by Ben Silke.


# Folder Structure
```


```


## Python/
This directory contains all the python related code.

#### Metaprogramming
These scripts and classes contain functionality to write markdown files in bulk, providing the necessary arguments once.
The script will create a file which contains the code to run, which then can be run or knitted together.

## R/
This directory contains all the R related code.\
*utils.R files contain reusable functions.
Every folder contains specific treatment scripts as artefacts and a general documented script which can be used.



# Trouble Shooting


### My script is not running?
If it is python, you can either use `python3 __file__name__.py --argument_1 "1" --argument_2 "2"`
In R, you can use, `Rscript __script_name__.py --install.packages("argparser")`

otherwise, you should be able to run all scripts directly with `__script__name__.* --arguments`.
However, you will need to make these scripts executable. This can be done with this command.
`chmod +x __script__name__.py`