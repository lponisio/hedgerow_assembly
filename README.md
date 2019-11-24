analysis code for "Opportunistic attachment assembles plant-pollinator networks" Ponisio, Gaiarsa, and Kremen. See write up within analysis folder for details on repository organization and purpose of each script. The data accompanies the code for reproducibility, but we ask to be contacted (lponisio@gmail.com and ckremen@zoology.ubc.ca) if authors wish to use the data for a publication.

In our study we examin the temporal dynamics of plant-pollinator network assembly using variety of different methods including 1) network change point detection 2) node/species-level position variation, 3) species and interaction turnover 4) network-level metrics and 5) extinction simulations. We are committed to reproducible science and all analytical code will be maintained on github, along with the write up.

The entire analysis is executable from the main.sh file. All of the packages needed to run the analyses are listed in the packages.sh file, including the specific versions of the python dendropy library, which is very finicky. We have had no problems with backward compatibility for the R scripts (in terms of R versions from v.~2-3.4 or packages) but the change point analysis must be run in python2 (not easily convertible to python3, though this is potentially in progress).

Navigate to the analysis folder within the github repo (hedgerow_assembly) then the main.sh file can be selected and run. I would not recommend doing this unless you are prepared to have two cores of your machine running for days to weeks (depending on how many of iterations of the change point analysis specific in mainChangePoint.sh loop), but you could run of the analyses in the study by running this line in BASH.

This will somewhat helpfully print the results of each analysis and re-create any accompanying figures.

I walk through each the main script for each analysis individually in the write-up.pdf
