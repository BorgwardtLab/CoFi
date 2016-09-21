## CoFi
Visualization and analysis of data sets containing time series with history information

**Coordination Finder (CoFi) is a program motivated by the recent work of the Schroeder group \cite{Schroeder}. CoFi takes protein expression levels of single-cell imaging as input and provides methods to extract information from the data.**
Coordination Finder does the following:
* Cluster time series data. It may consist of multiple, independent dimensions.
* Visualize history information from mouse embryonic stem cell lineages.
* Computation of changepoints in the time series.
* Provide expert feedback to measure the quality of the algorithm.
Coordination Finder runs on both Windows and Mac OSX. It can be run from the console. The requirements are having the newest version of Python installed as well as the PyQT4 package. By default, there is one test data set ready (*Simlation.csv*)

#Starting CoFi
Download and unzip the folder 'Coordination Finder' and navigate into it. Then run the following commands from the console:
For Mac:
'''
python -m pip install PyQt4
python CoordinationFinder.py
'''
For Windows:
'''
pip install PyQt4
start CoordinationFinder.py
'''

