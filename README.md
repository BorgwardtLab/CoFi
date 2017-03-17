# CoFi
Visualization and analysis of time series with history information

**Coordination Finder (CoFi) is a program motivated by the recent work of the Schroeder group of ETHZ. CoFi takes protein expression levels of single-cell imaging as input and provides methods to extract information about clusteirng and coordination.**

Coordination Finder does the following:
* Cluster time series data. It may consist of multiple, independent dimensions.
* Visualize history information from mouse embryonic stem cell lineages.
* Computation of changepoints in the time series.
* Provide expert feedback to measure the quality of the algorithm.

Coordination Finder runs on both Windows and Mac OSX. It can be run from the console. The requirements are Python and the PyQT4 package. By default, there is one test data set ready (*Simlation.csv*).

## Starting CoFi

Download and unzip the folder `Coordination Finder` and navigate into it. Then run the following commands from the console:

For **Mac**:
```
python -m pip install PyQt4
python CoordinationFinder.py
```

For Python 3, install MacPorts by downloading the correct package from https://guide.macports.org/chunked/installing.macports.html.

Then:
```
sudo port install py34-pyqt4
echo 'export PATH=/usr/local/bin:$PATH' >>~/.bash_profile
python CoordinationFinder.py
```

This will not work for MaxOSXSierra since PyQt4 is not supported anymore. Try PyQtX.  

For **Windows**:
```
pip install PyQt4
start CoordinationFinder.py
```

For Windows 7, an additional SIP installation is required.  


## Loading *Simulation.csv*

To load *Simulation.csv* it takes four steps:

1. Start CoFi. The main window will get displayed full screen.
2. Go to `File->Open` to show the Loading screen. A recommended sample selection of what to compute will appear. For later use, bear in mind that this configuration can be modified. For now, press `OK` to continue.
3. In the Opening-dialog, there is the option to select which columns are to be evaluated. The correct columns for the test data set have been preselected. To use a different data set, the selection has to be adjusted. Please press `OK` to start the computation.
4. *Simulation.csv* will be evaluated. Please be patient, this may take up to 5 min.

The test data set *Simulation.csv* consists of two cell lineages ('trees') with the simulated expression profiles of two Proteins: ProteinA and ProteinB. The first tree, *1nocp*, consists of time series with no changepoint. In *1cp*, every curve has a changepoint of fixed strength at a random position.
