# ClusterPatientsDrugTrajectories

For English, see below.

## Tutvustus
Repositoorium on loodud bakalaureusetöö “Patsientide ravimikasutuse klasterdamine ATC koodide alusel” raames.

Projekti eesmärk on leida patsientide ravimitrajektoorid ja klasterdada need sarnasuse alusel. Trajektoorid koostatakse ravimite ATC koodidest. Trajektooride sarnasus arvutatakse eukleidilise kauguse kui ka dünaamilise ajadeformatsiooni algoritmiga.

## Failid
Repositooriumis on kuus faili:
* R/functions.R, milles on kõik töövoo R-i funktsioonid;
* src/functions.cpp, milles on kõik töövoo C++-i funktsioonid;
* shiny/ui.R, milles on defineeritud rakenduse välimus;
* shiny/server.R, milles on defineeritud rakenduse käitumine;
* shiny/app.R, milles on defineeritud shiny äpp;
* tmp/datasets/example_data.csv, milles on näidisandmed, millega on võimalik rakendust kasutada.

## Töökeskkond
Töövoo käivitamiseks on vajalik R-i arenduskeskkond. Lisaks peavad olema installeeritud järgnevad R-i moodulid:
* Rcpp - 1.0.6
* dplyr - 1.1.4
* pheatmap - 1.0.12
* gridExtra - 2.3
* flashClust - 1.01-2
* shinyWidgets - 0.8.0
* shinycssloaders - 1.0.0
* plyr - 1.8.6
* grDevices - 4.0.3
* grid - 4.0.3
* shiny - 1.8.0

C++ päisefailidest kasutatakse järgnevaid:
* Rcpp.h
* map
* string
* limits

## Kasutamine
Töövoogu on võimalik käivitada käsurealt ja rakenduses.

Käsurealt käivitamiseks tuleb läbida järgmised sammud:
1. Lae alla repositoorium ja liigu kausta.
2. Defineeri failis R/functions.R endale sobivad muutujad kas funktsioonis *run_workflow* VÕI *only_cluster_patients*. Esimene funktsioon käivitab kogu töövoo aga teine teeb läbi ainult patsientide klasterdamise.
3. ```source(“./R/functions.R”)```
4. Terve töövoo käivitamiseks ```run_workflow()```
   VÕI
   ainult patsientide klastedamiseks ```only_cluster_patients()```

Rakenduses töövoo käivitamiseks tuleb läbida järgmised sammud:
1. Lae alla repositoorium ja liigu kausta.
2. ```source(“./R/functions.R”)```
3. ```runApp(“./shiny/app.R”)```
4. Avanenud rakenduses defineeri endale sobivad parameetrid ja vajuta nuppu “Salvesta”. See käivitab töövoo uuendatud parameetritega.

## Description
The repository was created for the bachelor's thesis "Clustering patients' drug usage based on ATC codes".

The project aims to find drug trajectories of patients and cluster them based on similarity. The trajectories are made from the ATC codes of the drugs. Trajectory similarity is calculated using Euclidean distance as well as dynamic time warping.

## Files
In the repository, there are six files:
* R/functions.R, which has all the R functions of the workflow;
* src/functions.cpp, which has all the C++ functions of the workflow;
* shiny/ui.R, which defines the client side of the application;
* shiny/server.R, which defines the server side of the application;
* shiny/app.R, which defines the shiny app;
* tmp/datasets/example_data.csv, which has example data that is used to run the app.

## Environment
An R development environment is required to run the workflow. In addition, the following R modules must be installed:
* Rcpp - 1.0.6
* dplyr - 1.1.4
* pheatmap - 1.0.12
* gridExtra - 2.3
* flashClust - 1.01-2
* shinyWidgets - 0.8.0
* shinycssloaders - 1.0.0
* plyr - 1.8.6
* grDevices - 4.0.3
* grid - 4.0.3
* shiny - 1.8.0

The following C++ header files are used:
* Rcpp.h
* map
* string
* limits

## Usage
Workflow can be run through the command line and application.

To run the workflow through the command line, the following steps must be followed:
1. Download the repository and navigate into the folder.
2. In the R/functions.R file, define the appropriate variables in the function *run_workflow* OR *only_cluster_patients*. The first function runs the entire workflow, while the second only performs patient clustering.
3. ```source(“./R/functions.R”)```
4. To run the entire workflow ```run_workflow()```
   OR
   to only cluster patients ```only_cluster_patients()```

To run the workflow through the applications, the following steps must be followed:
1. Download the repository and navigate into the folder.
2. ```source(“./R/functions.R”)```
3. ```runApp(“./shiny/app.R”)```
4. In the opened application, define the appropriate parameters and press the "Salvesta" button. This will start the workflow with the updated parameters.
