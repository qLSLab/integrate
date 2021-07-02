# INTEGRATE: model-based multi-omics data integration to characterize multi-level metabolic regulation

<img src="./images/workflow.pdf" width="100%"/>

## Overview
The study of metabolism and its regulation is finding increasing application in various fields, including health, wellness, and biotransformations. Complete characterization of regulatory mechanisms controlling metabolism requires knowledge of metabolic fluxes, whose direct determination lags behind other omic technologies, such as metabolomics and transcriptomics. In isolation, these methodologies do not allow accurate characterization of metabolic regulation. Hence, there is a need for integrated methodologies to disassemble the interdependence between different regulatory layers controlling metabolism.
To this aim, we propose a computational pipeline to characterize the landscape of metabolic regulation in different biological samples. The method integrates intracellular and extracellular metabolomics, and transcriptomics, using constraint-based stoichiometric metabolic models as a scaffold. We compute differential reaction expression from transcriptomic data and use constraint-based modeling to predict if the differential expression of metabolic enzymes directly originates differences in metabolic fluxes. In parallel, using metabolomic data, we predict how differences in substrate availability translate into differences in metabolic fluxes. By intersecting these two output datasets, we discriminate fluxes regulated at the metabolic and/or transcriptional level. This information is valuable to better inform targeted action planning in different fields, including personalized prescriptions in multifactorial diseases, such as cancer, and metabolic engineering.

## Installation

* Install Python (from release version 3.0 on)

* Install external libraries through the following commands:
  * for conda users: conda create --name myenv --file integrate_requirements.txt
  * for all other users: pip install -r integrate_requirements.txt

* Using terminal, navigate to the favourite installation directory and run the following Git command:
'git clone https://github.com/qLSLab/integrate.git'


## Usage

**Step 1: getGPRsFromModel**  
* Aim: Takes input model and returns the GPR rules associated to each reaction
* Usage: `python pipeline/getGPRsFromModel.py modelId`
* Inputs:
  * modelId: the input model name
   Users may also decided to set the following input files differently from the default imposed ones. Otherwise, these inputs can be left associated to the default values:
    * modelName: file name of the input generic model: Default value: 'ENGRO2_reversible_20210305'
* Output:
  * File containing for each reaction (column *id*) its GPR rule (column *rule*).

**Step 2: getRASscore**
* Aim: generate RAS starting from GPR rules and transcriptomics data
* Usage: `python pipeline/getRASscore.py rnaSeqFileName gprRule model regexOrgSpecific outputFileName`
* Inputs:
  * rnaSeqFileName: transcriptomics dataset
  * gprRule: output file of Step 1
  * model: the input generic model
  * regexOrgSpecific: define regex to extract genes from GPR rules. Dot is needed to also extract already computer scores, which are float
  * outputFileName: output file name
* Output:
  * File containing a column including the reactions ID (column *Rxn*), and a column for each cell line in the input transcriptomics dataset corresponding to the computer RAS score.

**Step 3: rasDfPreparation**
* Aim: compute the normalized RAS scores
* Usage: `python pipeline/rasDfPreparation.py inputFileName outputFileName`
* Inputs:
  * inputFileName: output file of Step 2
  * outputFileName: output file name
* Output:
  * File containing for each reaction (column *Rxn*) the mean (*mean_XXX* column) and normalized (*norm_XXX* column) RAS for each cell line *XXX* 

**Step 4: rasIntegration**
* Aim: integrate RAS scores within the input generic models to generate cell relative models
* Usage: `python pipeline/rasIntegration.py imposeYSI imposeMedium imposeRasConstraints rasNormFileName`
* Inputs:
  * imposeYSI: 'Y' (yes) or 'N' (no) according to whether extracellular flux ratio constraints have to be integrated
  * imposeMedium: 'Y' (yes) or 'N' (no) according to whether nutrients availability constraints have to be integrated
  * imposeRasConstraints: 'Y' (yes) or 'N' (no) according to whether transcriptomics derived constraints have to be integrated
  * rasNormFileName: output file name of Step 3
  Users may also decided to set the following input files differently from the default imposed ones. Otherwise, these files can be left associated to the default values:
    * ysiFileName: the input extracellular flux ratio dataset. This file consists of the following columns: "Ratio" listing the extracellular flux ratios to impose within the model; n columns for each cell line listing for each extracellular flux ratio the corresponding value in each replica.
    * mediumFileName: the input nutrients availability dataset. This file consists of the following columns: "Rxn" listing nutrients uptaken within the model; two columns for each cell line listing for each nutrient the corresponding lower bound (cell line + _LB) and upper bound (cell line + _UB)
    * modelName: the input generic metabolic model
    * prefixOutputModelName = the input model name
    * lcellLinesv list of input sample names
    * biomassRxn: biomass reaction name
    * lacRxn: the lactate demand reaction
    * glcRxn: the glucose uptake reaction
    * glnRxn: the glutamine uptake reaction
    * gluRxn: the glutamate demand reaction
    * lReplicas: list of the sample replicas
* Output:
  * A SBML model for each input cell line

**Step 5: Models splitting**
* Aim: Each model needs to be converted to a mat file in order to exploit the MATLAB function to convert model from the reversible into the irreversible format

**Step 6: randomSampling**
* Aim: sample the feasible flux region of each cell relative model
* Usage: `python pipeline/randomSampling.py nSamples`
* Inputs:
  * nSamples: number of solutions to sample
  Users may also decided to set the following input files differently from the default imposed ones. Otherwise, these files can be left associated to the default values:
    * biomassRxnName: biomass reaction name
    * lcellLines: list of input sample names
    * modelId: the input model name
    * timeStampInput: time stamp of the input models
    * ysiFileName: the input extracellular flux ratio dataset. This file consists of the following columns: *Ratio* listing the extracellular flux ratios to impose within the model; n columns for each cell line listing for each extracellular flux ratio the corresponding value in each replica.
    * lacRxn: the lactate demand reaction
    * glcRxn: the glucose uptake reaction
    * glnRxn: the glutamine uptake reaction
    * gluRxn: the glutamate demand reaction
    * lReplicas: list of the sample replicas
    * epsilon: lower bound imposed on the biomass synthesis reaction
* Output:
  * File containing the nSamples sampled solutions from each input sample model

**Step 7: mannWhitneyUTest**
* Aim: compute the Mann-Whitney U test
* Usage: `python pipeline/mannWhitneyUTest.py timeStampInput nSamples`
* Inputs:
  * timeStampInput: time stamp of the input random sampling datasets from Step 6
  * nSamples: number of sampled solutions in Step 6
  Users may also decided to set the following input files differently from the default imposed ones. Otherwise, these files can be left associated to the default values:
    * lcellLines: list of input sample names
    * modelId: the input model name
* Output:
  * For each pair of input cell line c_1 and c_2, a file returning the output of Mann-Whitney U test

## Cite Us
If you use INTEGRATE in your scientific work, please cite:
**AGGIUNGERE REF BIORXIVE**

## Getting Help
For support, please contact:
* marzia.difilippo@unimib.it
* chiara.damiani@unimib.it
* dario.pescini@unimib.it
* bruno.galuzzi@unimib.it
