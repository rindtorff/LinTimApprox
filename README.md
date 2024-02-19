## Abstract

In the field of planning a public transportation network, the concept of integrating and simultaneously considering two optimization steps for public transportation is well documented in the literature. This approach ensures that the solution achieves the highest quality, but it also comes with the drawback of increasing the complexity of the model, making it challenging to solve for large-scale instances. To address this, metamodels are commonly employed, as they provide a means to abstract complex systems, and in addition, it becomes easier to analyze and simulate various scenarios without the need to delve into the intricate details of the underlying system. In this thesis, we will develop an approximation algorithm for the integrated line planning and timetabling problem, with a metamodel playing a crucial role. We will then analyze and apply this algorithm to evaluate its effectiveness.

## Datasets and LinTim Project Overview

All data sets have been obtained by the LinTim-Project
LinTim is an open-source software toolbox that has the ability to solve various planning tasks in public transportation. The project was initiated by Anita Schöbel in 2007 at the University of Göttingen, then transferred to TU Kaiserslautern in 2019 and subsequently transferred to Philine Schiewe in 2022. It also offers already implemented algorithms for line planning, timetabling and more problems, but these have only been used to generate the data.

### References

For more information on LinTim and its contributions, the following references provide detailed documentation and insights:

- Schiewe, P., Schöbel, A., Jäger, S., Albert, S., Baumgart, U., Biedinger, C., Grafe, V., Roth, S., Schiewe, A., Spühler, F., et al. (2023). Documentation for LinTim 2023.12.

- Schiewe, P., Schöbel, A., Jäger, S., Albert, S., Baumgart, U., Biedinger, C., Grafe, V., Roth, S., Schiewe, A., Spühler, F., Stinzendörfer, M., & Urban, R. (2024). LinTim - Integrated Optimization in Public Transportation. Available online: [LinTim Website](https://www.lintim.net), Accessed in 2024.

- Goerigk, M., Schachtebeck, M., & Schöbel, A. (2013). Evaluating line concepts using travel times and robustness: Simulations with the LinTim toolbox. Public Transport, 5, 267-284.

## Getting Started

### Prerequisites

Before you begin, ensure you have the following installed:
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) (Anaconda or Miniconda)

### Installation

1. Clone the repository to your local machine:
   ```
   git clone https://github.com/rindtorff/LinTimApprox.git
   ```

2. Navigate to the repository directory:
   ```
   cd path/to/your-repository
   ```

3. Install and activate the provided Conda environment:
   ```
   conda env create -f environment.yml
   conda activate meta
   ```

### Configuration


Open `parameters.py` in your preferred text editor and modify the parameters as needed.

### Usage

To run the algorithm, execute the following command in your terminal:

```
python run.py
```

## Results Overview

Upon successful execution of the algorithm, the output will be organized into two main directories within the project folder:

### `results/`
This directory houses all the `.csv` files that contains the calculated results from the algorithm.

### `plots/`
When the `PLOT_RESULTS` parameter is set to `True` in the configuration file, this directory will be populated with visual representations of the results.