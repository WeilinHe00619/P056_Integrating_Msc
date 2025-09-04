# P056_Integrating_Msc

**Integrating Multi-Omics Datasets to Build Regulatory Networks and Identify Molecular Drivers of Micronutrient Deficiency in Pregnancy**

## Overview

This repository contains the complete codebase and analytical pipeline for an MSc research project focused on integrating multi-omics data with clinical features to investigate the molecular mechanisms underlying vitamin B12 deficiency during pregnancy. 

The project employs advanced statistical modeling, biological network reconstruction, and machine learning approaches to:

- Identify regulatory relationships across DNA methylation (CpG), mRNA, and miRNA molecular layers
- Develop interpretable multi-omics predictive models for micronutrient deficiency
- Discover consistent molecular drivers and biological pathways associated with B12 status

## Key Features

- 📊 **Multi-omics Integration**: Combines DNA methylation (RRBS), mRNA expression (RNA-seq), and miRNA expression (small RNA-seq) datasets
- 🧬 **Regulatory Network Modeling**: Infers CpG–mRNA and miRNA–mRNA regulatory interactions using biologically-informed statistical approaches  
- 🤖 **Machine Learning Pipeline**: Implements Logistic Regression, Random Forest, and Artificial Neural Networks with both data-driven and mechanism-driven feature selection strategies
- 🧠 **Model Interpretability**: Comprehensive feature importance analysis, cross-model consensus evaluation, and biological relevance assessment
- 📈 **Functional Analysis**: GO and KEGG pathway enrichment analysis of key genes identified across omics layers

## Repository Structure

```
P056_Integrating_Msc/
│
├── data/                    # Raw and processed multi-omics and clinical datasets
├── scripts/                 # Core analysis pipeline
│   ├── 1_preprocessing/     # Quality control, normalization, and ID mapping
│   ├── 2_differential/      # Differential analysis (DEGs, DMRs, DEmiRs)
│   ├── 3_network/           # Regulatory network construction and analysis
│   ├── 4_modeling/          # Machine learning model training and evaluation
│   └── 5_interpretation/    # Feature interpretation and pathway enrichment
├── results/                 # Analysis outputs: plots, tables, model results
├── figures/                 # Publication-ready figures and visualizations
├── requirements.txt         # Python package dependencies
└── README.md               # Project documentation
```

## Installation & Setup

### Prerequisites
- Python 3.8 or higher
- R 4.0 or higher (for certain statistical analyses)

### Installation Steps

1. **Clone the repository**
   ```bash
   git clone https://github.com/WeilinHe00619/P056_Integrating_Msc.git
   cd P056_Integrating_Msc
   ```

2. **Create virtual environment**
   ```bash
   # Using conda (recommended)
   conda create -n multiomics_env python=3.9
   conda activate multiomics_env
   
   # Or using venv
   python -m venv multiomics_env
   source multiomics_env/bin/activate  # On Windows: multiomics_env\Scripts\activate
   ```

3. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

4. **Verify installation**
   ```bash
   python scripts/4_modeling/run_model_rf.py --test
   ```

> **Note**: You may need to adjust data file paths in the configuration files based on your local directory structure.


### Individual Analysis Modules
Each analysis step can be run independently. See the `scripts/` directory for detailed usage instructions for each module.

## Main Findings

- **Predictive Performance**: Achieved robust classification of B12 deficiency status with interpretable multi-omics models
- **Key Molecular Features**: Identified consistent biomarkers across models including clinical features (BMI, Age) and molecular markers (CCNG1, IQCG)
- **Biological Insights**: Revealed enriched pathways related to metabolic regulation, immune response, and developmental processes
- **Methodological Innovation**: Demonstrated the effectiveness of both statistical and regulation-driven feature selection strategies


## Contributing

We welcome contributions to improve the analysis pipeline or extend the methodology. Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Weilin He**  
📧 Email: [te23071@bristol.ac.uk]  
🔗 GitHub: [@WeilinHe00619](https://github.com/WeilinHe00619)  
🏛️ Institution: [University of Bristol]

For questions about the methodology, data access, or collaboration opportunities, please don't hesitate to reach out.


*Last updated: September 2025*
