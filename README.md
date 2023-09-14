# Automated NDSRIs CPCA Digitalized Tool

We have developed a web-based application that can autonomously analyze the N-nitrosamine risk category of compounds from their SMILES notation, providing instant screening to identify high-risk formations of nitrosoamines, and the accuracy of the tool was validated using an FDA dataset of compounds with known nitrosoamine risks.

# Example Application

https://ndsris-cpca-in-silico-tool.streamlit.app/

## Installation

(0) Create anaconda environment

```
conda create --name nitro_env python=3.7.5
```

(1) Install required packages

```
conda activate nitro_env
pip install -r requirements.txt
```

(2) Once the pacakges are ready, you can set-up the streamlit web app

```
streamlit run .\test_app.py
```

### Batch Calculation

```
Directly use the combine_all_rules function in CPCA_rules.py
