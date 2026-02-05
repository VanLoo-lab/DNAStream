# Getting started

## Dependencies
`DNAStream` has the following dependencies. These will be automatically installed during installation. *Pinned versions to be determined later.*  
- `h5py`     
- `numpy`      
- `pandas`   


## Installation

Create a `conda/mamba` environment (recommended) and install the package from the Github tagged release. 

```bash
conda create -n dnastream python=3.11
conda activate dnastream 

# With optional docs dependencies (reccommended)
pip install "dnastream[docs] @ git+https://github.com/VanLoo-lab/DNAStream.git@v0.1.0-beta"

#Just the DNAStream package
pip install "dnastream @ git+https://github.com/VanLoo-lab/DNAStream.git@v0.1.0-beta"

```


Verify the installation. 

```
python -c "import dnastream; from dnastream import DNAStream; print('dnastream', dnastream.__version__)"
```

 Package is ready to use if no errors occurred!

### Optional dependencies

To view the documentation locally:  
-  `mkdocs>=1.5`  
- `mkdocs-material>=9`    
- `mkdocstrings[python]>=0.25`    

```bash
pip install ".[docs]"
```

To run the test suite:
  - `pytest>=7`

```bash
pip install ".[test]"
```

For developers, all of the above dependencies plus:
  - `black>=24`

```bash
pip install -e ".[dev]"
```




# Unit tests
If the optional dependencies are installed for the test suite, then the package can be tested with:

```bash
pyest
```