# Installation

`PySMACKS` is provided as a Python package and is freely available from [https://pypi.org/project/PySMACKS/](https://pypi.org/project/PySMACKS/).

## Using pip

The suggested installation method is to use pip within a virtual environment.

```
pip install PySMACKS
```

You can then run the PySMACKS GUI by typing `PySMACKS_GUI` directly inside the console:

```py
PySMACKS_GUI
```

### Setting a Python virtual environment

To create a virtual environment, you can type the following in your Python console:

1. Select a folder:
  
    ```
    cd PATH_TO_THE_FOLDER_WHERE_TO_INSTALL_VENV
    ```
  
    with `PATH_TO_THE_FOLDER_WHERE_TO_INSTALL_VENV` replaced by the path to the folder where you want to install your virtual environment.

2. Create your virtual environment:

    ```
    python -m venv NAME_OF_YOUR_VENV
    ```

    with `NAME_OF_YOUR_VENV` replaced by the name you want to give to your virtual environment (for instance `PySMACKS_venv`)

3. Activate your virtual environment:

    ```
    NAME_OF_YOUR_VENV\Scripts\activate
    ```