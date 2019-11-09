# attitude_control_system
Description:
Attitude Control and Reaction Control System

## **A. Installation:**
If you need to check path
run `get_path.py` \
If you want to install
run the following: `pip -m install -r requirements.txt`

## **B. Guideline:**
 **B0.  DO NOT:**
- Please NEVER use `sudo pip install` \
`pip install` gives you a permission error, it likely means you're trying to install/update packages in a system python, such as /usr/bin/python. Doing this can have bad consequences, as often the operating system itself depends on particular versions of packages within that Python installation. For day-to-day Python usage, you should isolate your packages from the system Python, using either virtual environments or Anaconda(conda) — I personally prefer conda for this. 
- If you are using python IDE or command lines \
use `python -m pip install <package>` rather than `$ pip install <package>`because the former is more explicit about where the package will be installed. 
- If you are using Juypter kernel \
use `import sys
!{sys.executable} -m pip install numpy` rather than `!pip install numpy` That bit of extra boiler-plate makes certain that you are running the pip version associated with the current Python kernel, so that the installed packages can be used in the current notebook. This is related to the fact that, even setting Jupyter notebooks aside, it's better to install packages using


 **B1. For code version control:**
-  Please use python 3's virtual environment in order to be consistent on packages version control

-  Please write down detailed description/messages when you commit; style: `"issue_number-description-purpose-file(s)_name_you_commited"`.  \
This makes it easier to trace code changes. You can find issue number here: https://github.com/kalufinnle/acs/issues

- Please remember to `git pull` every time before you start to edit any .py

**B2. For high level view of variable location and code structure:**
    `SUAVE 2.0` `UML`

**B3. For code management:**
    `Gitkraken` `Glo Board` `Atom`

**B4. For literature management and collaboration:**
     `Mendeley`
     
**B5. In-house function dependency:**
    `real_gas_model`

### **GitKraken signup below!**
https://www.gitkraken.com/invite/gN1gNQGw
    
## **C. Coding style:**
http://suave.stanford.edu/documentation/code_style.html
### **Naming Convention**
#### In terms of typography:

- any_variable_name - lower case with underscore. This includes working variables and instantiated objects.
- field_name - lower case with underscore Any field of an object should be lower case.
- function_name - lower case with underscore
- Class_Type - upper case with underscore. The underscores are chosen here to permit the inclusion of acronyms if needed and to maintain symmetry with field name styling.
- Package_Name - upper case with underscore. For example folders within the SUAVE package.

#### In terms of naming:
- Chunk similar field types under a containing field
- Bias names towards being specific
- Write out field name verbosely, but try to keep short


## **D. Anaconda:**
To list all discoverable environments, type: `conda info --envs`.

To activate enviroment, type: `conda activate <environmentName>`

To deactivate an environment, type: `conda deactivate`

To create an environment, type: `conda create --name <environmentName>`





## **E. Package and library dependencies before runing the scripts:**
- **Pint** \
-----Pint (General install)-----\
`pip install pint`\
-----Pint (Anaconda)-----\
`conda install -c conda-forge pint`
> You can check the installation with the following command:\
> pint.test()

- **CoolProp** \
-----Coolprop (General install)----- \
`pip install cmake`\
`pip install Cython`\
`pip install CoolProp`\
-----Coolprop (Anaconda install)-----\
`conda install pip six`\
`pip install cmake`\
`pip install -i https://pypi.anaconda.org/coolprop/simple coolprop`\
`pip install -i https://pypi.anaconda.org/coolprop/label/dev/simple coolprop`\
OR\
`conda install -c conda-forge coolprop` \
`pip install cmake`\
`conda install -c conda-forge/label/gcc7 coolprop` \
`conda install -c conda-forge/label/cf201901 coolprop`

- **SUAVE** \
`git clone https://github.com/suavecode/SUAVE.git` \
`cd SUAVE/trunk` \
`python setup.py install`
