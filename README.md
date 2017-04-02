# README

The **Iterative Placement Algorithm** is an algorithm designed for placing the most regular polygon modules inside a given border.

## ITPLA

**ITPLA** firstly generates the maximum number of modules and places them randomly.
Then the algorithm starts evolving the placement.
After reaching some kinds of stable states,
a decision will be made whether continue evolving with removing some modules or stop and return the result.

## 1. Usage

The **Iterative Placement Algorithm** has two different realization versions:

1. **C++**

    This version is the alogrithm's prototype.
    It is based on the **C++** language.
    To realize the GUI, It introduces the **Qt** library.
    It also uses **CGAL** library to do 2d geometric caculation and **Box2D** library to simulate the movement.

    Because the major purpose is to validate the algorithm, this version has many disadvantages:

    - the runtime enviroment is not easy to setup
    - the type definitions and variable names are confusable
    - it's difficult to add new features for the poor modular design
    - the minimum distance calculation between two objects is complex and just an approximation

2. **Python**

    For the main task of the algorithm is 2d calculation, physical simulation and GUI,
    **Python** maybe a better choice than C++ for
    
    1. less caring about underlying data types
    2. more powerful and free features help to reduce codes
    3. easily used package manager **pip** and libraries **numpy**, **box2d** and **pygame**
    
    This version's goal is for simplifying and enhancing, and the work is still undergoing.
    The finished part has already shown some advantages:
    
    - now the enviroment setup is clear and easy
    - the variables are less and their names are meaningful
    - many simple functions are defined and used to simplify and clear the code
    - the general minimum distance calculation procedure is modular and accurate now
    - the lines of code significiantly reduce from 1000+ to 200~, the compact code helps debugging and reading
    
## 2. Documentation

## 3. Licensing

For further license information, please contact the authors.

## 4. Authors contacts

If you want to be informed about algorithm updates and new releases, obtain further information about the provided code, or contribute to its development please write to:

Xuenan Guo - xug015@ucsd.edu

Fulvio Mastrogiovanni - fulvio.mastrogiovanni@unige.it
