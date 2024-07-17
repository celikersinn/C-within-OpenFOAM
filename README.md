# Cpp-within-OpenFOAM

# OpenFOAM C++ Studies

This repository contains my studies and implementations within the OpenFOAM environment. These studies are inspired by Håkan Nilsson's work at Chalmers University and guided by various YouTube tutorials.


## Introduction

This repository contains my OpenFOAM C++ studies, focusing on advanced implementations of boundary conditions and modified solvers. Future updates will include new turbulence modelling, function objects, and more.

## Getting Started

These instructions will help you set up your local copy of the project for development and testing purposes.

### Prerequisites

- OpenFOAM
- A C++ compiler (e.g., GCC)
- Git

### Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/your-username/Cpp-within-OpenFOAM.git
    ```

2. Navigate to the project directory:

    ```bash
    cd Cpp-within-OpenFOAM
    ```

3. Follow the specific setup instructions in each subdirectory's README.

## Files and Directories

- `boundaryConditions/` - Includes code and run files for boundary conditions.
  - `code/` - Contains several implementations of boundary conditions (e.g., `myvelocity`) and a make file.
  - `run/` - Contains test cases for the boundary conditions
- `Solvers/` - Modified solvers.
  - `code/` - Contains several modified solvers
  - `run/` - Contains test cases for new solvers
- More directories to be added soon...

## Usage

Detailed instructions on how to use the code can be found in each subdirectory's README. Here's a basic example to get you started:

1. Compile the boundary conditions library:

    ```bash
    cd boundaryConditions/code/myvelocity
    wmake libso
    ```
2. To use the new boundary conditions, go to one of the test cases located in `boundaryConditions/run`. Open one of the cases and add the following line to the `controlDict` file:

    ```plaintext
    libs (
        "libmyvelocity.so"
    );

## Future Work

- Implementation of new turbulence modelling
- Development of function objects
- Additional OpenFOAM utilities and enhancements

## References

- Håkan Nilsson's [website](https://sites.google.com/view/prof-hakan-nilsson/home)
- YouTube tutorials:
  - [Advance CFD with Pratyush](https://www.youtube.com/@advancecfdwithpratyush7162)
  - [CFD with OpenFOAM](https://www.youtube.com/@CFD-with-OpenFOAM)

## Contributing

Contributions are welcome! Please fork the repository and create a pull request with your changes.

## Acknowledgements

- Håkan Nilsson, Chalmers University
- Pratyush, Advance CFD with Pratyush
- CFD with OpenFOAM YouTube channel
