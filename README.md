This project explores Out-of-Time-Ordered Correlators (OTOCs) in the Bose-Hubbard model on a hexagonal lattice, which can help characterize quantum information scrambling in many-body quantum systems.

Table of Contents

About
Installation
Usage
Project Structure
Contributing
License
About

This repository contains code for computing OTOCs using Julia and Python. The hexagonal lattice setup highlights the scrambling of quantum information in the Bose-Hubbard model.

Installation

To set up the Julia environment:

Clone the repository:
bash
Copy code
git clone https://github.com/tripathydevjyoti/OTOC_Bose_Hubbard.git
Install dependencies specified in Project.toml:
julia
Copy code
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Usage

The main computations are in the src folder, and the Python code used for benchmarks is in the python files folder.

Project Structure

src: Julia scripts for OTOC computations.
python files: Legacy Python scripts for benchmarking.
benchmark: Code for plotting otocs.
Contributing

Contributions are welcome! Please open an issue or a pull request for any improvements or bug fixes.

License

This project is open source. See LICENSE for details.
