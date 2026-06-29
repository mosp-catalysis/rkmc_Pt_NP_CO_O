
### Kinetic Monte Carlo (KMC) Simulations for Pt Nanoparticle (NP) Systems under CO Oxidation Conditions  

This package provides the Fortran source code and input files required to perform KMC simulations of Pt nanoparticle systems under CO oxidation conditions. Two Pt NP models are included: Pt(110)-(1×1) and Pt(110)-(1×2).  

The `input` file defines the parameters required for the simulation, including:  
* Number of KMC steps  
* Lattice space size and spacing  
* System temperature  
* Partial pressures of CO and O<sub>2</sub> gas  

The structural files for the two NP models are provided in the following folders:
- `110-1x1-NP-model-files`: Pt(110)-(1×1) NP model
- `110-1x2-NP-model-files`: Pt(110)-(1×2) NP model

Each folder contains the following `.xyz` files:  

- `ini.xyz`: Initial NP structure file, which can be customized by the user. In this package, it corresponds to either the Pt(110)-(1×1) or Pt(110)-(1×2) NP structure.  
- `100.xyz`, `110.xyz`, `111.xyz`, and `edge.xyz`: Files defining the positions and numbers of the (100), (110), (111), and edge sites on the corresponding NP.  

To compile the KMC code, run:

```bash
gfortran rkmc_Pt_CO_O.f95 -o rkmc.exe
```

To run a KMC simulation, place the compiled executable `rkmc.exe`, the `input` file, and all `.xyz` files from one of the model folders into the same working directory. Then execute:

```bash
./rkmc.exe
```

---

### Output Files

After the simulation, the following files will be generated:  

- **`last_one.xyz`** This file provides the final structure. The last column records the state of the site: **0** represents an empty surface site, **1** represents a CO-occupied surface site, **2** represents a O-occupied surface site, and **5** represents a bulk site.  

- **`atom_str_00.xyz`** This file captures the structural evolution during the simulation with a recording interval of 10,000,000 steps. Due to the large storage requirements of trajectory files, only atoms with a coordination number (CN) < 12 are recorded here. 

- **`step_rec_00.dat`** This file provides time and event statistics during the simulation with a recording interval of 100,000 steps. The data columns represent:  
  1. Time (in s)  
  2. Steps  
  3. Placeholder (0)  
  4. Number of CO adsorption events  
  5. Number of CO desorption events  
  6. Number of O<sub>2</sub> adsorption events  
  7. Number of O<sub>2</sub> desorption events  
  8. Number of CO diffusion events  
  9. Number of O diffusion events    
  10. Number of CO oxidation reaction events  

- **`Ea_COO_00.dat`** This file provides detailed information regarding each CO oxidation reaction event. The data columns represent:  
  1. Steps   
  2. Time (in s)  
  3. Activation energy of the CO oxidation reaction event  
  4. Lattice ID of the involved CO  
  5. Lattice ID of the involved O  
  6. CN of the Pt atom where CO is adsorbed  
  7. CN of the Pt atom where O is adsorbed  
  8. Generalized coordination number (GCN) of the Pt atom where CO is adsorbed  
  9. GCN of the Pt atom where O is adsorbed  

- **`coverage_00.dat`, `coverage100_00.dat`, `coverage110_00.dat`, `coverage111_00.dat`, and `coveragedge_00.dat`** These files provide information on the CO and O coverage changes during the simulation with a recording interval of 100,000 steps. The different filenames correspond to the total coverage and specific coverage on the (100), (110), (111) facets, and edge sites. The data columns represent:  
  1. Time (in s)  
  2. Steps  
  3. Number of adsorbed CO  
  4. Number of adsorbed O  
  5. CO coverage  
  6. O coverage  

