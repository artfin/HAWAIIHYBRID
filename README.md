# Hawaii program for running trajectory-based calculations

# Quick Start 

### 1. Configure your build

Customize your personal build configuration
```
cp Makefile.config.example Makefile.config
```

Edit include and link paths in ```Makefile.config```.

### 2. Compile components 
Compile the driver and required modules with potential energy and dipole functions for your use case.
Let's taken CO2-Ar as an example:
```
make driver.exe
make build/ai_ids_co2ar.so # dipole moment dynamic lib
make build/ai_pes_co2ar.so # potential energy dynamic lib
```

### 3. Run simulation

Execute driver with the configuration file specifying the parameters of the the simulation:
```
./driver.exe ./examples/correlation_co2_ar.conf
```
