# Developer guide for CalibratedCPP GitHub Repo

## Read first

- [LPhy developer guide](https://github.com/LinguaPhylo/linguaPhylo/blob/master/DEV_NOTE.md)
- [LPhyBEAST developer guide](https://github.com/LinguaPhylo/LPhyBeast/blob/master/DEV_NOTE.md)


## For CalibratedCoalescentPointProcessing developers

### Project structure
CalibratedCPP contains 3 subproject:

1. calibratedcpp-beast
2. calibratedcpp-lphybeast
3. calibratedcpp-lphy

### Maven build

With IntelliJ. open the sidebar with maven icon:

<a href="./figures/sidebar.png">
<img src="./figures/sidebar.png" width="200"/>
</a>

1. Click `Reload all Maven Projects` on the top left.
2. Open the first tab `Profiles` and tick skipLphyTests.
3. Open the second tab `calibratedcpp-root/lphybeast-root/linguaphylo/`, run `clean` and then `install` in `Lifecycle`.
4. Go to `calibratedcpp-root/lphybeast-root`, run `clean` and then `install` in `Lifecycle`.
5. Go to `calibratedcpp-root`, run `clean` and then `install` in `Lifecycle`.
Then you can find all `jar` files in the `target` folder for each subproject.

### Run with configuration

The following can run as long as CalibratedCPP Beast2 package is correctly installed.

### Lphy Studio
<a href="./figures/studio.png">
<img src="./figures/studio.png"/>
</a>

#### Lphy Simulation
<a href="./figures/Slphy.png">
<img src="./figures/Slphy.png"/>
</a>

#### LphyBeast Generation
<a href="./figures/lphybeast.png">
<img src="./figures/lphybeast.png"/>
</a>

#### BEAST2 Run
<a href="./figures/beastRun.png">
<img src="./figures/beastRun.png"/>
</a>
