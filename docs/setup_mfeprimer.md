---
layout: default
nav_order: 1
parent: index
permalink: /setup-mfeprimer
---

# Setup MFEprimer
MFEprimer is used for dimer prediction in multiplex wormhole. This helper function automatically downloads, installs, and configures the MFEprimer binary when called, and subsequently finds and returns the binary location when the dependency is required. This *should* run behind-the-scenes, but help is provided below in case of failure.

## Usage
### Command-line syntax
```
mw-setup-mfeprimer
```
### Python syntax
```
import multiplex_wormhole as mw
mw.setup_mfeprimer()
```

## Manual Installation Instructions
1. Download the MFEprimer **v3.2.7** release that fits your operating system [here](https://github.com/quwubin/MFEprimer-3.0/releases).
2. Save the file to your multiplex_wormhole package directory (location can be found by running `pip show multiplex_wormhole`). If you cloned the repository directly from GitHub, save MFEprimer to `multiplex_wormhole/src/multiplex_wormhole`.
3. Unzip the download (if zipped).
4. From the command line, move into this location e.g. `cd multiplex_wormhole/` and change the permissions on the MFEprimer file to enable execution: `chmod +x mfeprimer*`.
5. Check that all of the above worked by running `mfeprimer dimer -h` from your command line interface. This should print out help information for MFEprimer.
