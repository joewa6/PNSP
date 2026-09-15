# Example PET model

This directory contains the canonical dry −50e PET nanoplastic model used by the example workflow.

- `NP_correct.lib`: Amber OFF residue library; the complete particle is stored as unit `NP`.
- `NP_charged.frcmod` and `NP_neut.frcmod`: PET force-field additions.
- `NP_dry.pdb`: dry coordinate reference.
- `pet.prmtop` and `pet.inpcrd`: standalone dry-particle topology and coordinates for inspection and validation. The system builder uses the OFF library, force-field additions, and PDB rather than importing this topology.

SHA-256 checksums:

```text
835f49920adc90ed6d3f68cfedb895a92d7852b585b8bc50bb267fc5d662999a  pet.inpcrd
ff0fa66076cc985376bda8fa31dde0d11bc5809dd2c69dc88f01c8a84f7dff31  NP_correct.lib
1710adc756374624fb1c6b1cf876765dd16ec09a6a37fbfa520ee78295883ede  NP_dry.pdb
5a72e9ffb67a1381eb344bf8b5be256c779b06dc0400d3ffb75edbad69f7eabd  pet.prmtop
e21b9f01fd9ea508be5bffeec44beddd183336d3ae5d543ab55d7ab8e57f3845  NP_charged.frcmod
40155faa8bba1b4f5fecd5a356dd14cb2aea5fb095d171915320714d3b94f297  NP_neut.frcmod
```
