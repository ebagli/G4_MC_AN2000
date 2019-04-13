# GECO_INFN_AN2000
The application simulates the AN2000 line with the experimental setup used for the GECO-INFN-CSNV experiment.
The application is based on the Geant4 10.03 release. Please cite the Geant4 articles (http://www.geant4.org/geant4/), the Geant4-Channeling article (http://link.springer.com/article/10.1140%2Fepjc%2Fs10052-014-2996-y) if you use the application and the Geant4-Crystal article (submitted to NIMB).
- http://link.springer.com/article/10.1140%2Fepjc%2Fs10052-014-2996-y#page-1

The charged particles impinging on a crystal may undergo the channeling effect. Because of the channeled particle traverse the crystal with a different nuclear encounter probability with respect to an amorphous material, the inelastic interaction rate under channeling condition is modified.

## Experimental setup
In the application, the experimental line is composed only by a crystal. The crystal has to have an Oxygen component.
The crystal has a crystalline structure which permits to exploit the channeling effect.
The crystal object is created in the DetectorConstructor class and used by the G4Channeling process.
The electric field, electron and nuclei density under continuum potential approximation have to be loaded into Geant4 in order to enable the channeling effect.
The nuclei density of Oxygen components have to be loaded to compute the average density experienced by a particle interacting with the sub-lattice and the average cross-section for the O18(p,alpha)N15 reaction.

## Input file
The standard input file, which can be generated using the generator in the filenameGenerator folder, is usually composed by the following lines:

- `/xtal/setMaterial G4_ALUMINUM_OXIDE` - The crystalline material composition
- `/xtal/setSize 50. 50. 0.01 mm` - The crystal size
- `/xtal/setAngle 0. 0. 0. degree` - The orientation angle of the crystal
- `/xtal/setEC data/Al2O3100ax` - The data files for channeling
- `/xtal/setECO data/Al2O3_O100ax` - The data files for the detector
- `/run/initialize` #Initialize the detector and the physics
- `/gps/ene/mono 642.5 keV` - The beam energy
- `/gps/particle proton` - The beam particle
- `/gps/pos/type Point` - The beam position distribution
- `/gps/pos/centre 0 0 -633. cm` - The beam initial position (not to be changed)
- `/gps/direction 0 0 1` - The beam direction
- `/gps/ang/type focused` - The beam angular distribution
- `/run/beamOn 1000` - Run the simulation with 1000 particles

## Output file
The files are stored in CSV files. The first 8 lines are header lines.
The four column with doubles are for:
- Average nuclei density encoutered by the channeled particle. For an amorphous crystal this value is 1.
- Average electron density encoutered by the channeled particle. For an amorphous crystal this value is 1.
- Average Oxygen nuclei density encoutered by the channeled particle. For an amorphous crystal this value is 1.
- Average O18(p,alpha)N15 cross section in arbitrary unit for the channeled particle.
