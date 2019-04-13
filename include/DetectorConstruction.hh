//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
// --------------------------------------------------------------
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#endif

#include "G4VUserDetectorConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "DetectorConstructionMessenger.hh"

#include "globals.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    
    DetectorConstruction();
    ~DetectorConstruction();
    
    void DefineMaterials();
    G4VPhysicalVolume* Construct();

private:
    void ConstructSDandField();

private:
    DetectorConstructionMessenger* fMessenger;

public:
    void SetEC(G4String aString) {fECfileName = aString;}
    void SetECO(G4String aString) {fECOfileName = aString;}
    void SetECAl(G4String aString) {fECAlfileName = aString;}
    void SetMaterial(G4String aString) {fMaterialName = aString;}
    void SetSizes(G4ThreeVector a3vec) {fSizes = a3vec;}
    void SetAngles(G4ThreeVector a3vec) {fAngles = a3vec;}

    G4String GetEC() {return fECfileName;}
    G4String GetECO() {return fECOfileName;}
    G4String GetECAl() {return fECAlfileName;}
    G4String GetMaterial() {return fMaterialName;}
    G4ThreeVector GetSizes() {return fSizes;}
    G4ThreeVector GetAngles() {return fAngles;}

private:
    G4String fECfileName;
    G4String fECOfileName;
    G4String fECAlfileName;
    G4String fMaterialName;
    G4ThreeVector fSizes;
    G4ThreeVector fAngles;
};


