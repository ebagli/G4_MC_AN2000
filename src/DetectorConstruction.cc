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
//
//

#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4CrystalExtension.hh"
#include "G4ExtendedMaterial.hh"
#include "G4LogicalCrystalVolume.hh"

#include "G4ChannelingMaterialData.hh"
#include "G4ChannelingOptrMultiParticleChangeCrossSection.hh"

#include "CrystalDetector.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():
fECfileName(""),
fECOfileName(""),
fECAlfileName(""),
fMaterialName(""),
fSizes(G4ThreeVector(0.,0.,0.)),
fAngles(G4ThreeVector(0.,0.,0.)){
    fMessenger = new DetectorConstructionMessenger(this);
;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::DefineMaterials(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){

    //** World **//
    G4Material* Galactic = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

    G4Element* elIn = G4NistManager::Instance()->FindOrBuildElement(49);
    G4Element* elP = G4NistManager::Instance()->FindOrBuildElement(15);
    G4Material* InP = new G4Material("GECO_InP",4.81 *CLHEP::g/CLHEP::cm3,2);
    InP->AddElement(elIn,1);
    InP->AddElement(elP,1);
    
    G4Element* elO = G4NistManager::Instance()->FindOrBuildElement(8);
    G4Element* elSi = G4NistManager::Instance()->FindOrBuildElement(14);
    G4Element* elLa = G4NistManager::Instance()->FindOrBuildElement(57);
    G4Element* elGa = G4NistManager::Instance()->FindOrBuildElement(31);
    G4Material* Langasite = new G4Material("GECO_Langasite",5.754 *CLHEP::g/CLHEP::cm3,4);
    Langasite->AddElement(elO,14);
    Langasite->AddElement(elSi,1);
    Langasite->AddElement(elLa,3);
    Langasite->AddElement(elGa,5);

    G4Element* elLi = G4NistManager::Instance()->FindOrBuildElement(3);
    G4Element* elNb = G4NistManager::Instance()->FindOrBuildElement(41);
    G4Material* LiNbO3 = new G4Material("GECO_LiNbO3",4.65 *CLHEP::g/CLHEP::cm3,3);
    LiNbO3->AddElement(elLi,1);
    LiNbO3->AddElement(elNb,1);
    LiNbO3->AddElement(elO,3);

    
    G4Isotope* O18 = new G4Isotope("O18",8,18);
    G4Element* elO18  = new G4Element("Oxygen 18","elO18",1);
    elO18->AddIsotope(O18,100.*CLHEP::perCent);
    
    G4Element* elAl = G4NistManager::Instance()->FindOrBuildElement(13);
    G4Material* Al2O3_O18 = new G4Material("GECO_Al2O3_O18",3.95 *CLHEP::g/CLHEP::cm3,2);
    Al2O3_O18->AddElement(elAl,2);
    Al2O3_O18->AddElement(elO,3);

    G4double worldSizeXY = 1. * CLHEP::meter;
    G4double worldSizeZ = 30. * CLHEP::meter;
    
    G4Box* worldSolid = new G4Box("world.solid",
                                  worldSizeXY/2.,
                                  worldSizeXY/2.,
                                  worldSizeZ/2.);
    
    G4LogicalVolume* worldLogic = new G4LogicalVolume(worldSolid,
                                                      Galactic,
                                                      "world.logic");
    
    G4PVPlacement* worldPhysical = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     worldLogic,
                                                     "world.physic",
                                                     0,
                                                     false,
                                                     0);
    
    
    //** Crystal solid parameters **//
    G4double crystalSizeXY = 50. * CLHEP::millimeter;
    G4double crystalSizeZ = fSizes.z();
    
    G4Box* crystalSolid = new G4Box("crystal.solid",
                                    crystalSizeXY/2.,
                                    crystalSizeXY/2.,
                                    crystalSizeZ/2.);

    //** Crystal Definition Start **//
    G4Material* mat;
    G4cout << "Material Name: " << fMaterialName << G4endl;
    if(fMaterialName == ""){
        mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");;
    }
    if(fMaterialName == "GECO_Langasite"){
        G4cout << "LANGASITE!!!!" << G4endl;
        mat = Langasite;
    }
    else if(fMaterialName == "GECO_LiNbO3"){
        G4cout << "LiNbO3!!!!" << G4endl;
        mat = LiNbO3;
    }
    else if(fMaterialName == "GECO_Al2O3_O18"){
        G4cout << "Al2O3_O18!!!!" << G4endl;
        mat = Al2O3_O18;
    }
    else{
        mat = G4NistManager::Instance()->FindOrBuildMaterial(fMaterialName);
    }
    
    
    G4ExtendedMaterial* Crystal = new G4ExtendedMaterial("crystal.material",mat);

    Crystal->RegisterExtension(std::unique_ptr<G4CrystalExtension>(new G4CrystalExtension(Crystal)));
    G4CrystalExtension* crystalExtension = (G4CrystalExtension*)Crystal->RetrieveExtension("crystal");
    crystalExtension->SetUnitCell(new G4CrystalUnitCell(5.43 * CLHEP::angstrom,
                                                        5.43 * CLHEP::angstrom,
                                                        5.43 * CLHEP::angstrom,
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi,
                                                        CLHEP::halfpi,
                                                        227));
  
    Crystal->RegisterExtension(std::unique_ptr<G4ChannelingMaterialData>(new G4ChannelingMaterialData("channeling")));
    G4ChannelingMaterialData* crystalChannelingData = (G4ChannelingMaterialData*)Crystal->RetrieveExtension("channeling");
    crystalChannelingData->SetFilename(fECfileName);

    G4LogicalCrystalVolume* crystalLogic = new G4LogicalCrystalVolume(crystalSolid,
                                                                      Crystal,
                                                                      "crystal.logic");
    crystalLogic->SetVerbose(1);
    //** Crystal Definition End **//

    
    G4RotationMatrix* rot = new G4RotationMatrix;
    if(fAngles.x()!=0.){
        rot->rotateX(fAngles.x());
    }
    if(fAngles.y()!=0.){
        rot->rotateY(fAngles.y());
    }
    
    new G4PVPlacement(rot,
                      G4ThreeVector(),
                      crystalLogic,
                      "crystal.physic",
                      worldLogic,
                      false,
                      0);
    
    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField(){
    G4LogicalVolume* crystalLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic");
    G4ChannelingOptrMultiParticleChangeCrossSection* testMany =
    new G4ChannelingOptrMultiParticleChangeCrossSection();
    testMany->AttachTo(crystalLogic);
    G4cout << " Attaching biasing operator " << testMany->GetName()
    << " to logical volume " << crystalLogic->GetName()
    << G4endl;
    
    G4VSensitiveDetector* crystaldetector = new CrystalDetector("/crystaldetector",fECOfileName,fECAlfileName);
    G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector);
    crystalLogic->SetSensitiveDetector(crystaldetector);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

