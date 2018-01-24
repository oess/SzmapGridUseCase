#!/usr/bin/env python
#############################################################################
# Copyright (C) 2017 OpenEye Scientific Software, Inc.
#############################################################################
# prep a mol for electrostatics: alts, hydrogens, charges, radii, split ligand
#############################################################################
import sys
from openeye.oechem    import *
from openeye.oequacpac import *

def WaterProcess(processName):
    if processName == "fullsearch":
        return OEPlaceHydrogensWaterProcessing_FullSearch
    elif processName == "focused":
        return OEPlaceHydrogensWaterProcessing_Focused
    return OEPlaceHydrogensWaterProcessing_Ignore

def main(argv=[__name__]):

    itf = OEInterface(InterfaceData)
    OEConfigureSplitMolComplexOptions(itf,
                                      OESplitMolComplexSetup_All              & \
                                  ~ ( OESplitMolComplexSetup_CovBondTreatment | \
                                      OESplitMolComplexSetup_CovCofactor ) )

    if not OEParseCommandLine(itf, argv):
        OEThrow.Fatal("Unable to interpret command line!")

    verbose = itf.GetBool("-verbose")
    if verbose:
        OEThrow.SetLevel(OEErrorLevel_Verbose)

    altProcess = itf.GetString("-alts")
    keepAlts     = not (altProcess == "a")
    highestOcc   = (altProcess == "occupancy")
    compareAlts  = (altProcess == "compare")

    siteNum    = itf.GetUnsignedInt("-bindingsitenum")
    allSites   = (siteNum == 0)
    otherModel = (itf.GetUnsignedInt("-modelnum") != 1)

    newHyd     = itf.GetBool("-newhyd")
    placeHyd   = itf.GetBool("-placehydrogens")
    addCharges = itf.GetBool("-charge")
    addRadii   = itf.GetBool("-radii")
    splitlig   = itf.HasString("-ligout")

    watProcessName = itf.GetString("-waterprocessing")
    waterProcess = WaterProcess(watProcessName)

    standardize = itf.GetBool("-standardizehyd")
    badclash  = itf.GetDouble("-clashcutoff")
    flipbias  = itf.GetDouble("-flipbias")
    maxStates = itf.GetDouble("-maxsubstates")

    flavor = OEIFlavor_PDB_Default | OEIFlavor_PDB_DATA
    if keepAlts:
        flavor = flavor | OEIFlavor_PDB_ALTLOC
    if otherModel:
        flavor = flavor & ~ OEIFlavor_PDB_ENDM

    ims = oemolistream()
    ims.SetFlavor(OEFormat_PDB, flavor)

    inputFile = itf.GetString("-in")
    if not ims.open(inputFile):
        OEThrow.Fatal("Unable to open %s for reading." % inputFile)

    if not OEIs3DFormat(ims.GetFormat()):
        OEThrow.Fatal("%s is not in a 3D format." % inputFile)

    inftype = OEGetFileType(OEGetFileExtension(inputFile))
    if (inftype == OEFormat_PDB) and not keepAlts:
        OEThrow.Verbose("Default processing of alt locations (keep just 'A' and ' ').")

    sopt = OESplitMolComplexOptions()
    OESetupSplitMolComplexOptions(sopt, itf)

    inmol = OEGraphMol()
    if not OEReadMolecule(ims, inmol):
        OEThrow.Fatal("Unable to read %s." % inputFile)

    ims.close()
    
    if newHyd:
      OEThrow.Verbose("Removing input hydrogens from %s." % inmol.GetTitle())
      for atom in inmol.GetAtoms(OEIsHydrogen()):
          inmol.DeleteAtom(atom)

    if inmol.NumAtoms() == 0:
        OEThrow.Fatal("Input molecule %s contains no atoms." % inputFile)

    if inmol.GetTitle() == "":
      inmol.SetTitle("input mol")

    OEThrow.Verbose("Processing %s." % inmol.GetTitle())

    if not OEHasResidues(inmol):
        OEPerceiveResidues(inmol, OEPreserveResInfo_All)

    if highestOcc:
        alf = OEAltLocationFactory(inmol)
        if not alf.GetGroupCount() == 0:
            OEThrow.Verbose("Dropping alternate locations from protein.")
            alf.MakePrimaryAltMol(inmol)
    elif compareAlts:
        inftype = OEGetFileType(OEGetFileExtension(inputFile))
        if (inftype == OEFormat_PDB):
            alf = OEAltLocationFactory(inmol)
            if not alf.GetGroupCount() == 0:
                OEThrow.Verbose("Repairing bonds to alts in %s." % inmol.GetTitle())
                inmol = alf.GetSourceMol()

    outmol = OEGraphMol()
    if allSites:
        outmol = inmol
    else:
        OEThrow.Verbose("Splitting out selected complex.")

        soptSiteSel = OESplitMolComplexOptions(sopt)
        soptSiteSel.SetSplitCovalent(False) # do any cov lig splitting later

        frags = OEAtomBondSetVector()
        if not OEGetMolComplexFragments(frags, inmol, soptSiteSel):
            OEThrow.Fatal("Unable to fragment %s." % inmol.GetTitle())

        howManySites = OECountMolComplexSites(frags)
        if howManySites < siteNum:
              OEThrow.Warning(("Binding site count (%d) " + \
                               "less than requested site (%d) in %s.") % \
                              (howManySites, siteNum, inmol.GetTitle()))
              exit(0)

        if not OECombineMolComplexFragments(outmol, frags, soptSiteSel):
            OEThrow.Fatal("Unable to collect fragments from %s." % inmol.GetTitle())

        if outmol.NumAtoms() == 0:
            OEThrow.Fatal("No fragments selected from %s." % inmol.GetTitle())

    if placeHyd:
        OEThrow.Verbose("Adding hydrogens to complex.")

        hopt = OEPlaceHydrogensOptions()
        hopt.SetAltsMustBeCompatible(compareAlts)
        hopt.SetStandardizeBondLen(standardize)
        hopt.SetWaterProcessing(waterProcess)
        hopt.SetBadClashOverlapDistance(badclash)
        hopt.SetFlipBiasScale(flipbias)
        hopt.SetMaxSubstateCutoff(maxStates)

        if verbose:
            details = OEPlaceHydrogensDetails()
            if not OEPlaceHydrogens(outmol, details, hopt):
                OEThrow.Fatal("Unable to place hydrogens and get details on %s." % inmol.GetTitle())
            OEThrow.Verbose(details.Describe())
        else:
            if not OEPlaceHydrogens(outmol, hopt):
                OEThrow.Fatal("Unable to place hydrogens on %s." % inmol.GetTitle())

    if addCharges:
        if not OEAssignCharges(outmol, OEMolComplexCharges()):
            OEThrow.Warning("Unable to assign mol complex charges to %s." % inmol.GetTitle())

    if addRadii:
        if not OEAssignRadii(outmol, OERadiiType_BondiHVdw, OERadiiType_HonigIonicCavity_Robust):
            OEThrow.Fatal("Unable to assign BondiH radii to %s." % inmol.GetTitle())

    oms1 = oemolostream()
    protFile = itf.GetString("-protout")
    if not oms1.open(protFile):
        OEThrow.Fatal("Unable to open %s for writing." % protFile)

    if splitlig:
        OEThrow.Verbose("Splitting ligand from complex.")

        frags = OEAtomBondSetVector()
        if not OEGetMolComplexFragments(frags, outmol, sopt):
            OEThrow.Fatal("Unable to fragment complex from %s." % inmol.GetTitle())

        lfilter = sopt.GetLigandFilter()
        wfilter = sopt.GetWaterFilter()
        lwfilter = OEOrRoleSet(lfilter, wfilter)

        protComplex = OEGraphMol()
        if not OECombineMolComplexFragments(protComplex, frags, sopt, OENotRoleSet(lwfilter)):
            OEThrow.Fatal("Unable to collect complex from %s." % inmol.GetTitle())

        if protComplex.NumAtoms() == 0:
            OEThrow.Warning("No complex identified in %s." % inmol.GetTitle())
        else:
            OEWriteMolecule(oms1, protComplex)

        lig = OEGraphMol()
        if not OECombineMolComplexFragments(lig, frags, sopt, lfilter):
            OEThrow.Fatal("Unable to collect ligand from %s." % inmol.GetTitle())

        if lig.NumAtoms() == 0:
            OEThrow.Warning("No ligand identified in %s." % inmol.GetTitle())
        else:
            oms2 = oemolostream()
            if splitlig:
                ligFile = itf.GetString("-ligout")
                if not oms2.open(ligFile):
                    OEThrow.Fatal("Unable to open %s for writing." % ligFile)

            OEThrow.Verbose("Ligand: %s" % lig.GetTitle())
            OEWriteMolecule(oms2, lig)
            oms2.close()
    else:
        OEWriteMolecule(oms1, outmol)

    oms1.close()

#############################################################################
# INTERFACE
#############################################################################

InterfaceData = '''
!BRIEF proteinprep.py [-options] <inmol> [<outcplx> [<outlig>]]

!CATEGORY "input/output options :" 1
   !PARAMETER -in 1
      !ALIAS -i
      !TYPE string
      !BRIEF Input molecule filename (must have 3D coordinates)
      !SIMPLE true
      !REQUIRED true
      !KEYLESS 1
   !END

   !PARAMETER -protout 2
      !ALIAS -p
      !TYPE string
      !DEFAULT proteinprep.oeb.gz
      !BRIEF Output protein filename
      !SIMPLE true
      !REQUIRED false
      !KEYLESS 2
   !END

   !PARAMETER -ligout 3
      !ALIAS -l
      !TYPE string
      !BRIEF Output ligand filename
      !SIMPLE true
      !REQUIRED false
      !KEYLESS 3
   !END
!END

!CATEGORY "Calculation options :" 2
    !PARAMETER -alts 1
       !TYPE string
       !LEGAL_VALUE occupancy
       !LEGAL_VALUE a
       !LEGAL_VALUE ignore
       !LEGAL_VALUE compare
       !DEFAULT occupancy
       !BRIEF Alternate location atom handling (affects atom:atom interactions)
       !SIMPLE true
       !REQUIRED false
       !DETAIL
         occupancy - keep just the highest average occupancy for each alt group
         a - keep only loc code A (and blank)
         ignore - assume alts already selected appropriately
         compare - keep all alts but only interact if same loc code (or blank)
    !END

    !PARAMETER -placehydrogens 2
       !TYPE bool
       !DEFAULT true
       !BRIEF If false, hydrogens will not be added
       !SIMPLE true
       !REQUIRED false
    !END

    !PARAMETER -charge 3
       !TYPE bool
       !DEFAULT true
       !BRIEF If false, OEMolComplexCharges will not be added
       !SIMPLE true
       !REQUIRED false
    !END

    !PARAMETER -radii 4
       !TYPE bool
       !DEFAULT true
       !BRIEF If false, BondiH radii will not be added
       !SIMPLE true
       !REQUIRED false
    !END

    !PARAMETER -waterprocessing 5
       !TYPE string
       !LEGAL_VALUE ignore
       !LEGAL_VALUE focused
       !LEGAL_VALUE fullsearch
       !DEFAULT fullsearch
       !BRIEF How waters are processed
       !SIMPLE true
       !REQUIRED false
       !DETAIL
         ignore - leave water hydrogens in a random orientation
         focused - search orientations based on neighboring polar groups
         fullsearch - do an extensive search of water orientations
    !END

    !PARAMETER -standardizehyd 6
       !ALIAS -stdhyd
       !TYPE bool
       !DEFAULT true
       !BRIEF If false, bonds for hydrogens are not adjusted to standard lengths
       !SIMPLE false
       !REQUIRED false
    !END

    !PARAMETER -clashcutoff 7
        !TYPE double
        !DEFAULT 0.4
        !BRIEF Van der Waals overlap (in Angstroms) defined to be a bad clash
        !SIMPLE false
        !REQUIRED false
    !END

    !PARAMETER -flipbias 8
        !TYPE double
        !DEFAULT 1.0
        !BRIEF Scale factor for the bias against flipping sidechains such as HIS
        !SIMPLE false
        !REQUIRED false
    !END

    !PARAMETER -maxsubstates 9
        !TYPE double
        !DEFAULT 1.0e8
        !BRIEF Maximum number of substates in a single step of hydrogen placement optimization
        !SIMPLE false
        !REQUIRED false
    !END

    !PARAMETER -newhyd 10
       !TYPE bool
       !DEFAULT false
       !BRIEF If true, hydrogens will be removed as a first step
       !SIMPLE true
       !REQUIRED false
    !END
!END

!CATEGORY "Display options :" 3
   !PARAMETER -verbose 1
      !ALIAS -v
      !TYPE bool
      !DEFAULT false
      !BRIEF Display more information about the process
      !SIMPLE true
      !REQUIRED false
   !END
!END
'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))
