// Sample C++ program for Ex 4A: Extract the residue sequence from an
//                               ATM-file and write it to a SEQ-file

// INCLUDES
// @ add here the "#include" for the library functions you will need
//   in your program (probably at least: iostream, fstream, and string)

// NAMESPACE
// @ add here the usual "namespace" declaration

// GLOBAL CONSTANTS TO DEFINE ARRAY SIZES & FILE PATHS
// @ already written for you below: you can use these throughout the code
const int max_atoms        = 5000;      // max number of atoms that can be stored
const int max_residues     = 500;       // max number of residues that can be stored
const string atm_directory = "../atm/"; // directory where to search for ATM-files
const string atm_extension = ".atm";    // filename extension of ATM-files
const string seq_directory = "../seq/"; // directory where to write SEQ-files
const string seq_extension = ".seq";    // filename extension of SEQ-files

// PROTOTYPES
// @ this time, you have to write the prototypes of the functions yourself;
//   the functions you will need are:
//
//   - bool GetPdbCode(...):    read the PDB-entry code from the standard input;
//                              return true if all went fine, false if problem
//   - int  ReadAtmFile(...):   read the ATM-file and fill the appropriate arrays;
//                              return the number of atoms read, or zero if problem
//   - int  ListResidues(...):  fill the array with the residue sequence based on the 
//                              arrays from the ATM-file; return the number of
//                              residues found, or zero if problem 
//   - bool WriteResidues(...): write the residue sequence to the SEQ-file;
//                              return true if all went fine, false if problem

// MAIN

int main() {

  // VARIABLE DECLARATIONS
  // @ declare your variables; to avoid that everyone uses different variable names
  //   (which would make the correction difficult), we already defined for you the 
  //   main variables and arrays that your program will need...
  
  string pdb_code;                   // PDB-entry code
  int    num_atoms;                  // number of atoms in the protein structure
  int    num_residues;               // number of residues in the protein structure
  int    atom_numbers[max_atoms];    // int array with the atom numbers of atoms
  string atom_types[max_atoms];      // string array with the atom types of atoms
  int    residue_numbers[max_atoms]; // int array with the residue numbers of atoms
  string residue_types[max_atoms];   // string array with the residue types of atoms
  string residue_list[max_residues]; // string array with residue types of residues

  // CODE FOR MAIN
  // @ add here your code for main: call in turn GetPdbCode, 
  //   ReadAtmFile, ListResidues and WriteResidues

  return 0;

}

// FUNCTION DEFINITIONS FOR THE 4 REQUIRED FUNCTIONS
// @ add here the function definitions for GetPdbCode,
//   ReadAtmFile, ListResidues and WriteResidues

