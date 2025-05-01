// C++ program for Ex 4A: Extract the residue sequence from an sequence file
/*
Student name    : Bente Moebius
Username	    : bmoebius
Date		    : 23.11.2024
Task 		    : A
Time required	: 4 hours
Program Description:
This programm extracts information from an atm file about the structure of a protein. The user has to put in the
code for the protein and after that the number of atoms, residues as well as the type and numbers of residues are
extracted to separate arrays. When this is done succesfully then the residues are extracted from the arrays and it is
checked, that the sequence does not have any residues that get named two times in a row. This new list is then put out
into a new SEQ file and stored under the code for the protein.
Short report:
A:
P.      N. o. A.    N. o M.     Location and Nature
2ntg    1291        1           T115C
2ou8    1292        1           T115C
2ou9    1285        2           T115C, R119A
2q9e    1289        2           S44C, N55A
1ssy    1291        3           G28A, I29A, G30A
1ssw    1273        4           Y24A, Y25A, T26A, I27A
1t8f    1269        6           R14A, K16A, I17A, L19A, T21A, E22A

B:
I could check if the file is not empty before reading it out.
Futher I could check that the needed directory for seq already exists, or that it is empty etc.
I also could add an extra clause for the first residue, so that it is safely formatted.
I could also make the variable max_atoms accessible to guarantee that the file that is read is not bigger than the
arrays.
*/
// INCLUDES
#include <iostream>
#include  <iomanip>
#include <fstream>
#include  <string>

// NAMESPACE
using namespace std;

// GLOBAL CONSTANTS TO DEFINE ARRAY SIZES & FILE PATHS
// @ already written for you below: you can use these throughout the code
const int max_atoms        = 5000;      // max number of atoms that can be stored
const int max_residues     = 500;       // max number of residues that can be stored
const string atm_directory = "../atm/"; // directory where to search for ATM-files
const string atm_extension = ".atm";    // filename extension of ATM-files
const string seq_directory = "../seq/"; // directory where to write SEQ-files
const string seq_extension = ".seq";    // filename extension of SEQ-files

// PROTOTYPES
// read the PDB-entry code from the standard input;
bool GetPdbCode(string& pdb_code);
// Reat the ATM File and extract the needed information
int  ReadAtmFile(string& pdb_code, string atm_directory, string atm_extension, int atom_numbers[], string atom_types[], int residue_numbers[], string residue_types[]);
// fill the array with the residue sequence based on the atm file
int ListResidues(int num_atoms, string residue_types[], int residue_numbers[], string residue_list[]);
//  WRITE THE RESIDUE SEQUENCE TO A SEPARATE FILE
bool WriteResidues(string& pdb_code, string seq_directory, string seq_extension, int num_residues, string residue_list[]);

// MAIN
int main() {

  // VARIABLE DECLARATIONS
  string pdb_code;                   // PDB-entry code
  int    num_atoms;                  // number of atoms in the protein structure
  int    num_residues;               // number of residues in the protein structure
  int    atom_numbers[max_atoms];    // int array with the atom numbers of atoms
  string atom_types[max_atoms];      // string array with the atom types of atoms
  int    residue_numbers[max_atoms]; // int array with the residue numbers of atoms
  string residue_types[max_atoms];   // string array with the residue types of atoms
  string residue_list[max_residues]; // string array with residue types of residues

  // CODE FOR MAIN
    // GET THE INPUT FROM USER AND END PROGRAMM IF PDB CODE IS NOT UP TO STANDARD
    if (!GetPdbCode(pdb_code)) {
        cout << "# PDB code not found" << endl;
        return 1;
    }
    cout << "# Got PDB code: " << pdb_code << endl;
    // EXTRACT DATA FROM ATM FILE AND RETURN NUMBER OF ATOMS
    num_atoms = ReadAtmFile(pdb_code, atm_directory, atm_extension, atom_numbers, atom_types, residue_numbers, residue_types);
    // END PROGRAMM IF AN ERROR OCCURED
    if (!num_atoms) {
        cout << "# An Error occcured while reading the file." << endl;
        return 1;
    }
    cout << "# Number of atoms: " << num_atoms << endl;
    // GET THE DIFFERENT AMINO ACIDS AND LIST IT INTO A NEW ARRAY
    // RETURN NUMBER OF AMINOC ACIDS / RESIDUES
    num_residues = ListResidues(num_atoms, residue_types, residue_numbers, residue_list);
    // END PROGRAMM IF AN ERROR ORCCURED
    if (!num_residues) {
        cout << "# An Error occcured while calculating the residues." << endl;
        return 1;
    }
    cout << "# Number of residues in the sequence: " << num_residues << endl;
    // WRITE THE EXTRACTED SEQUENCE OF RESIDUES INTO A NEW SEQ-FILE
    if (!WriteResidues(pdb_code, seq_directory, seq_extension, num_residues, residue_list)) {
        cout << "# An Error occcured while writing the file." << endl;
        return 1;
    }
    cout << "# Residue sequence written to file." << endl;
  return 0;

}

// FUNCTION DEFINITIONS FOR THE 4 REQUIRED FUNCTIONS
// GET THE PDB CODE FROM THE USER INPUT
bool GetPdbCode(string& pdb_code) {
  cout << "# Enter the PDB code (use lower case letters): " << endl;
  cin >> pdb_code;
  // FOR LOOP TO CHECK EACH CHARACTER IF IT IS THE CORRECT FORMAT
  for (int i = 0; i < pdb_code.length(); i++) {
     if (!(islower(pdb_code[i]) || isdigit(pdb_code[i]))) { return false; }
  }
  // RETURN TRUE IF EVERYTHIN IS UP TO STANDARD
  return pdb_code.size()==4;
}

// READ THE ATM FILE AND SAFE THE ATOM NUMBERS, TYPES AND RESIDUE TYPES TO AN ARRAY
int  ReadAtmFile(string& pdb_code, string atm_directory, string atm_extension, int atom_numbers[], string atom_types[], int residue_numbers[], string residue_types[]) {
    string filename = atm_directory + pdb_code + atm_extension;
    int current_atom = 0;
    // open the file
    ifstream inp(filename);
    // check if file exists and is readable
    if (!inp.is_open()) {
        cout << "Error: Unable to open file: " << filename << endl;
        return 0; // Return 0 to indicate failure
    }
    // begin reading from file
    cout << "# Reading data from ATM file: " << filename << endl;
    string word;
    //
    while (!inp.eof() && current_atom < max_atoms) {
        // read everythin valuable from the input
        inp >> word >> atom_numbers[current_atom] >> atom_types[current_atom] >> residue_numbers[current_atom] >> residue_types[current_atom];
        // skip the coordinates at the end
        inp.ignore(max_atoms, '\n');
        // increase the counter of the atoms
        ++current_atom;
    }
    inp.close();
    // return the number of atoms in the protein
    return current_atom;
}

// EXTRACT THE DIFFERENT RESIDUES TO A SEPARATE LIST
int ListResidues(int num_atoms, string residue_types[], int residue_numbers[], string residue_list[]) {
    cout << "# Counting Residues... " << endl;
    int current_residue = 0;
    // first expression to ensure that nothing goes wrong
    residue_list[current_residue] = residue_types[0];
    current_residue++;
    // for loop to go through all of the residues
    for (int i = 1; i < num_atoms; i++) {
        // check if the residue in front is the same as the one looked at right now
        if (residue_numbers[i] != residue_numbers[i-1]) {
            residue_list[current_residue] = residue_types[i];
            // increase counter of residues
            current_residue++;
        }
    }
    // return the number of residues
    return current_residue;
}

// WRITE THE SEQUENCE INTO A SEPARATE SEQ FILE
bool WriteResidues(string& pdb_code, string seq_directory, string seq_extension, int num_residues, string residue_list[]) {
    string filename = seq_directory + pdb_code + seq_extension;
    ofstream out(filename);
    // check if the output file is really open
    if (!out.is_open()) {
        cout << "Error: Unable to open file: " << filename << endl;
        return false;
    }
    cout << "# Writing residues to seq-file: " << filename << endl;
    // for loop to list each residue
    for (int i = 0; i < num_residues-1; i++) {
        out << i+1 << " " << residue_list[i] << endl;
    }
    out.close();
    return true;
}
