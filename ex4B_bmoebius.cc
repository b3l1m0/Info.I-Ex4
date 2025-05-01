// C++ program for Ex 4B: calculate structural differences of a mutation to the normal type
/*
Student name    : Bente Moebius
Username	    : bmoebius
Date		    : 23.11.2024
Task 		    : B
Time required	: 3 hours
Program Description:
This programm takes the input from the user, which consists of 2 codes for the protein. After that the number of atoms,
residues as well as the type and numbers of residues and the coordinates of the atoms of each of these proteins are
extracted from their atm file to separate arrays.
When this is done succesfully, then for the two proteins, the residues are extracted from the arrays and it is checked,
that the sequences do not have any residues that get named two times in a row. This new list is then put out into a
new SEQ file for each protein and stored under the code for the protein.
After that, the list of the two proteins are compared to find out the number of mutations between these two, which is
stored accordingly in the variable.
Further the distance of each of the atoms at the same number of each protein is calculated and added up to a vector,
which descripes the difference in structures mathematically (rmsd).
At the end, the number of mutations between the two proteins as well as their rmsd value is put out.
Short report:
A.
        Num. Mut.     RMSD
2oe4:       0       7.37996e-15
2ntg:       1       0.216464
2ou8:       1       0.197609
2ou9:       2       0.237204
2q9e:       2       0.442787
1ssy:       3       1.26215
1ssw:       4       0.2302
1t8f:       6       1.15409

B:
There is a small correlation, that with larger structural differences, also the RMSD is pretty large.
And i think that it also plays a role, when in the protein the mutation takes place. In 2ou9 and 2q9e
it is visible, that the RMSD of 2q9e is larger, which could correlate to the fact that in this protein
the mutation took place at the place 44, which is pretty early in the sequence whereas in 2ou9 the
mutation takes place at 119, which is in comparison much later.
Further it seems that some amino acids are more similar than others regarding the structure of the proten.
As the protein 1ssw has more mutations than 1ssy but the RMSD is lower, which could maybe stem from the fact
that there are different amino acids, which mutated.
I think that we can expect some correlation of the strcture to the sequence in regard to the placing of the
mutation and also to some extent to the number of mutations.
*/
// INCLUDES
#include <iostream>
#include  <iomanip>
#include <fstream>
#include  <string>
#include <cmath>
#include "fit.h"

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
const int num_coords = 3;          // number of coordinates
const int num_structures = 2;      // number of structures analysed

// PROTOTYPES
// read the PDB-entry code from the standard input;
bool GetPdbCode(string pdb_code[]);
// Reat the ATM File and extract the needed information
int  ReadAtmFile(string pdb_code, string atm_directory, string atm_extension, int atom_numbers[], string atom_types[], int residue_numbers[], string residue_types[], double atom_coordinates[][num_coords]);
// fill the array with the residue sequence based on the atm file
int ListResidues(int num_atoms, string residue_types[], int residue_numbers[], string residue_list[]);
//  WRITE THE RESIDUE SEQUENCE TO A SEPARATE FILE
bool WriteResidues(string pdb_code, string seq_directory, string seq_extension, int num_residues, string residue_list[]);
// CALCULATE THE NUMBER OF MUTATIONS
int CalculateMutations(int num_residues[], string residue_list[][max_residues]);
// CALCULATE THE BACKBONE RMSD BETWEEN REFERENCE AND TARGET
double CalculateRMSD(int num_atoms, double atom_coordinates[][max_atoms][num_coords]);

// MAIN
int main() {
    // VARIABLE DECLARATIONs
    // also changed to arrays of num_structures
    string pdb_code[num_structures];                   // PDB-entry code
    int    num_atoms[num_structures];                  // number of atoms in the protein structure
    int    num_residues[num_structures];               // number of residues in the protein structure
    int    atom_numbers[num_structures][max_atoms];    // int array with the atom numbers of atoms
    string atom_types[num_structures][max_atoms];      // string array with the atom types of atoms
    int    residue_numbers[num_structures][max_atoms]; // int array with the residue numbers of atoms
    string residue_types[num_structures][max_atoms];   // string array with the residue types of atoms
    string residue_list[num_structures][max_residues]; // string array with residue types of residues
    double atom_coordinates[num_structures][max_atoms][num_coords]; // double array to store the coordinates of each atom
    int num_mutations; // number of mutations between two proteins
    double rmsd; //measure of the difference between two structures

    // CODE FOR MAIN
    // END PROGRAMM IF PDB CODE IS NOT UP TO STANDARD
    if (!GetPdbCode(pdb_code)) {
        cout << "# PDB code not found" << endl;
        return 1;
    }
    cout << "# Got PDB code: " << pdb_code << endl;

    // loop to loop all of the tasks for each structure
    for (int n = 0; n < num_structures; n++) {

        cout << "## NOW PERFOMING TASKS FOR STRUCTURE: " << pdb_code[n] << endl;

        // EXTRACT DATA FROM ATM FILE AND RETURN NUMBER OF ATOMS
        num_atoms[n] = ReadAtmFile(pdb_code[n], atm_directory, atm_extension, atom_numbers[n], atom_types[n], residue_numbers[n], residue_types[n], atom_coordinates[n]);
        // END PROGRAMM IF AN ERROR OCCURED
        if (!num_atoms[n]) {
            cout << "# An Error occcured while reading the file." << endl;
            return 1;
        }
        cout << "# Number of atoms: " << num_atoms[n] << endl;

        // LIST THE RESIDUES AND FIND THE NUMBER
        num_residues[n] = ListResidues(num_atoms[n], residue_types[n], residue_numbers[n], residue_list[n]);
        // END PROGRAMM IF AN ERROR ORCCURED
        if (!num_residues[n]) {
            cout << "# An Error occcured while calculating the residues." << endl;
            return 1;
        }
        cout << "# Number of residues in the sequence: " << num_residues[n] << endl;

        // WRITE RESIDUES TO FILE
        if (!WriteResidues(pdb_code[n], seq_directory, seq_extension, num_residues[n], residue_list[n])) {
            cout << "# An Error occcured while writing the file." << endl;
            return 1;
        }
        cout << "# Residue sequence written to file." << endl;
    }
    // check that the different proteins have the same number of of atoms
    if (num_atoms[0] != num_atoms[1]) {
    	cout << "# The number of atoms in the different sequences is not the same." << endl;
    	return 1;
	}

    // superimpose the target structure onto the reference one
    cout << "# Superimposing the target structure on the reference structure." << endl;
    StructureFit(num_atoms[0], atom_coordinates[0], atom_coordinates[1]);

    // calculate number of mutations;
    cout << "# Calculating the Number of Mutations..." << endl;
    num_mutations = CalculateMutations(num_residues, residue_list);

    // calculate the backbone
    cout << "# Calculating the rmsd" << endl;
    rmsd = CalculateRMSD(num_atoms[0], atom_coordinates);

    // output for the user
    cout << "# Number of mutations between the two proteins: " << endl;
    cout << num_mutations << endl;
	cout << "# Root Mean Square Deviation (RMSD) of backbone atoms: " << endl;
    cout << rmsd << endl;

    return 0;

}

// FUNCTION DEFINITIONS FOR THE 4 REQUIRED FUNCTIONS
// GET THE PDB CODE FROM THE USER INPUT
bool GetPdbCode(string pdb_code[]) {
    // for loop to loop through the inputs
    for (int input = 0; input < num_structures; input++) {
        // get different pdb codes from user
        cout << "# Enter the PDB code (use lower case letters) for the " << input+1 << ". structure: " << endl;
        cin >> pdb_code[input];
        // check if the input is the needed length
        if (pdb_code[input].length() != 4) {
            return false;
        }
        // loop through all of the characters of the pdb code to check if the inout is valid
        for (int character = 0; character < 4; character++) {
            if (!(islower(pdb_code[input][character]) || isdigit(pdb_code[input][character]))) { return false; }
        }
    }
    return true;
}

// READ THE ATM FILE AND SAFE THE ATOM NUMBERS, TYPES AND RESIDUE TYPES TO AN ARRAY
int  ReadAtmFile(string pdb_code, string atm_directory, string atm_extension, int atom_numbers[], string atom_types[], int residue_numbers[], string residue_types[], double atom_coordinates[][num_coords]) {
    string filename = atm_directory + pdb_code + atm_extension;
    int current_atom = 0;
    // open the file
    ifstream inp(filename);
    // check if file exists and is readable
    if (!inp.is_open()) {
        cout << "Error: Unable to open file: " << filename << endl;
        return false; // Return 0 to indicate failure
    }
    // begin reading from file
    cout << "# Reading data from ATM file: " << filename << endl;
    string word;
    inp >> word;
    while (!inp.eof() && current_atom < max_atoms) {
        // Skip the coordinates at the end
        inp >> atom_numbers[current_atom] >> atom_types[current_atom] >> residue_numbers[current_atom] >> residue_types[current_atom] >> atom_coordinates[current_atom][0] >> atom_coordinates[current_atom][1] >> atom_coordinates[current_atom][2];
        // only record the atoms with CA N C or O
        if (atom_types[current_atom] == "N" || atom_types[current_atom] == "CA" || atom_types[current_atom] == "C" || atom_types[current_atom] == "O") { ++current_atom; }
        inp >> word;
    }
    inp.close();
    return current_atom;
}

// EXTRACT THE DIFFERENT RESIDUES TO A SEPARATE LIST
int ListResidues(int num_atoms, string residue_types[], int residue_numbers[], string residue_list[]) {
    cout << "# Counting Residues... " << endl;
    int current_residue = 0;
    // DIFFERENT CASE FOR THE FIRST COMPARISON
    residue_list[current_residue] = residue_types[0];
    current_residue++;
    // LOOP THROUGH ALL OF THE RESIDUES IN THE ARRAY
    for (int i = 1; i < num_atoms; i++) {
        // CHECK IF THE RESIDUE IN FRONT IS THE SAME AS THE ONE THE COUNTER IS ON RIGHT NOW
        if (residue_numbers[i] != residue_numbers[i-1]) {
            residue_list[current_residue] = residue_types[i];
            current_residue++;
        }
    }
    return current_residue;
}

// WRITE THE SEQUENCE INTO A SEPARATE SEQ FILE
bool WriteResidues(string pdb_code, string seq_directory, string seq_extension, int num_residues, string residue_list[]) {
    string filename = seq_directory + pdb_code + seq_extension;
    ofstream out(filename);
    // check if the output is open
    if (!out.is_open()) {
        cout << "Error: Unable to open file: " << filename << endl;
        return false;
    }
    cout << "# Writing residues to seq-file: " << filename << endl;
    // for loop to go through every residue and write it into the file
    for (int i = 0; i < num_residues; i++) {
        out << i+1 << " " << residue_list[i] << endl;
    }
    out.close();
    return true;
}

// calculation of the number of mutations between two proteins
int CalculateMutations(int num_residues[], string residue_list[][max_residues]) {
    int num_mutations = 0;
    // for loop to go through each of the residues
    for (int index = 0; index < num_residues[0]; index++) {
        // if the structures do not match, increase counter
        if (!(residue_list[0][index] == residue_list[1][index])) {
            num_mutations++;
        }
    }
    // return number of mutations
    return num_mutations;
}

// calculate the root mean square deviation between two samples
double CalculateRMSD(int num_atoms, double atom_coordinates[][max_atoms][num_coords]) {
    // vector for the difference
    double sum = 0.0;
    double diff = 0.0;
    for (int counter = 0; counter < num_atoms; counter++) {
        for (int coor_counter = 0; coor_counter < num_coords; coor_counter++) {
            // multiply each of the differences between the coordinates with itself
            diff = atom_coordinates[1][counter][coor_counter] - atom_coordinates[0][counter][coor_counter];
            sum += diff * diff;
        }
    }
    // return the total rmsd
    return sqrt(sum/double(num_atoms));
}