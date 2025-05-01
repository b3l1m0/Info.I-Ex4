// C++ program for Ex 4C: Compare 8 structures to each other and calculate each number of mutations and rmsd
/*
Student name    : Bente Moebius
Username	    : bmoebius
Date		    : 23.11.2024
Task 		    : C
Time required	: 2 hours
Program Description:
This programm takes the input from the user, which consists of 8 codes for 8 proteins. After that the number of atoms,
residues as well as the type and numbers of residues and the coordinates of the atoms of each of these proteins are
extracted from their atm file to separate two-dimensional arrays.
When this is done succesfully, then for all of the proteins, the residues are extracted from the arrays and it is checked,
that the sequences do not have any residues that get named twice in a row. This new list is then put out into a
new SEQ file for each protein and stored under the code for the protein.
After that, each of the proteins is compared to each other protein to find out the number of mutations between these two,
which is stored accordingly in the variable and also in the storage array for each pair. This is also done for the codes
of each pair, which is stored in separate arrays.
Further the distance of each of the atoms at the same number of each protein is calculated and added up to a vector,
which descripes the difference in structures mathematically (rmsd). This is done for each individual pair and is also
stored in a new array for each of these values.
The arrays are getting sorted using a simple selection sort of the array of the rmsd values. This array, as well as the
others containing the pairs and the number of mutations is then sorted at the end from lowest rmsd value to highest
rmsd value.
At the end the codes for each pair of proteins, the number of mutations between these proteins as well as their rmsd
value is put out.
Short report:
A.
protein 1 | protein 2 | mutations | RMSD
2oe4 | 2oe4 | 0 | 0
2ou8 | 2ou8 | 0 | 0
2ntg | 2ntg | 0 | 0
2ou9 | 2ou9 | 0 | 0
1ssw | 1ssw | 0 | 0
2q9e | 2q9e | 0 | 0
1ssy | 1ssy | 0 | 0
1t8f | 1t8f | 0 | 0
2oe4 | 2ou8 | 1 | 0.197609
2oe4 | 2ntg | 1 | 0.216464
2ntg | 2ou9 | 1 | 0.220935
2oe4 | 1ssw | 4 | 0.2302
2oe4 | 2ou9 | 2 | 0.237204
2ou8 | 2ou9 | 1 | 0.241263
2ou8 | 1ssw | 5 | 0.281267
2ntg | 2ou8 | 0 | 0.293572
2ntg | 1ssw | 5 | 0.303237
2ou9 | 1ssw | 6 | 0.307217
2oe4 | 2q9e | 2 | 0.442787
2q9e | 1ssw | 6 | 0.484165
2ntg | 2q9e | 3 | 0.506763
2ou8 | 2q9e | 3 | 0.507254
2ou9 | 2q9e | 4 | 0.538306
2ntg | 1t8f | 7 | 1.11802
2ou9 | 1t8f | 8 | 1.13076
2ou8 | 1t8f | 7 | 1.14394
2oe4 | 1t8f | 6 | 1.15409
1ssw | 1t8f | 10 | 1.17375
2q9e | 1t8f | 8 | 1.19238
2oe4 | 1ssy | 3 | 1.26215
2q9e | 1ssy | 5 | 1.30078
2ou8 | 1ssy | 4 | 1.26486
2ou9 | 1ssy | 5 | 1.30503
1ssy | 1ssw | 7 | 1.31654
2ntg | 1ssy | 4 | 1.36257
1ssy | 1t8f | 9 | 2.01851

B.
It is visislbe that there is a slight increase in RMSD as the number of mutations increase. But again, it is
also visible, that not all mutations impact the RMSD the same, as some differences have many mutations, but a
relatively small RMSD value.
Also the location of the mutation again, impacts the RMSD value greatly. This is probably because this impacts
the location of the rest of the atoms greatly, even if the rest is structured regularly.
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
const int num_structures = 8;      // number of structures analysed
const int max_combi = num_structures*(num_structures-1)/2+num_structures; // for initialization of arrays, the needed data storage

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
int CalculateMutations(int index_structure1, int index_structure2, int num_residues[], string residue_list[][max_residues]);
// CALCULATE THE BACKBONE RMSD BETWEEN REFERENCE AND TARGET
double CalculateRMSD(int index_structure1, int index_structure2, int num_atoms, double atom_coordinates[][max_atoms][num_coords]);
// Sort the arrays using selection sort
void SortData(int first_structure[], int second_structure[], int num_mutations[], double rmsd[]);
// function to swap the pointers of an array consisting of doubles, needed for the sort function
void SwapDouble(double *pointer1, double *pointer2);
// function to swap the pointers of an array consisting of integers
void SwapInt(int *pointer1, int *pointer2);

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
    int num_mutations[max_combi]; // number of mutations between two proteins
    double rmsd[max_combi]; //measure of the difference between two structures
    int first_structure[max_combi]; // store the first protein of the pair
    int second_structure[max_combi]; // store the second protein of the pair

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
    for (int protein = 1; protein < num_structures; protein++) {
        if (num_atoms[protein] != num_atoms[protein-1]) { return 1; }
    }

    // loop through all of the proteins
    int index = 0;
    for (int protein1 = 0; protein1 < num_structures; protein1++) {
        // second loop to loop through all of the proteins it still has to match with
        for (int protein2 = protein1; protein2 < num_structures; protein2++) {
            // superimposing the structures
            cout << "# Superimposing the structures: " << pdb_code[protein1] << " & " << pdb_code[protein2] << endl;
            StructureFit(num_atoms[0], atom_coordinates[protein1], atom_coordinates[protein2]);

            // calculate number of mutations;
            cout << "# Calculating the Number of Mutations..." << endl;
            num_mutations[index] = CalculateMutations(protein1, protein2, num_residues, residue_list);
            // calculate the backbone
            cout << "# Calculating the rmsd" << endl;
            rmsd[index] = CalculateRMSD(protein1, protein2, num_atoms[0], atom_coordinates);
            // store in arrays
            first_structure[index] = protein1;
            second_structure[index] = protein2;
            ++index;
        }
    }

    //Sort data in the arrays
    SortData(first_structure, second_structure, num_mutations, rmsd);

    // put out the data from each of the protein pairs
    cout << "# protein 1 | protein 2 | number of mutations | RMSD" << endl;
    for (int n = 0; n < max_combi; n++) {
        cout << pdb_code[first_structure[n]] << " | ";
        cout << pdb_code[second_structure[n]] << " | ";
        cout << num_mutations[n] << " | ";
        cout << rmsd[n] << endl;
    }
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
    // different first clause to not get problems
    residue_list[current_residue] = residue_types[0];
    current_residue++;
    // for loop to compare each of the residues to the one in front to not get doublings
    for (int i = 1; i < num_atoms; i++) {
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
    // check that the file is open
    if (!out.is_open()) {
        cout << "Error: Unable to open file: " << filename << endl;
        return false;
    }
    // write each of the residues to a seq file
    cout << "# Writing residues to seq-file: " << filename << endl;
    for (int i = 0; i < num_residues; i++) {
        out << i+1 << " " << residue_list[i] << endl;
    }
    out.close();
    return true;
}

// calculation of the number of mutations between two proteins
int CalculateMutations(int index_structure1, int index_structure2, int num_residues[], string residue_list[][max_residues]) {
    int num_mutations = 0;
    for (int index = 0; index < num_residues[0]; index++) {
        // if the structures do not match, increase counter
        if (!(residue_list[index_structure1][index] == residue_list[index_structure2][index])) {
            num_mutations++;
        }
    }
    // return number of mutations
    return num_mutations;
}

// calculate the root mean square deviation between two samples
double CalculateRMSD(int index_structure1, int index_structure2, int num_atoms, double atom_coordinates[][max_atoms][num_coords]) {
    // vector for the difference
    double sum = 0.0;
    double diff = 0.0;
    for (int counter = 0; counter < num_atoms; counter++) {
        for (int coor_counter = 0; coor_counter < num_coords; coor_counter++) {
            // multiply each of the differences between the coordinates with itself
            diff = atom_coordinates[index_structure2][counter][coor_counter] - atom_coordinates[index_structure1][counter][coor_counter];
            sum += diff * diff;
        }
    }
    // return the rmsd value
    return sqrt(sum/double(num_atoms));
}

// sort the data of the arrays from lowest backbone to highest
void SortData(int first_structure[], int second_structure[], int num_mutations[], double rmsd[]) {
    // variable to store the last index which is sorted
    int help;
    for (int i = 0; i < max_combi; i++) {
        help = i;
        for (int j = i+1; j < max_combi; j++) {
            if (rmsd[j] < rmsd[help])
                help = j;
            // swap the two pointers as this also swaps the sequence of the values
            SwapDouble(&rmsd[help], &rmsd[i]);
            SwapInt(&first_structure[help], &first_structure[i]);
            SwapInt(&second_structure[help], &second_structure[i]);
            SwapInt(&num_mutations[help], &num_mutations[i]);

        }
    }
}

// function to swap the pointers of two doubles
void SwapDouble(double *pointer1, double *pointer2) {
    double help = *pointer1;
    *pointer1 = *pointer2;
    *pointer2 = help;
}

// function to swap the pointers of two ints
void SwapInt(int *pointer1, int *pointer2) {
    int help = *pointer1;
    *pointer1 = *pointer2;
    *pointer2 = help;
}