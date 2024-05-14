#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <cmath>
using namespace std;
 
// Global constants 
#define size 100  // size of 2D lattice 
#define N 1e9  // number of steps
#define e_ff -2.5   // interaction energy 
#define e_sf -5.0   // solid-fluid adsorption potential 


// Functions declaration:
void writePositions( long array[size][size], long size_lat, std::string file_name ) {
    /*
        Writes coordinates of particles to a specified file *file_name*
    */

    ofstream myfile;
    myfile.open(file_name);
    for (int i = 0; i < size_lat; ++i) {
        for (int j = 0; j < size_lat; ++j) {
            if ( array[i][j] == 1 ){
              myfile << i << " " << j << "\n";
            }
        }
    }         
 
    myfile.close();
}


void writeEnergies( double eff_loc, double esf_loc, long N_loc, std::string file_name ) {
    /*
        Writes Eff and Esf energies to a specified file *file_name*
        Variables: 
            eff_loc - total eff energy
            esf_loc - total esf energy
            N - "time" step
    */

    ofstream myfile;
    myfile.open(file_name);
    
    myfile << eff_loc << " " << esf_loc << " " << N_loc << "\n";
    
    myfile.close();
}


double calcEloc( long array[size][size], long size_lat, long i, long j  ) {
    /*
        Calculates local energy for i,j coordinates
    */
    double eloc = 0.0;  
    int i_u, i_b, j_l, j_r; 

    j_r = j + 1;  // right particle 
    j_l = j - 1;  // left particle
    i_u = i + 1;  // upper particle
    i_b = i - 1;  // bottom particle
    if ( j_r >= size ) { j_r = 0; }
    if ( j_l < 0 ) { j_l = size-1; }
    if ( array[i][j_r] == 1 ) { eloc += e_ff; }
    if ( array[i][j_l] == 1 ) { eloc += e_ff; }
    if ( i_u >= size ) { eloc += e_sf; }
    else if ( i_b < 0 ) { eloc += e_sf; }
    else {
        if ( array[i_u][j] == 1 ) { eloc += e_ff; }
        if ( array[i_b][j] == 1 ) { eloc += e_ff; }
    }

    return eloc;
}


double calcEff( long array[size][size], long size_lat ) {
    /*
        Calculates a configuration energy
    */
    double eff = 0.0;
    double eff_cent = 0.0; 

    long i_u, i_b, j_l, j_r;  

    for (int i = 0; i < size_lat; ++i) {
        for (int j = 0; j < size_lat; ++j) {
            if (array[i][j] == 1){
                j_r = j + 1;  // right particle 
                j_l = j - 1;  // left particle
                i_u = i + 1;  // upper particle
                i_b = i - 1;  // bottom particle

                if ( j_r >= size ) { j_r = 0; }
                if ( j_l < 0 ) { j_l = size-1; }
                if ( array[i][j_r] == 1 ) { eff += e_ff/2.0 ; }
                if ( array[i][j_l] == 1 ) { eff += e_ff/2.0 ; }
                if ( i_u < size ) {
                    if ( array[i_u][j] == 1 ) { eff += e_ff/2.0 ; }
                } 
                if ( i_b >= 0 ) {
                    if ( array[i_b][j] == 1 ) { eff += e_ff/2.0 ; }
                }
            }
        }
    }  

    eff /= size*size*0.4*0.4;
    return eff;
}


double calcEsf( long array[size][size], long size_lat ) {
    /*
        Calculates a configuration energy
    */
    double esf = 0.0;

    for (int i = 0; i < size_lat; ++i) {
        for (int j = 0; j < size_lat; ++j) {
            if ( (i == 0) || (i == size-1) ) {
            if ( array[i][j] == 1) { esf += e_sf; } 
            } 
        }
    }  

    esf /= size*size*0.4*0.4;
    return esf;
}


int main () {
    // Declaration of variables 
    // Files, where print positions and energies
    std::string file_positions = "positions_all_";
    file_positions += to_string(e_sf);
    file_positions += ".txt";
    std::string file_energies = "energies_all_";
    file_energies += to_string(e_sf);
    file_energies += ".txt";

    long size_droplet = int(size * 0.2);  // droplet half-size equal to 20% of lattice's size
    long sites[size][size] = {{0}};   // initialize an empty 2D lattice 
    long filled_i[4*size_droplet*size_droplet] = {0};   // i coords filled positions particles 
    long filled_j[4*size_droplet*size_droplet] = {0};   // j coords filled positions particles 
    long empty_i[size*size-4*size_droplet*size_droplet] = {0};   // i coords empty positions  
    long empty_j[size*size-4*size_droplet*size_droplet] = {0};   // j coords empty positions  
    
    // Random number generator 
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist_filled(0, 4*size_droplet*size_droplet-1); // distribution in range [1, 6]
    std::uniform_int_distribution<std::mt19937::result_type> dist_empty(0, size*size-4*size_droplet*size_droplet-1); // distribution in range [1, 6]
    std::uniform_real_distribution<double> unif(0, 1);

    float ksi; 
    double eff, esf, eloc; 
    double en_old, en_new, en_del;  
    long count_filled = 0;  
    long count_empty = 0; 
    long ran_filled, ran_empty; 
    long i_fill, j_fill, i_emp, j_emp, swap;  

    // Init of a square with the size equal to 40% of 2D lattice at the center of left size  
    for (int i = size-2*size_droplet; i < size; ++i) {
        for (int j = int(size / 2)-size_droplet; j < int(size / 2)+size_droplet; ++j) {
            sites[i][j] = 1;
        }
    }        

    // Creating two list with filled and empty positions
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if ( sites[i][j] == 1 ) {
                filled_i[count_filled] = i;
                filled_j[count_filled] = j ;
                count_filled += 1;
            } else {
                empty_i[count_empty] = i;
                empty_j[count_empty] = j;
                count_empty += 1;
            }
        }
    }             

    eff = calcEff(sites, size);
    esf = calcEsf(sites, size);
    writeEnergies(eff, esf, 10, file_energies);

    long en_steps = long(N / 100);  // 100 timesteps for plots with energies 

    for ( int i = 0; i < N; ++i) {
        // Selection of random indices for arrays with filled and empty positions 

        ran_filled = dist_filled(rng);
        ran_empty = dist_empty(rng);

        i_fill = filled_i[ran_filled];
        j_fill = filled_j[ran_filled];
        i_emp = empty_i[ran_empty];
        j_emp = empty_j[ran_empty];

        en_old = calcEloc(sites, size, i_fill, j_fill);
        sites[i_emp][j_emp] = 1;    // swap of filled and 
        sites[i_fill][j_fill] = 0;  // empty positions
        
        en_new = calcEloc(sites, size, i_emp, j_emp);
        en_del = en_new - en_old;

        filled_i[ran_filled] = i_emp;
        filled_j[ran_filled] = j_emp;
        empty_i[ran_empty] = i_fill;
        empty_j[ran_empty] = j_fill;
           
        if ( en_del > 0 ) {
            ksi = unif(rng);
            if ( ksi >= exp(-en_del) ) { 
                sites[i_emp][j_emp] = 0;    // swap of filled and 
                sites[i_fill][j_fill] = 1;  // empty positions
                filled_i[ran_filled] = i_fill;   
                filled_j[ran_filled] = j_fill;   
                empty_i[ran_empty] = i_emp;
                empty_j[ran_empty] = j_emp;
            }
        }
    }

    writePositions(sites, size, file_positions);

    /*
    eff = calcEff(sites, size);
    esf = calcEsf(sites, size);
    writeEnergies(eff, esf, 10, file_energies);
    */

    long count_particles = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if ( sites[i][j] == 1 ) { count_particles += 1; }
        }
    }             

    printf("%lu \n", count_particles);

   return 0;
}
