#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <cmath>
using namespace std;
 
// Global constants 
#define size 30             // size of 2D lattice 
#define steps 50000000          // number of energies in a file 
#define gamma_steps 100     // 
#define e_ff -2.0           // interaction energy 
#define e_sf -2.0           // solid-fluid adsorption potential 
#define rho_vle 0.0185      // vapor-liquid equilibrium equilibrium density 
#define rho_steps 30     // number of rho_bulk between rho_min and rho_max
#define V size*size         // system's volume 


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

double calcEloc( long array[size][size], long size_lat, long i, long j  ) {
    /*
        Calculates local energy for i,j coordinates
        Only one contact layer y=1
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
    else {
        if ( array[i_u][j] == 1 ) { eloc += e_ff; }
    }
    if ( i_b > 0 ) { 
        if ( array[i_b][j] == 1 ) { eloc += e_ff; }
    }
    return eloc;
}

double calcGamma( long array[size][size], long size_lat ) {
    /*
        Calculates gamma 
        Variables: 
            array - list with positions  
            size_lat - size of lattice 
    */
    double gamma = 0.0; 
    
    for ( int j = 0; j < size_lat; ++j ) {
        for ( int i = size_lat-1; i > 0; --i ) {
            if (array[i][j] == 1) { gamma += 1; }
            else { break; }
        }
    }

    gamma /= size; 
    return gamma; 
}

void writeGamma( double list_1[gamma_steps+1], long list_2[gamma_steps+1], std::string file_name ) {
    /*
        Writes Eff and Esf energies to a specified file *file_name*
        Variables: 
            list_1 - list with gamma 
            list_2 - list with step
            file_name - name of file where the data is written 
    */

    ofstream myfile;
    myfile.open(file_name);
    
    for ( int i = 0; i < gamma_steps; ++i ) {
    myfile << list_1[i] << "\t" << list_2[i] << "\n";
    }

    myfile.close();
}


int main () {
    // Declaration of variables 
    // Files, where print positions and energies
    std::string file_positions = "positions_";
    file_positions += to_string(e_sf);
    std::string file_gamma = "gamma_";
    file_gamma += to_string(e_sf);

    std::string file_positions_cur;
    std::string file_gamma_cur;

    long sites[size][size] = {{0}};   // initialize an empty 2D lattice 
    long filled_i[V] = {0};   // i coords filled positions particles 
    long filled_j[V] = {0};   // j coords filled positions particles 
    
    // Random number generator 
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist_size(0, size-1); // distribution in range [0, size]
    std::uniform_real_distribution<double> unif(0, 1);

    // Variables for the MC simulation
    float ksi; 
    long i_in, j_in, i_rem, j_rem;
    double en_new, comp;  
    long count_filled = 0;  

    // Variables for the file with energies
    long N;          
    long idx_rem;           // index to remove 
    double rho_b;           // bulk density

    // Variables for adsorption values
    long  gamma_inter = long(steps / gamma_steps);    // number of gamma steps 
    long gamma_idx = 0; 
    long gamma_step[gamma_steps+3] = {0};
    double gamma_list[gamma_steps+3] = {0}; 
    double rho_min = 0.01;   // minimum rho
    double rho_max = 0.95;   // maximum rho
    double rho_step; 

    rho_step = (rho_max - rho_min) / (rho_steps-1); 
    for (int iter = 0; iter < rho_steps; ++iter ){
        // initial configuration for each rho 
        long sites[size][size] = {{0}};    
        filled_i[V] = {0};    
        filled_j[V] = {0};    
        gamma_idx = 0; 
        gamma_step[gamma_steps+3] = {0};
        gamma_list[gamma_steps+3] = {0}; 
   

        rho_b = rho_vle * (rho_min + iter*rho_step); 
        N = 0;
        for ( int k = 0; k < steps; ++k) {
            // Insertion of a particle 
            i_in = dist_size(rng);
            j_in = dist_size(rng);
            if ( sites[i_in][j_in] == 0 ) {
                // en_old = calcEloc(sites, size, i_in, j_in);
                sites[i_in][j_in] = 1;
                en_new = calcEloc(sites, size, i_in, j_in);
                ksi = unif(rng);
                comp = (V*rho_b)/(N+1)*exp(-en_new); 
                if ( ksi >= comp ) {
                    sites[i_in][j_in] = 0;
                } else {
                    filled_i[N] = i_in;
                    filled_j[N] = j_in;
                    N += 1; 
                } 
            }

            //Removal a particle from an occupied position 
            if ( N > 0 ) {
                // If there is at least one particle in the system
                idx_rem = long( unif(rng)*N );
                i_rem = filled_i[idx_rem];
                j_rem = filled_j[idx_rem];
                sites[i_rem][j_rem] = 0;
                en_new = calcEloc(sites, size, i_rem, j_rem);
                comp = (N)/(rho_b*V)*exp(en_new);
                ksi = unif(rng);
                if ( ksi < comp ) {
                    // Shift to indices left
                    for ( int j = idx_rem; j < N-1; ++j ) {
                        filled_i[j] = filled_i[j+1];
                        filled_j[j] = filled_j[j+1];
                    }
                    filled_i[N] = 0;
                    filled_j[N] = 0;
                    N -= 1;
                } else {
                    sites[i_rem][j_rem] = 1;
                }
            }

            // Data for file with adsorption coefficients
            if ( (k % gamma_inter) == 0 ) {
                // printf("%ld %ld \n", idx, gamma_idx );
                gamma_list[gamma_idx] = calcGamma(sites, size);
                gamma_step[gamma_idx] = k;  
                gamma_idx += 1;
            }
        }

        // Positions and gammas after the iteration 
        file_positions_cur = file_positions + "_" + to_string(round( rho_b * 10000000.0 ) / 10000000.0) + ".txt"  ;
        writePositions(sites, size, file_positions_cur);

        file_gamma_cur = file_gamma + "_" + to_string(round( rho_b * 10000000.0 ) / 10000000.0) + ".txt"  ;
        gamma_list[gamma_idx] = calcGamma(sites, size);
        gamma_step[gamma_idx] = steps+1;          
        writeGamma(gamma_list, gamma_step, file_gamma_cur );
     
        long count_particles = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if ( sites[i][j] == 1 ) { count_particles += 1; }
            }
        }             

        printf("For rho_bulk = %f, Total particles number = %lu \n", rho_b, count_particles);
    }    
    return 0;
}

