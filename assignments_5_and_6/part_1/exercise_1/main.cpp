#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>


// Definition der Lane-Emden-Gleichung
void lane_emden_eq(double xi, double w, double z, double n, double& dw_dxi, double& dz_dxi){
    dw_dxi = (xi==0.0) ? 0.0 : z;
    dz_dxi = (xi==0.0) ? 0.0 : - (2.0/xi) * z - pow(w, n);
}

// ein Schritt mit Heun-Verfahren maachen
void heun_step(double xi, double& w, double& z, double delta_xi, double n){

    double dw1, dw2, dz1, dz2; // Steigungen vor (1) und nach (2) Schritt

    lane_emden_eq(xi, w, z, n, dw1, dz1);

    // w und z nach einem Schritt
    double w_new = w + delta_xi*dw1;
    double z_new = z + delta_xi*dz1;

    lane_emden_eq(xi+delta_xi, w_new, z_new, n, dw2, dz2); // neue Ableitungen dw2, dz2 nach einem Schritt delta_xi

    // neue w und z
    w += 0.5 * delta_xi*(dw1 + dw2);
    z += 0.5 * delta_xi*(dz1 + dz2);

}


int main(){
    std::vector<double> n_values = {0.0, 1.0, 5.0};

    double delta_xi = 0.001;

    for (double n : n_values){
        // Anfangsbedingungen
        double xi = 0;
        double w = 1.0;
        double z = 0.0;

        // Daten in Datei schreiben
        std::ofstream file("lane_emden_n" + std::to_string((int)n) + ".dat");
        file << "# xi;w\n"; // schreibt Kommentarzeile mit Spaltenüberschrift

        while (w>=0.0 && xi < 15){
            file << xi << ";" << w << "\n";
            heun_step(xi, w, z, delta_xi, n);
            xi += delta_xi;
        }

        file.close();
        std::cout << "Für n = " << n << " ist der Radius r = " << xi << std::endl;

    }

    return 0;
}