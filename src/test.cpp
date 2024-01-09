#include <iostream>
#include <string>
#include <fstream>
#include "board.cpp"

void run(Board_2IdealGases* b, int n, std::string results_folser, int when_to_make_shots) {
    double tau;
    double time = 0;
    int shot = 0;

    std::ofstream timestems("./" + results_folser + "/timesteps.txt", std::ios::app);

    for (int i = 0; i < n; i++) {
        tau = b->get_tau_max() / 3;
        time += tau;        
        std::cout << "### " << i << "\n";
        std::cout << "  # time step :" << tau << "\n";
        
        b->make_one_step(tau);

        if (i % when_to_make_shots == 0) {
            timestems << time << '\n';            
            b->write_everything(results_folser, shot);
            shot++;
        }
    }
    timestems.close();
}

int main() {
    Board_2IdealGases b(200, 0.1);

    b.experiment_02_initial_state();

    run(&b, 101, "results", 10);

    return 0;
}
