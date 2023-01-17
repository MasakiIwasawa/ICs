#include"mkpl.hpp"

void check_status(const PS::F64 * mass,
		  const PS::F64vec * pos,
		  const PS::F64vec * vel,
		  const int n){
  PS::F64 eng_kin = 0.0;
  PS::F64 eng_pot = 0.0;
  PS::F64vec pos_cm = 0.0;
  PS::F64vec vel_cm = 0.0;
  PS::F64 mass_tot = 0.0;
  for(int i=0; i<n; i++){
    eng_kin += 0.5*mass[i]*vel[i]*vel[i];
    pos_cm += mass[i]*pos[i];
    vel_cm += mass[i]*vel[i];
    mass_tot += mass[i];
    for(int j=i+1; j<n; j++){
      PS::F64vec rij = pos[i] - pos[j];
      PS::F64 over_r = 1.0 / sqrt(rij*rij);
      eng_pot -= mass[i]*mass[j] * over_r;
    }
  }

  std::cerr<<"eng_kin= "<<eng_kin<<" eng_pot= "<<eng_pot<<std::endl;
  pos_cm /= mass_tot;
  vel_cm /= mass_tot;
  std::cerr<<"pos_cm= "<<pos_cm<<" vel_cm= "<<vel_cm<<std::endl;
}

int main(){
  PS::F64 * mass = nullptr;
  PS::F64vec * pos = nullptr;
  PS::F64vec * vel = nullptr;
  PS::F64 mass_glb = 1.0;
  PS::S64 n_glb = 1024;
  PS::S64 n_loc = 1024;
  MakePlummerModel(mass_glb, n_glb, n_loc, mass, pos, vel);
  //for(int i=0; i<n_loc; i++){
  //std::cerr<<"i="<<i<<" pos[i]= "<<pos[i]<<std::endl;
  //}
  check_status(mass, pos, vel, n_glb);
  return 0;
}
