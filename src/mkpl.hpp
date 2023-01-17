#include<iostream>
#include<random>
#include<cassert>
#include<ps_types.hpp>
//#include<particle_simulator.hpp>
namespace PS = ParticleSimulator;
void MakePlummerModel(const PS::F64 mass_glb,
                      const PS::S64 n_glb,
                      const PS::S64 n_loc,
                      PS::F64 *& mass,
                      PS::F64vec *& pos,
                      PS::F64vec *& vel,
                      const PS::F64 eng = -0.25,
                      const PS::S32 seed = 0){

    assert(eng < 0.0);
    static const PS::F64 PI = atan(1.0) * 4.0;
    const PS::F64 r_cutoff = 22.8 / (-3.0 * PI * mass_glb * mass_glb / (64.0 * -0.25)); // 22.8 is cutoff in Nbody units
    //const PS::F64 r_cutoff = 22.8 * 0.25;
    mass = new PS::F64[n_loc];
    pos  = new PS::F64vec[n_loc];
    vel  = new PS::F64vec[n_loc];
    //PS::MTTS mt;
    //mt.init_genrand( PS::Comm::getRank() );
    std::mt19937_64 mt(seed);
    std::uniform_real_distribution<PS::F64> dist(0.0, 1.0);
    for(int i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
        PS::F64 r_tmp = 9999.9;
        while(r_tmp > r_cutoff){
	    //PS::F64 m_tmp = mt.genrand_res53();
	    PS::F64 m_tmp = dist(mt);
            r_tmp = 1.0 / sqrt( pow(m_tmp, (-2.0/3.0)) - 1.0);
        }
        PS::F64 phi = 2.0 * PI * dist(mt);
        PS::F64 cth = 2.0 * (dist(mt) - 0.5);
        PS::F64 sth = sqrt(1.0 - cth*cth);
        pos[i][0] = r_tmp * sth * cos(phi);
        pos[i][1] = r_tmp * sth * sin(phi);
        pos[i][2] = r_tmp * cth;
        while(1){
            const PS::F64 v_max = 0.1;
            const PS::F64 v_try = dist(mt);
            const PS::F64 v_crit = v_max * dist(mt);
            if(v_crit < v_try * v_try * pow( (1.0 - v_try * v_try), 3.5) ){
                const PS::F64 ve = sqrt(2.0) * pow( (r_tmp*r_tmp + 1.0), -0.25 );
                phi = 2.0 * PI * dist(mt);
                cth = 2.0 * (dist(mt) - 0.5);
                sth = sqrt(1.0 - cth*cth);
                vel[i][0] = ve * v_try * sth * cos(phi);
                vel[i][1] = ve * v_try * sth * sin(phi);
                vel[i][2] = ve * v_try * cth;
                break;
            }
        }
    }
    PS::F64vec cm_pos = 0.0;
    PS::F64vec cm_vel = 0.0;
    PS::F64  cm_mass = 0.0;
    for(auto i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(auto i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }
    const PS::F64 r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
    const PS::F64 coef = 1.0 / sqrt(r_scale);
    for(auto i=0; i<n_loc; i++){
        pos[i] *= r_scale;
        vel[i] *= coef;
    }
    PS::F64 r_max_sq = -1.0;
    for(auto i=0; i<n_loc; i++){
        if(r_max_sq < pos[i] * pos[i]){
            r_max_sq = pos[i] * pos[i];
        }
    }
}


/*
template<class Tpsys>
void SetParticlesPlummer(Tpsys & psys,
                         const PS::S64 n_glb,
                         PS::S32 & n_loc,  
                         PS::F32 & t_sys){

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    t_sys = 0.0;

    const PS::F64 m_tot = 1.0;
    const PS::F64 eng = -0.25;
    MakePlummerModel(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(size_t i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}
*/
