//
// Created by Bo on 4/5/2016.
//

#ifndef REMUSIMM_CONTROLLER_H
#define REMUSIMM_CONTROLLER_H

#include <iostream>
#include <cmath>
#include "remus.h"

using namespace std;

// control the position of the internal moving mass to achieve heading stabilization
class MovingMassController
{
public:
    // constructor
    MovingMassController(const Remus&, const double&, const double&, const double&, const double&, const double&, const double&);
    MovingMassController(const Remus& vehicle): MovingMassController(vehicle, 0.5, 1.0, 1.0, 1.0, 1.0, 0.0) {}

private:
    // controller parameters
    const double k0_, k1_, k2_, k3_, k4_;
    // reference heading direction
    const double heading_reference_;
    // intermediate parameters
    double A_, B_, C_, E_, F_, L_, M_, N_, Nq_, epsilon1_, epsilon2_, epsilon3_;
    double beta1_, beta2_, beta3_, beta4_, gamma1_, gamma2_, gamma3_, gamma4_;
    double G_, dG_, H_, dH_, Gamma1_, Gamma2_, I_, dI_, J_, dJ_, Gamma4_;
    double R1_, R2_, G2_, H2_, I2_, J2_, N2_, N2q_, Gamma12_, Gamma22_, Gamma42_, epsilon12_, epsilon32_;
    double U0_, mv_;

};
// constructor
MovingMassController::MovingMassController(const Remus &vehicle, const double &k0, const double &k1, const double &k2,
                                           const double &k3, const double &k4, const double &heading_ref):
        k0_(k0), k1_(k1), k2_(k2), k3_(k3), k4_(k4), heading_reference_(heading_ref) {
    A_ = (vehicle.mass_+vehicle.a22_)/(vehicle.mass_+vehicle.a11_);
    B_ = -(vehicle.mass_+vehicle.a11_)/(vehicle.mass_+vehicle.a22_);
    C_ = -(vehicle.a22_-vehicle.a11_)/(vehicle.Izz_+vehicle.a66_);
    E_ = (vehicle.mass_*vehicle.cog_[0]+vehicle.a26_)/(vehicle.mass_+vehicle.a22_);
    F_ = (vehicle.mass_*vehicle.cog_[0]+vehicle.a26_)/(vehicle.Izz_+vehicle.a66_);
    L_ = -vehicle.mass_*vehicle.cog_[2]/(vehicle.mass_+vehicle.a22_);
    M_ = -vehicle.mass_*vehicle.cog_[2]/(vehicle.Ixx_+vehicle.a44_);
    N_ = vehicle.Kp_/(vehicle.Ixx_+vehicle.a44_);
    Nq_ = vehicle.Kpp_/(vehicle.Ixx_+vehicle.a44_);
    epsilon1_ = (vehicle.mass_-vehicle.displaced_mass_)*vehicle.gravity_/(vehicle.mass_+vehicle.a22_);
    epsilon2_ = vehicle.mass_*vehicle.gravity_*vehicle.cog_[0]/(vehicle.Izz_+vehicle.a66_);
    epsilon3_ = -vehicle.mass_*vehicle.gravity_*vehicle.cog_[2]/(vehicle.Ixx_+vehicle.a44_);
    beta1_ = vehicle.Yvv_/(vehicle.mass_+vehicle.a22_);
    beta2_ = vehicle.Yrr_/(vehicle.mass_+vehicle.a22_);
    beta3_ = vehicle.Yuv_/(vehicle.mass_+vehicle.a22_);
    beta4_ = vehicle.Yur_/(vehicle.mass_+vehicle.a22_);
    gamma1_ = vehicle.Nvv_/(vehicle.Izz_+vehicle.a66_);
    gamma2_ = vehicle.Nrr_/(vehicle.Izz_+vehicle.a66_);
    gamma3_ = vehicle.Nuv_/(vehicle.Izz_+vehicle.a66_);
    gamma4_ = vehicle.Nur_/(vehicle.Izz_+vehicle.a66_);
    G_ = -C_*E_/(1-E_*F_); dG_ = (beta3_-E_*gamma3_)/(1-E_*F_);
    H_ = (B_+E_*F_)/(1-E_*F_); dH_ = (beta4_-E_*gamma4_)/(1-E_*F_);
    Gamma1_ = (beta1_-E_*gamma1_)/(1-E_*F_); Gamma2_ = (beta2_-E_*gamma2_)/(1-E_*F_);
    I_ = C_/(1-E_*F_); dI_ = (gamma3_-F_*beta3_)/(1-E_*F_);
    J_ = -F_*(1+B_)/(1-E_*F_); dJ_ = (gamma4_-F_*beta4_)/(1-E_*F_);
    Gamma4_ = (gamma2_-F_*beta2_)/(1-E_*F_);
    R1_ = L_/(1-E_*F_); R2_ = (epsilon1_-epsilon2_*E_)/(1-E_*F_);
    G2_ = (G_+dG_)/(1-M_*R1_); H2_ = (H_+dH_+M_*R1_)/(1-M_*R1_);
    I2_ = (I_+dI_-M_*R1_*(C_+gamma3_))/(1-M_*R1_);
    J2_ = (J_+dJ_-M_*R1_*gamma4_)/(1-M_*R1_);
    N2_ = N_/(1-M_*R1_); N2q_ = Nq_/(1-M_*R1_);
    Gamma12_ = Gamma1_/(1-M_*R1_); Gamma22_ = Gamma2_/(1-M_*R1_);
    Gamma42_ = (Gamma4_-M_*R1_*gamma2_)/(1-M_*R1_);
    epsilon12_ = (R2_-epsilon3_*R1_)/(1-M_*R1_);
    epsilon32_ = (epsilon3_-M_*R2_)/(1-M_*R1_);
    U0_ = sqrt(-9/vehicle.Xuu_); mv_ = vehicle.mass_/6;
}

#endif //REMUSIMM_CONTROLLER_H
