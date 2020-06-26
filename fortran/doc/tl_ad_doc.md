#  Modular arbitrary-order ocean-atmosphere model: The Tangent Linear and Adjoint model #

## Description : ##

 The Tangent Linear and Adjoint model model are implemented in the same way as the nonlinear model, with a tensor storing the different terms. The Tangent Linear (TL) tensor \f$\mathcal{T}_{i,j,k}^{TD}\f$ is defined as:

\f[ \mathcal{T}_{i,j,k}^{TL} = \mathcal{T}_{i,k,j} + \mathcal{T}_{i,j,k} \f]

while the Adjoint (AD) tensor \f$\mathcal{T}_{i,j,k}^{AD}\f$ is defined as:

\f[ \mathcal{T}_{i,j,k}^{AD} = \mathcal{T}_{j,k,i} + \mathcal{T}_{j,i,k} . \f]

where \f$ \mathcal{T}_{i,j,k}\f$ is the tensor of the nonlinear model.

These two tensors are used to compute the trajectories of the models, with the equations

\f[  \frac{d\delta y_i}{dt} = \sum_{j=1}^{ndim}\sum_{k=0}^{ndim} \, \mathcal{T}_{i,j,k}^{TL} \, y^{\ast}_k \; \delta y_j . \f]

\f[   -\frac{d\delta y_i}{dt} = \sum_{j=1}^{ndim} \sum_{k=0}^{ndim} \, \mathcal{T}_{i,j,k}^{AD} \, y^{\ast}_k \; \delta y_j . \f]

where \f$\boldsymbol{y}^{\ast}\f$ is the point where the Tangent model is defined (with \f$y_0^{\ast}=1\f$).

## Implementation : ##

The two tensors are implemented in the module tl_ad_tensor and must be initialized inside a given model_def::model object with the method model_def::init_tl_model and model_def::init_ad_model. The tendencies are then given by the routine model_def::tl_tendencies and model_def::ad_tendencies. Integrators with the Heun method (RK2) or the 4th-order Runge-Kutta method are available with the classes rk2_tl_integrator, rk2_ad_integrator, rk4_tl_integrator and rk4_ad_integrator. An example on how to use it can be found in the test file test_tl_ad.f90
