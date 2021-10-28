# HPW_MembraneDiffusion


## Circular inclusions: Hughes, Pailthorpe, and White model

MATLAB Code for implementing the series solutions to the membrane hydrodynamics equations of Hughes, Pailthorpe, and White (1981):

Hughes, B.D., B.A. Pailthorpe, and White L. R., "The translational and rotational drag on a cylinder moving in a membrane," *J. Fluid. Mech.* **110:** 349–372  (1981).

The function **HPW_drag.m** calculates the translational and rotational drag and diffusion coefficients of a membrane inclusion, given input viscosities and the inclusion radius.  

**HPW_fit.m** calculates the best-fit membrane viscosity and inclusion radius, given (measured) rotational and translational diffusion coefficients.

See https://eighteenthelephant.com/2018/09/04/membrane-diffusion-software/ for discussion. 

## Rod-like inclusions: Levine

Functions applying the theory of: 

\A. J. Levine, T. B. Liverpool, F. C. MacKintosh, 
"Mobility of extended bodies in viscous films and membranes,"  *Phys. Rev. E* **69,** 021503 (2004).

**LLM_calcD.m** calculates the rotational diffusion coefficient and the translational diffusion coefficients parallel and perpendicular to the rod axis as a function of membrane viscosity (eta) and rod length (L).

**LLM_fit_viscosity_L.m** calulates membrane viscosity and effective rod length for a rod-like membrane inclusion. Uses both D_parallel and D_R, and D_perpendicular and D_R, combining the results from each pair. (October 2021) Note: we've extracted points from the graphs in the Levine et al. 2004 paper, and made simple polynomial interpolations. See notes, 2021, and Jahl and Parthasarathy 2021.



