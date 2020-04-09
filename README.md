# Scatter parameter extraction and optimization
`scatter_opt` is a MatLab package for processing electromagnetic scattering parameter data.

## Features
`scatter_opt` currently implements the following with a user-friendly programmatic interface:
* Nicholson-Ross-Weir (NRW) [1,2] and new non-iterative (NNI) [3] methods for permittivity and permeability extraction from S parameters
* Group delay method for branch selection as described by Weir [2]
* Kramers-Kronig method for branch selection as described by Szabó et al. [4]
* A new method for branch selection that is valid at high frequencies (unlike group delay) and does not suffer from truncation error (unlike the conventional Kramers-Kronig method)
* Nonlinear least-squares optimization of permittivity and permeability as described by  Domich, Baker-Jarvis, and Geyer [5,6]

<figure>
  <img src="https://github.com/jdhuang-csm/scatter-opt/blob/master/images/branch_rationalfit.jpg" width="700">
  <figcaption><i>Automatic branch determination using new rational fit method</i></figcaption>
</figure>

<figure>
  <img src="https://github.com/jdhuang-csm/scatter-opt/blob/master/images/tef_nrw.jpg" width="700">
  <figcaption><i>Permittivity and permeability extracted via the NRW method</i></figcaption>
</figure>
  
<figure>
  <img src="https://github.com/jdhuang-csm/scatter-opt/blob/master/images/tef_PPfit.jpg" width="700">
  <figcaption><i>Optimized permittivity and permeability</i></figcaption>
</figure>

<figure>
  <img src="https://github.com/jdhuang-csm/scatter-opt/blob/master/images/tef_Sfit.jpg" width="700">
  <figcaption><i>Optimized fit of measured S parameters</i></figcaption>
</figure>

References
1. Nicolson, A. M., & Ross, G. F. (1970). Measurement of the Intrinsic Properties Of Materials by Time-Domain Techniques. *IEEE Transactions on Instrumentation and Measurement, 19*(4), 377–382. https://doi.org/10.1109/TIM.1970.4313932
1. Weir, W. B. (1974). Automatic Measurement of Complex Dielectric Constant and Permeability at Microwave Frequencies. *Proceedings of the IEEE, 62*(1), 33–36. https://doi.org/10.1109/PROC.1974.9382
1. Boughriet, A. H., Legrand, C., & Chapoton, A. (1997). Noniterative stable transmission/reflection method for low-loss material complex permittivity determination. *IEEE Transactions on Microwave Theory and Techniques, 45*(1), 52–57. https://doi.org/10.1109/22.552032
1. Szabó, Z., Park, G. H., Hedge, R., & Li, E. P. (2010). A unique extraction of metamaterial parameters based on Kramers-Kronig relationship. *IEEE Transactions on Microwave Theory and Techniques, 58*(10), 2646–2653. https://doi.org/10.1109/TMTT.2010.2065310
1. Domich, P. D., Baker-Jarvis, J., & Geyer, R. G. (1991). Optimization techniques for permittivity and permeability determination. *Journal of Research of the National Institute of Standards and Technology, 96*(5), 565–575. https://doi.org/10.6028/jres.096.033
1. Baker-Jarvis, J., Geyer, R. G., & Domich, P. D. (1992). A Nonlinear Least-Squares Solution with Causality Constraints Applied to Transmission Line Permittivity and Permeability Determination. *IEEE Transactions on Instrumentation and Measurement, 41*(5), 646–652. https://doi.org/10.1109/19.177336
