---
title: "Desiccation model"
author: "Remko"
date: "23 June 2017"
output:
  pdf_document: default
  html_document: default
---



| Parameter |  Definition         | Units     |
|-----------|-------------------- | --------- |
| J~rs~     | Root to stem flow   | mol s^-1^   |
| J~sl~     | Stem to leaf flow   | mol s^-1^   |
| E~t~      | Transpiration       | mol s^-1^   |
| c~s~      | Specific stem capacitance | RWC MPa^-1^ |
| C~s~      | Absolute stem capacitance | mol MPa^-1^ |
| c~l~      | Specific leaf capacitance | RWC MPa^-1^ |
| C~l~      | Absolute leaf capacitance | mol MPa^-1^ |
| k~S~      | soil to root hydr. conductance  | mol m^-2^ s^-1^ MPa^-1^ |
| k~P~      | plant hydr. conductance | mol m^-2^ s^-1^ MPa^-1^ |
| A~L~      | Plant leaf area | m^2^ |
| g~min~    | Minimum conductance | mol m-2 s-1
|  $\Psi_{st}$ | stem water potential | MPa |
|  $\Psi_{S}$ | soil water potential | MPa |
|  $\Psi_{L}$ | leaf water potential | MPa |
| $\Psi_{min}$ | Minimum leaf water potential (turgor loss point) | MPa |
| $\rho_w$   | wood density | kg m^-3^ |
| $\rho_s$   | non-lumen wood density (constant) | kg m^-3^ |
| V~s~       | sapwood volume | m^3^ |


$$J_{rs} = A_L k_S (\Psi_S - \Psi_{st})$$

$$J_{sl} = A_L k_P (\Psi_{st} - \Psi_L)$$

$$E_t = A_L f(D,PPFD, \Psi_{min}, g_{min}, etc)$$

$$C_s = c_s V_s (1 - \rho_w/\rho_s)$$ 

$$C_l = c_l A_L f(LMA?)$$ 

$$\frac{d\Psi_{st}}{dt} = \frac{J_{rs} - J_{sl}}{C_s}$$

$$\frac{d\Psi_{l}}{dt} = \frac{J_{sl} - E_t}{C_s}$$




