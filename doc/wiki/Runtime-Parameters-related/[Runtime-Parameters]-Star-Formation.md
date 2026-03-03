Parameters described on this page:
[SF_CREATE_STAR_SCHEME](#SF_CREATE_STAR_SCHEME), &nbsp;
[SF_CREATE_STAR_RSEED](#SF_CREATE_STAR_RSEED), &nbsp;
[SF_CREATE_STAR_DET_RANDOM](#SF_CREATE_STAR_DET_RANDOM), &nbsp;
[SF_CREATE_STAR_MIN_LEVEL](#SF_CREATE_STAR_MIN_LEVEL), &nbsp;
[SF_CREATE_STAR_MIN_GAS_DENS](#SF_CREATE_STAR_MIN_GAS_DENS), &nbsp;
[SF_CREATE_STAR_MASS_EFF](#SF_CREATE_STAR_MASS_EFF), &nbsp;
[SF_CREATE_STAR_MIN_STAR_MASS](#SF_CREATE_STAR_MIN_STAR_MASS), &nbsp;
[SF_CREATE_STAR_MAX_STAR_MFRAC](#SF_CREATE_STAR_MAX_STAR_MFRAC), &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="SF_CREATE_STAR_SCHEME"></a>
* #### `SF_CREATE_STAR_SCHEME` &ensp; (0=off, 1=AGORA) &ensp; [0]
    * **Description:**
Star formation schemes.
See Sec.3.2 in [Kim et al. 2016](https://iopscience.iop.org/article/10.3847/1538-4357/833/2/202)
and Sec.2.4 in [Goldbaum et al. 2015](https://iopscience.iop.org/article/10.1088/0004-637X/814/2/131)
for the AGORA star formation scheme.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | [Installation]-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_RSEED"></a>
* #### `SF_CREATE_STAR_RSEED` &ensp; (&#8805;0) &ensp; [123]
    * **Description:**
Random seed used by star formation.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | [Installation]-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_DET_RANDOM"></a>
* #### `SF_CREATE_STAR_DET_RANDOM` &ensp; (0=off, 1=on; <0 &#8594; set by [[--bitwise_reproducibility | [Installation]-Option-List#--bitwise_reproducibility]]) &ensp; [-1]
    * **Description:**
Make random numbers deterministic (i.e., independent of OpenMP and MPI).
Needed for the [[bitwise reproducibility | Bitwise Reproducibility]].
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | [Installation]-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_MIN_LEVEL"></a>
* #### `SF_CREATE_STAR_MIN_LEVEL` &ensp; (&#8805;0; <0 &#8594; [[MAX_LEVEL | [Runtime Parameters]-Refinement#MAX_LEVEL]]) &ensp; [0]
    * **Description:**
Minimum AMR level allowed to form stars.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | [Installation]-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_MIN_GAS_DENS"></a>
* #### `SF_CREATE_STAR_MIN_GAS_DENS` &ensp; (&#8805;0.0) &ensp; [1.0e1]
    * **Description:**
Minimum gas density allowed to form stars.
See $\rho_{\rm gas,\ thres}$ in Eq.(4) in [Kim et al. 2016](https://iopscience.iop.org/article/10.3847/1538-4357/833/2/202).
Note that the input value should always be in units of HI count/cm^3,
and it will be converted internally to the gas mass density
as $m_H\times$ HI count/cm^3 (i.e., assuming the gas is composed of only HI and the mean molecular weight $\mu=1$).
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | [Installation]-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_MASS_EFF"></a>
* #### `SF_CREATE_STAR_MASS_EFF` &ensp; (0.0 < input &#8804; 1.0) &ensp; [1.0e-2]
    * **Description:**
Gas-to-star mass conversion efficiency.
See $\epsilon_*$ in Eq.(4) in [Kim et al. 2016](https://iopscience.iop.org/article/10.3847/1538-4357/833/2/202).
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | [Installation]-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_MIN_STAR_MASS"></a>
* #### `SF_CREATE_STAR_MIN_STAR_MASS` &ensp; (&#8805;0.0) &ensp; [0.0]
    * **Description:**
Minimum star particle mass for the stochastical star formation.
See $m_{\rm sf}$ in Eq.(5) in [Goldbaum et al. 2015](https://iopscience.iop.org/article/10.1088/0004-637X/814/2/131).
Note that the input value should always be in units of Msun.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | [Installation]-Option-List#--star_formation]].

<a name="SF_CREATE_STAR_MAX_STAR_MFRAC"></a>
* #### `SF_CREATE_STAR_MAX_STAR_MFRAC` &ensp; (0.0 < input &#8804; 1.0) &ensp; [0.5]
    * **Description:**
Maximum gas mass fraction in a cell allowed to be converted to stellar mass per substep.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--star_formation | [Installation]-Option-List#--star_formation]].



## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Star Formation | Star-Formation]]
