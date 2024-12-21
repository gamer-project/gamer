Parameters described on this page:
[FB_LEVEL](#FB_LEVEL), &nbsp;
[FB_RSEED](#FB_RSEED), &nbsp;
[FB_SNE](#FB_SNE), &nbsp;
[FB_USER](#FB_USER) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**


<a name="FB_LEVEL"></a>
* #### `FB_LEVEL` &ensp; (0 &#8804; input < [[--nlevel | Installation:-Option-List#--nlevel]]; <0 &#8594; set to [[MAX_LEVEL | Runtime Parameters:-Refinement#MAX_LEVEL ]]) &ensp; [-1]
    * **Description:**
AMR level to apply feedback.
    * **Restriction:**
Must be [[MAX_LEVEL | Runtime Parameters:-Refinement#MAX_LEVEL ]] for now.

<a name="FB_RSEED"></a>
* #### `FB_RSEED` &ensp; (&#8805;0) &ensp; [456]
    * **Description:**
Random seed used by feedback.
    * **Restriction:**

<a name="FB_SNE"></a>
* #### `FB_SNE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Supernova explosion feedback.
    * **Restriction:**
Not supported yet.

<a name="FB_USER"></a>
* #### `FB_USER` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
User-defined feedback.
See [Add User-defined Feedback](#add-user-defined-feedback) for details.
    * **Restriction:**


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of Feedback | Feedback]]
