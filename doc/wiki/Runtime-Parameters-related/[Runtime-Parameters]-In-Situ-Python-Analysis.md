Parameters described on this page:
[YT_SCRIPT](#YT_SCRIPT), &nbsp;
[YT_VERBOSE](#YT_VERBOSE), &nbsp;
[YT_FIG_BASENAME](#YT_FIG_BASENAME) &nbsp;
[YT_JUPYTER_USE_CONNECTION_FILE](#YT_JUPYTER_USE_CONNECTION_FILE) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="YT_SCRIPT"></a>
* #### `YT_SCRIPT` &ensp; (string) &ensp; [none]
    * **Description:**
The Python script name to be loaded.
Do not include the file extension `.py`.
For example, set it to `inline_script` if the script name is `inline_script.py`.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--libyt | [Installation]-Option-List#--libyt]].

<a name="YT_VERBOSE"></a>
* #### `YT_VERBOSE` &ensp; (0=off, 1=info, 2=warning, 3=debug) &ensp; [1]
    * **Description:**
Log level of libyt.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--libyt | [Installation]-Option-List#--libyt]].

<a name="YT_FIG_BASENAME"></a>
* #### `YT_FIG_BASENAME` &ensp; (string) &ensp; [none]
    * **Description:**
Figure basename of the outputs from yt.
For example, set it to `Fig` will get `Fig000000000_Projection_z_density.png`
for the first figure of density z-projection.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--libyt | [Installation]-Option-List#--libyt]].

<a name="YT_JUPYTER_USE_CONNECTION_FILE"></a>
* #### `YT_JUPYTER_USE_CONNECTION_FILE` &ensp; (0=off, 1=on) &ensp; [0]
    * **Description:**
Use user-provided connection file when using libyt Jupyter UI.
    * **Restriction:**
Only applicable when enabling the compilation option
[[--libyt | [Installation]-Option-List#--libyt]] and
[[--libyt_jupyter | [Installation]-Option-List#--libyt_jupyter]].


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of In Situ Python Analysis | In-Situ-Python-Analysis]]
