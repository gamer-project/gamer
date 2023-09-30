#include "Macro.h"
#include "CUPOT.h"

#if ( defined GRAVITY  &&  defined GPU  &&  POT_SCHEME == DST )


#if (RHO_NXT==24)
__device__ real M[24*24] = {0.12533323356430426 , 0.2486898871648548 , 0.3681245526846779 , 0.4817536741017153 , 0.5877852522924731 , 0.6845471059286886 , 0.7705132427757891 , 0.8443279255020151 , 0.9048270524660196 , 0.9510565162951535 , 0.9822872507286886 , 0.9980267284282716 , 0.9980267284282716 , 0.9822872507286887 , 0.9510565162951536 , 0.9048270524660195 , 0.844327925502015 , 0.7705132427757893 , 0.6845471059286888 , 0.5877852522924732 , 0.4817536741017152 , 0.36812455268467814 , 0.24868988716485524 , 0.12533323356430454 , 0.2486898871648548 , 0.4817536741017153 , 0.6845471059286886 , 0.8443279255020151 , 0.9510565162951535 , 0.9980267284282716 , 0.9822872507286887 , 0.9048270524660195 , 0.7705132427757893 , 0.5877852522924732 , 0.36812455268467814 , 0.12533323356430454 , -0.1253332335643043 , -0.3681245526846779 , -0.5877852522924727 , -0.7705132427757894 , -0.9048270524660198 , -0.9822872507286887 , -0.9980267284282716 , -0.9510565162951536 , -0.8443279255020151 , -0.684547105928689 , -0.4817536741017161 , -0.24868988716485535 , 0.3681245526846779 , 0.6845471059286886 , 0.9048270524660196 , 0.9980267284282716 , 0.9510565162951536 , 0.7705132427757893 , 0.4817536741017152 , 0.12533323356430454 , -0.24868988716485457 , -0.5877852522924727 , -0.8443279255020147 , -0.9822872507286887 , -0.9822872507286886 , -0.8443279255020151 , -0.587785252292474 , -0.24868988716485535 , 0.12533323356430506 , 0.48175367410171493 , 0.7705132427757887 , 0.9510565162951532 , 0.9980267284282716 , 0.90482705246602 , 0.6845471059286884 , 0.368124552684678 , 0.4817536741017153 , 0.8443279255020151 , 0.9980267284282716 , 0.9048270524660195 , 0.5877852522924732 , 0.12533323356430454 , -0.3681245526846779 , -0.7705132427757894 , -0.9822872507286887 , -0.9510565162951536 , -0.684547105928689 , -0.24868988716485535 , 0.2486898871648549 , 0.6845471059286886 , 0.9510565162951532 , 0.9822872507286886 , 0.7705132427757886 , 0.368124552684678 , -0.12533323356430318 , -0.5877852522924728 , -0.9048270524660197 , -0.9980267284282716 , -0.8443279255020161 , -0.48175367410171627 , 0.5877852522924731 , 0.9510565162951535 , 0.9510565162951536 , 0.5877852522924732 , 1.2246467991473532e-16 , -0.5877852522924727 , -0.9510565162951535 , -0.9510565162951536 , -0.5877852522924732 , -2.4492935982947064e-16 , 0.5877852522924722 , 0.9510565162951532 , 0.9510565162951536 , 0.5877852522924734 , 3.6739403974420594e-16 , -0.5877852522924728 , -0.9510565162951534 , -0.9510565162951538 , -0.5877852522924735 , -4.898587196589413e-16 , 0.5877852522924728 , 0.9510565162951529 , 0.9510565162951543 , 0.5877852522924751 , 0.6845471059286886 , 0.9980267284282716 , 0.7705132427757893 , 0.12533323356430454 , -0.5877852522924727 , -0.9822872507286887 , -0.8443279255020151 , -0.24868988716485535 , 0.48175367410171493 , 0.9510565162951532 , 0.90482705246602 , 0.368124552684678 , -0.3681245526846789 , -0.9048270524660197 , -0.9510565162951543 , -0.48175367410171627 , 0.24868988716485635 , 0.8443279255020146 , 0.9822872507286889 , 0.5877852522924751 , -0.1253332335643047 , -0.7705132427757879 , -0.9980267284282716 , -0.6845471059286887 , 0.7705132427757891 , 0.9822872507286887 , 0.4817536741017152 , -0.3681245526846779 , -0.9510565162951535 , -0.8443279255020151 , -0.12533323356430467 , 0.6845471059286886 , 0.9980267284282716 , 0.5877852522924734 , -0.24868988716485302 , -0.9048270524660197 , -0.9048270524660185 , -0.2486898871648556 , 0.5877852522924714 , 0.9980267284282716 , 0.6845471059286886 , -0.1253332335643047 , -0.8443279255020135 , -0.9510565162951538 , -0.36812455268467664 , 0.4817536741017121 , 0.9822872507286884 , 0.7705132427757888 , 0.8443279255020151 , 0.9048270524660195 , 0.12533323356430454 , -0.7705132427757894 , -0.9510565162951536 , -0.24868988716485535 , 0.6845471059286886 , 0.9822872507286886 , 0.368124552684678 , -0.5877852522924728 , -0.9980267284282716 , -0.48175367410171627 , 0.4817536741017155 , 0.9980267284282716 , 0.5877852522924751 , -0.3681245526846787 , -0.982287250728689 , -0.6845471059286887 , 0.2486898871648527 , 0.9510565162951533 , 0.7705132427757888 , -0.12533323356430268 , -0.9048270524660179 , -0.8443279255020163 , 0.9048270524660196 , 0.7705132427757893 , -0.24868988716485457 , -0.9822872507286887 , -0.5877852522924732 , 0.48175367410171493 , 0.9980267284282716 , 0.368124552684678 , -0.684547105928688 , -0.9510565162951538 , -0.12533323356430578 , 0.8443279255020146 , 0.8443279255020152 , -0.1253332335643047 , -0.9510565162951534 , -0.6845471059286887 , 0.3681245526846786 , 0.9980267284282714 , 0.48175367410171666 , -0.5877852522924724 , -0.9822872507286887 , -0.24868988716485782 , 0.7705132427757877 , 0.9048270524660202 , 0.9510565162951535 , 0.5877852522924732 , -0.5877852522924727 , -0.9510565162951536 , -2.4492935982947064e-16 , 0.9510565162951532 , 0.5877852522924734 , -0.5877852522924728 , -0.9510565162951538 , -4.898587196589413e-16 , 0.9510565162951529 , 0.5877852522924751 , -0.5877852522924726 , -0.9510565162951538 , -7.347880794884119e-16 , 0.9510565162951533 , 0.5877852522924738 , -0.5877852522924724 , -0.9510565162951539 , -9.797174393178826e-16 , 0.9510565162951532 , 0.5877852522924769 , -0.5877852522924694 , -0.951056516295155 , 0.9822872507286886 , 0.36812455268467814 , -0.8443279255020153 , -0.684547105928689 , 0.5877852522924729 , 0.9048270524660192 , -0.24868988716485474 , -0.9980267284282716 , -0.12533323356430404 , 0.9510565162951534 , 0.4817536741017165 , -0.7705132427757901 , -0.7705132427757888 , 0.4817536741017152 , 0.951056516295155 , -0.12533323356430268 , -0.9980267284282717 , -0.24868988716485438 , 0.9048270524660179 , 0.587785252292474 , -0.68454710592869 , -0.8443279255020164 , 0.3681245526846748 , 0.9822872507286882 , 0.9980267284282716 , 0.12533323356430454 , -0.9822872507286887 , -0.24868988716485535 , 0.9510565162951532 , 0.368124552684678 , -0.9048270524660197 , -0.48175367410171627 , 0.8443279255020146 , 0.5877852522924751 , -0.7705132427757879 , -0.6845471059286887 , 0.6845471059286902 , 0.7705132427757888 , -0.5877852522924695 , -0.8443279255020163 , 0.4817536741017181 , 0.9048270524660202 , -0.3681245526846749 , -0.951056516295155 , 0.24868988716485566 , 0.9822872507286895 , -0.1253332335643057 , -0.9980267284282716 , 0.9980267284282716 , -0.1253332335643043 , -0.9822872507286887 , 0.2486898871648549 , 0.9510565162951536 , -0.36812455268467725 , -0.90482705246602 , 0.4817536741017155 , 0.8443279255020152 , -0.5877852522924726 , -0.770513242775791 , 0.6845471059286876 , 0.6845471059286887 , -0.7705132427757878 , -0.5877852522924739 , 0.8443279255020153 , 0.48175367410171377 , -0.9048270524660194 , -0.3681245526846804 , 0.9510565162951532 , 0.24868988716485124 , -0.9822872507286876 , -0.125333233564305 , 0.9980267284282713 , 0.9822872507286887 , -0.3681245526846779 , -0.8443279255020151 , 0.6845471059286886 , 0.5877852522924734 , -0.9048270524660197 , -0.2486898871648556 , 0.9980267284282716 , -0.1253332335643047 , -0.9510565162951538 , 0.4817536741017121 , 0.7705132427757888 , -0.7705132427757924 , -0.48175367410171677 , 0.9510565162951522 , 0.12533323356430465 , -0.9980267284282716 , 0.24868988716485566 , 0.9048270524660219 , -0.587785252292472 , -0.6845471059286866 , 0.8443279255020112 , 0.36812455268468075 , -0.9822872507286889 , 0.9510565162951536 , -0.5877852522924727 , -0.5877852522924732 , 0.9510565162951532 , 3.6739403974420594e-16 , -0.9510565162951538 , 0.5877852522924728 , 0.5877852522924751 , -0.9510565162951534 , -7.347880794884119e-16 , 0.951056516295155 , -0.5877852522924724 , -0.587785252292471 , 0.9510565162951532 , 1.102182119232618e-15 , -0.951056516295155 , 0.587785252292475 , 0.5877852522924741 , -0.9510565162951521 , -1.4695761589768238e-15 , 0.951056516295153 , -0.5877852522924661 , -0.5877852522924774 , 0.951056516295153 , 0.9048270524660195 , -0.7705132427757894 , -0.24868988716485535 , 0.9822872507286886 , -0.5877852522924728 , -0.48175367410171627 , 0.9980267284282716 , -0.3681245526846787 , -0.6845471059286887 , 0.9510565162951533 , -0.12533323356430268 , -0.8443279255020163 , 0.8443279255020153 , 0.12533323356430465 , -0.951056516295155 , 0.6845471059286898 , 0.3681245526846739 , -0.9980267284282716 , 0.4817536741017115 , 0.5877852522924744 , -0.9822872507286889 , 0.24868988716485174 , 0.770513242775794 , -0.9048270524660176 , 0.844327925502015 , -0.9048270524660198 , 0.12533323356430418 , 0.7705132427757886 , -0.9510565162951534 , 0.24868988716485463 , 0.68454710592869 , -0.982287250728689 , 0.3681245526846786 , 0.5877852522924738 , -0.9980267284282714 , 0.481753674101715 , 0.48175367410171377 , -0.9980267284282718 , 0.5877852522924693 , 0.3681245526846739 , -0.9822872507286888 , 0.6845471059286896 , 0.2486898871648584 , -0.9510565162951541 , 0.770513242775794 , 0.12533323356430887 , -0.9048270524660207 , 0.8443279255020147 , 0.7705132427757893 , -0.9822872507286887 , 0.48175367410171493 , 0.368124552684678 , -0.9510565162951538 , 0.8443279255020146 , -0.1253332335643047 , -0.6845471059286887 , 0.9980267284282714 , -0.5877852522924724 , -0.24868988716485782 , 0.9048270524660202 , -0.9048270524660194 , 0.24868988716485566 , 0.5877852522924741 , -0.9980267284282716 , 0.6845471059286896 , 0.12533323356430864 , -0.8443279255020167 , 0.951056516295153 , -0.3681245526846776 , -0.48175367410172076 , 0.9822872507286896 , -0.770513242775787 , 0.6845471059286888 , -0.9980267284282716 , 0.7705132427757887 , -0.12533323356430318 , -0.5877852522924735 , 0.9822872507286889 , -0.8443279255020155 , 0.2486898871648527 , 0.48175367410171666 , -0.9510565162951539 , 0.9048270524660179 , -0.3681245526846749 , -0.3681245526846771 , 0.9048270524660189 , -0.9510565162951521 , 0.4817536741017115 , 0.24868988716485152 , -0.8443279255020167 , 0.9822872507286882 , -0.5877852522924716 , -0.12533323356430548 , 0.7705132427757941 , -0.9980267284282716 , 0.684547105928684 , 0.5877852522924732 , -0.9510565162951536 , 0.9510565162951532 , -0.5877852522924728 , -4.898587196589413e-16 , 0.5877852522924751 , -0.9510565162951538 , 0.9510565162951533 , -0.5877852522924724 , -9.797174393178826e-16 , 0.5877852522924769 , -0.951056516295155 , 0.9510565162951532 , -0.587785252292472 , -1.4695761589768238e-15 , 0.5877852522924744 , -0.9510565162951541 , 0.951056516295153 , -0.5877852522924716 , -1.959434878635765e-15 , 0.5877852522924748 , -0.9510565162951564 , 0.9510565162951506 , -0.5877852522924655 , 0.4817536741017152 , -0.8443279255020151 , 0.9980267284282716 , -0.9048270524660197 , 0.5877852522924728 , -0.1253332335643047 , -0.36812455268467664 , 0.7705132427757888 , -0.9822872507286887 , 0.9510565162951532 , -0.6845471059286847 , 0.24868988716485566 , 0.24868988716485124 , -0.6845471059286866 , 0.9510565162951552 , -0.9822872507286889 , 0.770513242775794 , -0.3681245526846776 , -0.12533323356430548 , 0.5877852522924748 , -0.9048270524660177 , 0.9980267284282709 , -0.8443279255020126 , 0.4817536741017169 , 0.36812455268467814 , -0.684547105928689 , 0.9048270524660192 , -0.9980267284282716 , 0.9510565162951534 , -0.7705132427757901 , 0.4817536741017152 , -0.12533323356430268 , -0.24868988716485438 , 0.587785252292474 , -0.8443279255020164 , 0.9822872507286882 , -0.9822872507286889 , 0.8443279255020151 , -0.5877852522924661 , 0.24868988716485174 , 0.12533323356430182 , -0.48175367410171455 , 0.7705132427757941 , -0.9510565162951543 , 0.9980267284282718 , -0.9048270524660174 , 0.6845471059286837 , -0.36812455268468347 , 0.24868988716485524 , -0.4817536741017161 , 0.6845471059286884 , -0.8443279255020161 , 0.9510565162951532 , -0.9980267284282716 , 0.9822872507286884 , -0.9048270524660179 , 0.7705132427757877 , -0.5877852522924751 , 0.3681245526846748 , -0.1253332335643057 , -0.12533323356429796 , 0.36812455268468075 , -0.5877852522924774 , 0.770513242775794 , -0.9048270524660207 , 0.9822872507286896 , -0.9980267284282716 , 0.9510565162951551 , -0.8443279255020126 , 0.6845471059286837 , -0.48175367410171366 , 0.24868988716485765 , 0.12533323356430454 , -0.24868988716485535 , 0.368124552684678 , -0.48175367410171627 , 0.5877852522924751 , -0.6845471059286887 , 0.7705132427757888 , -0.8443279255020163 , 0.9048270524660202 , -0.951056516295155 , 0.9822872507286895 , -0.9980267284282716 , 0.9980267284282718 , -0.9822872507286889 , 0.9510565162951509 , -0.9048270524660176 , 0.8443279255020186 , -0.770513242775787 , 0.684547105928684 , -0.5877852522924655 , 0.4817536741017169 , -0.36812455268467026 , 0.24868988716485765 , -0.12533323356430429};

#elif (RHO_NXT==16)
__device__ real M[16*16]= {0.18374951781657034 , 0.3612416661871529 , 0.5264321628773557 , 0.6736956436465572 , 0.7980172272802395 , 0.8951632913550623 , 0.961825643172819 , 0.9957341762950345 , 0.9957341762950346 , 0.961825643172819 , 0.8951632913550626 ,0.7980172272802396 , 0.6736956436465571 , 0.5264321628773561 , 0.3612416661871533 , 0.18374951781657037 , 0.3612416661871529 , 0.6736956436465572 , 0.8951632913550623 , 0.9957341762950345 , 0.961825643172819 , 0.7980172272802396 , 0.5264321628773561 , 0.18374951781657037 , -0.18374951781657015 , -0.5264321628773558 , -0.7980172272802388 , -0.961825643172819 , -0.9957341762950345 , -0.8951632913550626 , -0.6736956436465578 , -0.361241666187153 , 0.5264321628773557 , 0.8951632913550623 , 0.9957341762950346 , 0.7980172272802396 , 0.3612416661871533 , -0.18374951781657015 , -0.6736956436465572 , -0.961825643172819 , -0.961825643172819 , -0.6736956436465578 , -0.18374951781657092 , 0.3612416661871526 , 0.7980172272802399 , 0.9957341762950345, 0.8951632913550635 , 0.5264321628773563 , 0.6736956436465572 , 0.9957341762950345 , 0.7980172272802396 , 0.18374951781657037 , -0.5264321628773558 , -0.961825643172819 , -0.8951632913550626 , -0.361241666187153 , 0.3612416661871526 , 0.8951632913550623 , 0.9618256431728196 , 0.5264321628773563 , -0.1837495178165712 , -0.7980172272802387 , -0.9957341762950347 , -0.6736956436465573 , 0.7980172272802395 , 0.961825643172819 , 0.3612416661871533 , -0.5264321628773558 , -0.9957341762950346 , -0.6736956436465578 , 0.18374951781656956 , 0.8951632913550623 , 0.8951632913550627 , 0.18374951781657017 , -0.6736956436465568 , -0.9957341762950347 , -0.5264321628773548 , 0.3612416661871515 , 0.9618256431728189 , 0.7980172272802394 , 0.8951632913550623 , 0.7980172272802396 , -0.18374951781657015 , -0.961825643172819 , -0.6736956436465578 , 0.3612416661871526 , 0.9957341762950345 , 0.5264321628773563 , -0.5264321628773557 , -0.9957341762950347 , -0.3612416661871541 , 0.6736956436465567 , 0.9618256431728187 , 0.18374951781657042 , -0.7980172272802365 , -0.8951632913550628 , 0.961825643172819 , 0.5264321628773561 , -0.6736956436465572 , -0.8951632913550626 , 0.18374951781656956 , 0.9957341762950345 , 0.361241666187154 , -0.7980172272802387 ,-0.7980172272802393 , 0.3612416661871515 , 0.9957341762950347 , 0.18374951781657042 , -0.8951632913550638 , -0.6736956436465589 , 0.5264321628773523 , 0.9618256431728197 , 0.9957341762950345 , 0.18374951781657037 , -0.961825643172819 , -0.361241666187153 , 0.8951632913550623 , 0.5264321628773563 , -0.7980172272802387 , -0.6736956436465573 , 0.6736956436465567 , 0.7980172272802394 , -0.5264321628773524 , -0.8951632913550628 , 0.3612416661871546 , 0.9618256431728197 , -0.18374951781656723 , -0.9957341762950346 , 0.9957341762950346 , -0.18374951781657015 , -0.961825643172819 , 0.3612416661871526 , 0.8951632913550627 , -0.5264321628773557 , -0.7980172272802393 , 0.6736956436465567 , 0.6736956436465588 , -0.7980172272802386 , -0.5264321628773566 , 0.8951632913550622 , 0.36124166618715275 , -0.9618256431728193 , -0.1837495178165725 , 0.9957341762950344 , 0.961825643172819 , -0.5264321628773558 , -0.6736956436465578 , 0.8951632913550623 , 0.18374951781657017 , -0.9957341762950347 , 0.3612416661871515 , 0.7980172272802394 , -0.7980172272802386 , -0.36124166618715264 , 0.9957341762950344 , -0.18374951781656723 , -0.8951632913550613 , 0.6736956436465549 , 0.526432162877357 , -0.9618256431728192 , 0.8951632913550626 , -0.7980172272802388 , -0.18374951781657006 , 0.9618256431728196 , -0.673695643646558 , -0.3612416661871524 , 0.9957341762950345 , -0.5264321628773524 , -0.5264321628773566 , 0.9957341762950347 , -0.3612416661871512 , -0.6736956436465564 , 0.9618256431728203 , -0.18374951781657048 , -0.7980172272802418 , 0.8951632913550588 , 0.7980172272802396 , -0.961825643172819 , 0.3612416661871526 , 0.5264321628773563 , -0.9957341762950347 , 0.6736956436465567 , 0.18374951781657042 , -0.8951632913550628 , 0.8951632913550622, -0.18374951781656723 , -0.673695643646559 , 0.9957341762950344 , -0.5264321628773581 , -0.3612416661871531 , 0.9618256431728218 , -0.7980172272802382 , 0.6736956436465571 , -0.9957341762950345 , 0.7980172272802394 , -0.1837495178165712 , -0.5264321628773564 , 0.9618256431728192 , -0.8951632913550622 , 0.3612416661871546 , 0.36124166618715275 , -0.8951632913550629 , 0.9618256431728183 , -0.5264321628773551 , -0.18374951781656926 , 0.7980172272802398 , -0.9957341762950344 , 0.6736956436465599, 0.5264321628773561 , -0.8951632913550626 , 0.9957341762950345 , -0.7980172272802387 , 0.3612416661871515 , 0.18374951781657042 , -0.6736956436465589 , 0.9618256431728197 , -0.9618256431728193 , 0.6736956436465549 , -0.18374951781656698 , -0.3612416661871531 , 0.7980172272802355 , -0.9957341762950349 , 0.8951632913550587 , -0.5264321628773516 , 0.3612416661871533 , -0.6736956436465578 , 0.8951632913550627 , -0.9957341762950347 , 0.9618256431728189 , -0.7980172272802386 , 0.5264321628773553 ,-0.18374951781656723 , -0.1837495178165725 , 0.526432162877357 , -0.7980172272802418 , 0.9618256431728199 , -0.995734176295035 , 0.8951632913550618 , -0.6736956436465571 , 0.36124166618714704 , 0.18374951781657037 , -0.361241666187153 , 0.5264321628773563 , -0.6736956436465573 , 0.7980172272802394 , -0.8951632913550628 , 0.9618256431728197 , -0.9957341762950346 , 0.9957341762950344 , -0.9618256431728192 , 0.8951632913550588 , -0.7980172272802382 , 0.6736956436465599 , -0.5264321628773516 , 0.36124166618714704 , -0.18374951781656976};

#endif
#endif
