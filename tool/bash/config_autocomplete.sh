#!/bin/bash

__gamer_configure_autocomplete() {

    local options option_array option_array2
    local not_set=true
    local subsub="${COMP_WORDS[COMP_CWORD-2]}"
    local sub="${COMP_WORDS[COMP_CWORD-1]}"
    local cur="${COMP_WORDS[COMP_CWORD]}"
    local all_options=$(./configure.py --autocomplete_info=all)
    IFS=' ' read -r -a option_array <<< "${all_options}"

    COMPREPLY=() # NOTE: please add a space when ending the option

    # option choices
    for opt in "${option_array[@]}"
    do
        # --option=xx
        if [[ "$opt" == "$subsub=" && "=" == "$sub" ]]; then
            options=$(./configure.py --autocomplete_info="$opt")
            IFS=' ' read -r -a option_array2 <<< "${options}"
            for opt2 in "${option_array2[@]}"
            do
                if [[ "$opt2" != "$cur"* ]]; then continue;fi
                COMPREPLY+=("$opt2 ")
            done
            not_set=false
            break
        # --option, --option xxx, or --option=
        elif [[ "$opt" == "$sub=" ]]; then
            options=$(./configure.py --autocomplete_info="$opt")
            IFS=' ' read -r -a option_array2 <<< "${options}"
            for opt2 in "${option_array2[@]}"
            do
                # --option= or --option
                if [[ "$cur" == "" || "$cur" == "=" ]]; then
                    COMPREPLY+=("$opt2 ")
                # --option xxx
                elif [[ "$opt2" == "$cur"* ]]; then
                    COMPREPLY+=("$opt2 ")
                fi
            done
            not_set=false
            break
        # special case without trailing `=` (-h and -lh)
        elif [[ "$cur" == "-h" || "$cur" == "-lh" ]]; then
            COMPREPLY+=("$cur ")
            not_set=false
            break
        fi
    done

    if ${not_set}; then
        COMPREPLY=( $(compgen -W "${all_options}" -- "$cur") )
    fi

} # __gamer_configure_autocomplete()

complete -o nospace -F __gamer_configure_autocomplete ./configure.py
