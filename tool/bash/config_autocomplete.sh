#!/bin/bash

_gamer_configure_autocomplete() {
    local options
    local subcommand="${COMP_WORDS[COMP_CWORD-1]}"
    all_options=$(./configure.py --autocomplete_info=all)
    IFS=' ' read -r -a option_array <<< "${all_options}"

    not_set=true
    for opt in "${option_array[@]}"
    do
        if [[ "$subcommand" != "$opt" ]]; then continue; fi

        options=$(./configure.py --autocomplete_info="$opt")
        not_set=false
        break
    done

    if ${not_set}; then options=${all_options}; fi

    COMPREPLY=( $(compgen -W "${options}" -- "${COMP_WORDS[COMP_CWORD]}") )
}

complete -F _gamer_configure_autocomplete ./configure.py
