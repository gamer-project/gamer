#!/bin/bash

__gamer_check_gamer_info() {
    # Check if the current directory is GAMER and the `configure.py` for generating `Makefile`
    # $1 : configure.py filename

    local isGitDir=`git rev-parse --is-inside-work-tree 2>/dev/null`
    [ ! "$isGitDir" = true ] && return 1

    local firstGAMERCommit=`git rev-list --max-parents=0 HEAD`
    [ ! "$firstGAMERCommit" = "2705887745ac1fd254774c512af4a05b0f7f7fb6" ] && return 1

    local configurePyFirstCommit=`git log --diff-filter=A -- $1 | head -n 1 | awk '{print $2}'`
    [ ! "$configurePyFirstCommit" = "a93978a39b4388d169a0f45d1987f55486578fab" ] && return 1

    return 0

} # __gamer_check_gamer_info()


__gamer_configure_autocomplete() {

    local all_options sub_options all_option_array sub_option_array configure_filename
    local not_set=true
    local first=${COMP_WORDS[0]}
    local second=${COMP_WORDS[1]}
    local subsub="${COMP_WORDS[COMP_CWORD-2]}"
    local sub="${COMP_WORDS[COMP_CWORD-1]}"
    local cur="${COMP_WORDS[COMP_CWORD]}"

    if [[ "$first" == "python"* ]]; then
        configure_filename=$second
    else
        configure_filename=$first
    fi

    __gamer_check_gamer_info $configure_filename

    # Not the GAEMR `configure.py`
    if [ $? -ne 0 ]; then
        return 0
    fi

    all_options=$(./configure.py --autocomplete_info=all)
    IFS=' ' read -r -a all_option_array <<< "${all_options}"

    COMPREPLY=() # NOTE: please add a space when ending the option

    # option choices
    for opt in "${all_option_array[@]}"
    do
        # --option=xx
        if [[ "$opt" == "$subsub=" && "=" == "$sub" ]]; then
            sub_options=$(./configure.py --autocomplete_info="$opt")
            IFS=' ' read -r -a sub_option_array <<< "${sub_options}"
            for opt2 in "${sub_option_array[@]}"
            do
                if [[ "$opt2" != "$cur"* ]]; then continue; fi
                COMPREPLY+=("$opt2 ")
            done
            not_set=false
            break
        # --option, --option xxx, or --option=
        elif [[ "$opt" == "$sub=" ]]; then
            sub_options=$(./configure.py --autocomplete_info="$opt")
            IFS=' ' read -r -a sub_option_array <<< "${sub_options}"
            for opt2 in "${sub_option_array[@]}"
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

complete -o default -o nospace -F __gamer_configure_autocomplete ./configure.py
complete -o default -o nospace -F __gamer_configure_autocomplete python
complete -o default -o nospace -F __gamer_configure_autocomplete python3
