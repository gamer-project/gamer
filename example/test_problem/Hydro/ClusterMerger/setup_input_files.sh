#!/bin/bash
set -e

if [ "$#" -gt 1 ]; then
    echo "Error: Too many arguments."
    echo "Usage: $0 [1-9]"
    exit 1
fi


if [ "$#" -eq 1 ]; then
    if [[ ! "$1" =~ ^[1-9]$ ]]; then
        echo "Error: Argument must be a single digit from 1 to 9."
        echo "Usage: $0 [1-9]"
        exit 1
    fi

    printf "Linking default case ${1} input files ..."
    ln -fs "./Input_Options/case${1}/Input__Flag_Lohner"
    ln -fs "./Input_Options/case${1}/Input__Flag_NParPatch"
    ln -fs "./Input_Options/case${1}/Input__Parameter"
    ln -fs "./Input_Options/case${1}/Input__TestProb"
    printf " Done!\n"
else
    printf "Linking default case input files ..."
    ln -fs "./Input_Options/default/Input__Flag_Lohner"
    ln -fs "./Input_Options/default/Input__Flag_NParPatch"
    ln -fs "./Input_Options/default/Input__Parameter"
    ln -fs "./Input_Options/default/Input__TestProb"
    printf " Done!\n"
fi
