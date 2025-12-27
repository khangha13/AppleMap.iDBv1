#!/bin/bash

prompt_dataset_name() {
    local dataset
    while true; do
        read -p "Dataset name: " dataset
        if [ -n "${dataset}" ]; then
            echo "${dataset}"
            return 0
        fi
        echo "❌ Dataset name cannot be empty. Please try again."
    done
}

prompt_directory() {
    local prompt_msg="$1"
    local default_path="$2"
    local path
    while true; do
        read -p "${prompt_msg} [${default_path}]: " path
        path=${path:-$default_path}
        if [ -d "${path}" ]; then
            echo "${path}"
            return 0
        fi
        echo "❌ Path does not exist: ${path}"
    done
}

confirm_action() {
    local message="$1"
    local response
    while true; do
        read -p "${message} [y/n]: " response
        case "${response}" in
            [Yy]*) return 0 ;;
            [Nn]*) return 1 ;;
            *) echo "Please answer y or n." ;;
        esac
    done
}
