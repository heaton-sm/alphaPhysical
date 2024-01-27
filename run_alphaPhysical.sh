#!/bin/bash

# This script takes protein sequence(s) in #
# FASTA format and applies them to the     #
# Alphafold DNN structure inference        #
# software. It then extracts various 	   #
# biophysical information from the         #
# resulting best/relaxed structure and     #
# generates various plots.	               #

# Steven Heaton
# 2024-01-18

# define version info
    name="alphaPhysical"
    version="v0.1"
    build_date="2024/01/27"

# define default args
    alphafold_path="/mnt/d/analyze/alphafold/run"
    alphafold_db="/mnt/e/alphafold/database"
    fasta_in="/mnt/d/analyze/alphafold/target/active.done"
    output="/mnt/d/analyze/alphafold/target"

# define colours
    NC=$'\e[0;0m' RED=$'\e[0;31m' PURPLE=$'\e[0;35m'

# check for Docker
	if ! docker info > /dev/null 2>&1; then
		echo "${RED}[ERROR] This software uses Docker, but Docker it isn't running."
		echo "        Please start Docker and try again.${NC}"
		exit 1
	fi

# set pipeline to run
	pipeline_run="default"

# program suite information
Help(){
    cat <<HELP

Program: "${name}"
Version: "${version}"
Build:   "${build_date}"

alphaPhysical takes protein sequences in FASTA format and iterates these through Alphafold,
then extracts and plots various biophysical properties from the resulting best/relaxed structure.

Usage: run_alphaPhysical --[required arguments] (value) --[optional arguments] (value)
       Arguments may be all defined by either whitespace ' ' or equalities '='

       Alternatively, set required arguments under 'define default args' in run_alphaPhysical.sh.

Required arguments:
    -a, --alphafold_path [str]     Path to Alphafold installation.
	                                 (e.g., "${alphafold_path}")
    -d, --alphafold_db [str]       Path to Alphafold protein/genetic databases.
	                                 (e.g., "${alphafold_db}")
    -i, --fasta_in [str]           Path to protein sequences in FASTA format.
	                                 (e.g., "${fasta_in}")
    -o, --output [str]             Directory to write pipeline outputs
	                                 (e.g., "${output}")

Optional arguments:
    -h, --help                     Display this help information and exit.
    -v, --version                  Display version information and exit.
    -z, --license                  Display licence information and exit.

HELP
}

# interpret user args
while [[ $# -gt 0 ]]; do
    key="$1"
    case "${key}" in
        --alphafold_path|-a|--alphafold_db|-d|--fasta_in|-i|--output|-o)
            value="${2#*[= ]}"
            
            if [[ -z "${value}" || "${value}" == -* ]]; then
                Help
                echo -e "${RED}[ERROR] Missing or invalid value for option '${key}'${NC}"
                exit 1
            fi
            
            case "${key}" in
                --alphafold_path|-a) alphafold_path="${value}";;
                --alphafold_db|-d) alphafold_db="${value}";;
                --fasta_in|-i) fasta_in="${value}";;
                --output|-o) output="${value}";;
            esac
            
            shift 2;;
        
        --version|-v)
            echo "${name} ${version} ${build_date}"; echo; exit 0 ;;
        --license|-z)
            cat "${PWD}/LICENSE.md"; echo; exit 0 ;;
        --help|-h)
            Help; exit 0 ;;
        -*|--*)
            Help; echo -e "${RED}[ERROR] Unknown option '${key}'${NC}"; exit 1 ;;
    esac
done

# define preflight args
	alphafold_runfile="$( find "${alphafold_path}" -type f -name "run_docker.py" )"
	files="$( find "${fasta_in}" -type f -name *.fasta )"

# display args
	echo -e "${PURPLE}[INFO] Arguments interpreted:"
	for var in alphafold_path alphafold_db fasta_in output; do
		echo -e "           ${var}=${!var}"
	done
	echo -e "${NC}"


## DEFINE PIPELINE INSTRUCTIONS ##
	# function for checking for valid Alphafold install path
		check_alphafold_install(){
			echo ${alphafold_runfile}
			if [[ ! -f "${alphafold_runfile}" ]]; then
				echo -e "${RED}[ERROR] Alphafold installation not detected. Please check installation path.${NC}"
				exit 1
			else
				echo -e "${PURPLE}[INFO] Alphafold installation detected.${NC}"
			fi
		}

	# function for checking for valid Alphafold database
		check_database_install(){
			if [[ ! -d "${alphafold_db}" ]]; then
				echo -e "${RED}[ERROR] Alphafold database not detected. Please check database path.${NC}"
				exit 1
			else
				echo -e "${PURPLE}[INFO] Alphafold database detected.${NC}"
			fi
		}

	# function to execute Alphafold on input FASTA files
		run_alphafold(){
			activePDB="$( basename -- "${file}" | cut -d '.' -f 1 )"
			if [[ -d "${output}/complete/${activePDB}/msas" ]]; then
				echo -e "${PURPLE}[INFO] ${activePDB} structure already inferred. Skipping.${NC}"
			else
				echo -e "${PURPLE}[INFO] Running Alphafold DNN inference on ${activePDB}...${NC}"
				[[ ! -d "${output}/complete/${activePDB}" ]] && mkdir -p "${output}/complete/${activePDB}"
				sudo python3 "${alphafold_runfile}" \
					--fasta_paths="${file}" \
					--output_dir="${output}/complete/" \
					--max_template_date="2022-12-17" \
					--model_preset="monomer" \
					--data_dir="${alphafold_db}" \
					&& echo "${file}" >> "${output}/inferred_proteins.txt"
			fi
		}

	# function for extracting biophysical parameters from inferred structures
		get_structure_parameters(){
			echo -e "${PURPLE}[INFO] Extracting biophysical properties of inferred structures...${NC}"
			docker run -it --rm --name "alphaphysical_extract" \
				-v /var/run/docker.sock:/var/run/docker.sock \
				-v /usr/bin/com.docker.cli:/usr/bin/docker \
				-v "${output}":/project/target \
					alphaphysical:latest python3 /project/scripts/get_ss_charge.py
		}

	# function for plotting extracted biophysical parameters
		plot_structure_parameters(){
			echo -e "${PURPLE}[INFO] Creating structure plots...${NC}"
			docker run -it --rm --name "alphaphysical_plot" \
				-v /var/run/docker.sock:/var/run/docker.sock \
				-v /usr/bin/com.docker.cli:/usr/bin/docker \
				-v "${output}":/project/target \
					alphaphysical:latest python3 /project/scripts/plot_ss_results.py
		}

	# function for extracting charge states
		extract_charge_values(){
			echo -e "${PURPLE}[INFO] Extracting charge states...${NC}"
			docker run -it --rm --name "alphaphysical_chargeextract" \
				-v /var/run/docker.sock:/var/run/docker.sock \
				-v /usr/bin/com.docker.cli:/usr/bin/docker \
				-v "${output}":/project/target \
					alphaphysical:latest bash /project/scripts/extract_charge_values.sh
		}

	# function for plotting surface charge maps
		create_charge_maps(){
			echo -e "${PURPLE}[INFO] Creating charge maps...${NC}"
			docker run -it --rm --name "alphaphysical_chargemap" \
				-v /var/run/docker.sock:/var/run/docker.sock \
				-v /usr/bin/com.docker.cli:/usr/bin/docker \
				-v "${output}":/project/target \
					alphaphysical:latest python3 /project/scripts/create_charge_map.py
		}

	# function for plotting extracted charge states
		plot_charge_curves(){
			echo -e "${PURPLE}[INFO] Creating charge plots...${NC}"
			docker run -it --rm --name "alphaphysical_chargeplot" \
				-v /var/run/docker.sock:/var/run/docker.sock \
				-v /usr/bin/com.docker.cli:/usr/bin/docker \
				-v "${output}":/project/target \
					alphaphysical:latest python3 /project/scripts/plot_charge_curves.py
		}

	# function for performing structural alignments and calculating RMSD
		align_all_structures(){
			echo -e "${PURPLE}[INFO] Performing structural alignment...${NC}"
			docker run -it --rm --name "alphaphysical_chargeplot" \
				-v /var/run/docker.sock:/var/run/docker.sock \
				-v /usr/bin/com.docker.cli:/usr/bin/docker \
				-v "${output}":/project/target \
					alphaphysical:latest python3 /project/scripts/align_all_structures.py
		}

## EXECUTION PIPELINES ##
	# Default execution pipeline for Alphafold inference
	# followed by biophysical analysis and plotting.
	# Commented modules are currently not functional.
		if [[ "${pipeline_run}" == "default" ]]; then
			check_alphafold_install
			check_database_install
			for file in ${files}; do
				run_alphafold
			done

			get_structure_parameters
			extract_charge_values
			plot_charge_curves
			# create_charge_maps
			plot_structure_parameters
			# align_all_structures

			echo -e "${PURPLE}Analysis finished."
		fi
